
## this script is derived in part from the RforMassSpectrometry 'Metabonoaut' efforts
## modified to make it compatible with LC or GC Q-TOF data using Waters MSe based DIA MS/MS.  

### install and load all packages
install.packages(c('readxl', 'BiocManager', 'remotes', 'pander', 'pheatmap', 'vioplot', 'ggfortify', 'gridExtra', 'ggVennDiagram', 'UpSetR', 'fastmatch'), update = FALSE)
BiocManager::install(c('xcms', "CAMERA", "BiocStyle", "alabaster.se",
                       "RforMassSpectrometry/RforMassSpectrometry",
                       "RforMassSpectrometry/MsIO", "msdata",
                       "RforMassSpectrometry/MsBackendMetaboLights",
                       "RforMassSpectrometry/MsBackendMgf",
                       "AnnotationHub", "CompoundDb", "MetaboAnnotation"),
                     update = FALSE)

remotes::install_github(c('cbroeckl/RAMClustR', 'cbroeckl/csu.pmf.tools', 'cbroeckl/csu.arcbio.tools'))

## load libraries
## Data Import and handling
library(readxl)
library(MsExperiment)
library(MsIO)
# library(alabaster.se)
# library(MsBackendMetaboLights)
library(SummarizedExperiment)
library(xcms)
library(Spectra)
library(MetaboCoreUtils)
library(limma) # Differential abundance
library(matrixStats) # Summaries over matrices
library(pander)
library(RColorBrewer)
library(pheatmap)
library(vioplot)
library(ggfortify)   # Plot PCA
library(gridExtra)   # To arrange multiple ggplots into single plots
library(knitr)
library(AnnotationHub) # Annotation resources
library(CompoundDb)    # Access small compound annotation data.
library(MetaboAnnotation) # Functionality for metabolite annotation.


#####
# VARIABLES
#####

## how many parallel processes.  Budget about 6-8 GB RAM per process
## i.e. if you have 48 threads available, but only 64 GB RAM, 
## do not use n.threads <- 48, but rather n.threads <- 8

n.threads <- 2


## set up parallel processing: 
## for linux/unix, use MulticoreParam, and for everything else SnowParam
if (.Platform$OS.type == "unix") {  
  register(MulticoreParam(n.threads))
} else {
  register(SnowParam(n.threads))
}

## files to process
## your mzML files and a file called sequence.csv should be in the 'in.dir' location.  
## output files will be written to the 'out.dir' location. 
in.dir  <- 'R:/RSTOR-PMF/Projects/20250628-MH-6116/qTof/G2-XS_Data/mzML/'  #R:\RSTOR-PMF\Projects\20250628-MH-6116\qTof
out.dir <- 'R:/RSTOR-PMF/Projects/20250628-MH-6116/qTof/G2-XS_Data/'

# 
# files <- list.files(in.dir, pattern = ".mzML", full.names = TRUE)
# files <- sort(files[grep("_Org_", files)])
# basename.files <- basename(files)
# cat(basename.files, sep = '\n')

## read in sequence with sample metadata/phenotype
sd <- DataFrame(readxl::read_xlsx(paste0(in.dir, 'sequence.xlsx')))
sd <- sd[order(sd$run.order),]
filenames <- paste0(sd[,'file'], ".mzML", sep = '')
filepaths <- paste0(in.dir, "/", filenames, sep = '')
spec <- Spectra::Spectra(filepaths, backend = MsBackendMzR())

## set up MsExperiment object called 'mse'
msefiles <- MsExperimentFiles(mzML_files = paste0(sd[,'file'], ".mzML", sep = ''))
mse <- MsExperiment()
experimentFiles(mse) <- msefiles
mse <- linkSampleData(mse, with = "experimentFiles.mzML_files", sampleIndex = 1:length(sd$file), withIndex = 1:length(sd$file))
spectra(mse) <- spec
sampleData(mse) <- sd
sampleData(mse)$raw_file <- normalizePath(filepaths)
mse <- linkSampleData(mse, with = "sampleData.raw_file = spectra.dataOrigin")



#####################################
### this section does some basic data visualization
### it does not represent any high level processing
### but can be useful for some simple visual QC


## set sample.type colors
sample.types <- unique(sampleData(mse)[,"sample.type"])
col_phenotype <- brewer.pal(length(sample.types), name = "Set1")
names(col_phenotype) <- sample.types
col_sample <- col_phenotype[sampleData(mse)[,"sample.type"]]


## get base peak ion chromatograms
bpc <- chromatogram(mse, aggregationFun = "max")
plot(bpc, col = paste0(col_sample, 80), main = "BPC", lwd = 1.5, xlim = c(30, ))
grid()
legend("topright", col = col_phenotype,
       legend = names(col_phenotype), lty = 1, lwd = 2, horiz = TRUE,
       bty = "n")

#' Filter the data based on retention time
# mse2 <- filterRt(mse, c(200, 240))
# bpc <- chromatogram(mse2, aggregationFun = "max")
# plot(bpc, col = paste0(col_sample, 80), main = "BPC", lwd = 1.5)
# grid()
# legend("topright", col = col_phenotype,
#        legend = names(col_phenotype), lty = 1, lwd = 2, horiz = TRUE,
#        bty = "n")

## print sampleData table
# sampleData(mse) |>
#   kable(format = "pipe")


## BPC heatmap
bpc <- bpc |>
  bin(binSize = 2)

#' Calculate similarity (Pearson correlation) between BPCs
bpcmap <- do.call(cbind, lapply(bpc, intensity)) |>
  cor()

rownames(bpcmap) <- colnames(bpcmap) <- paste0(sampleData(mse)$sample, ".", sampleData(mse)$rep, sampleData(mse)$msms.type)
ann <- data.frame(phenotype = sampleData(mse)$sample.type)
rownames(ann) <- rownames(bpcmap)

#' Plot heatmap
pheatmap(bpcmap, annotation_col = ann,
         annotation_colors = list(phenotype = col_phenotype))

### end simple visual QC
#####################################


#####################################
### now we do the computationally demanding parts, processing the dat with XCMS
### the parameters below are quite reasonable for our LC-Q-TOF systems
### that doesn't mean they are optimal, and determining optimal is quite difficult
### i suggest starting here, and if we see any anomolies along the way, lets discuss


#####################################
### XCMS


## detect peaks with CentWave
param <- CentWaveParam(peakwidth = c(1.5, 12), ppm = 20, integrate = 1, prefilter = c(3,200),
                       extendLengthMSW = TRUE, snthresh = 10,verboseBetaColumns = TRUE)

xd <- findChromPeaks(mse, param = param, chunkSize = 2)

## detect DIA (MSe) MS/MS peaks with CentWave 
param.2 <- CentWaveParam(peakwidth = c(1.5, 12), ppm = 20, integrate = 1, prefilter = c(3,80),
                         extendLengthMSW = TRUE, snthresh = 10, verboseBetaColumns = TRUE)
xd <- findChromPeaksIsolationWindow(xd, param = param.2)

## remove lower quality peaks via the beta_cor metric.  
cp <- chromPeaks(xd)
# ggplot(cp, aes(x = beta_cor, y = beta_snr)) +
#   geom_point() + 
#   geom_density_2d()
xd <- filterChromPeaks(xd, keep = cp[,"beta_cor"] > 0.8)

## refine peaks, to remove some peak boundry artifacts
param <- MergeNeighboringPeaksParam(expandRt = 2.5, expandMz = 0.0015,
                                    minProp = 0.75)
xd <- refineChromPeaks(xd, param = param, chunkSize = n.threads, msLevel = c(1L))
xd <- refineChromPeaks(xd, param = param, chunkSize = n.threads, msLevel = c(2L))

## group peaks with PeakDensity
## minFraction is calcluated by sample.type.  you can change this to any other variable you wish
## just keep in mind that you want it to be a variable that is replicated
param <- PeakDensityParam(sampleGroups = sampleData(xd)$sample.type,
                          minFraction = 0.4,
                          binSize = 0.02, ppm = 10,
                          bw = 2)
xd <- groupChromPeaks(xd, param = param, msLevel = 1L)

## Rt alignment
param <- PeakGroupsParam(minFraction = 0.4, extraPeaks = 50, span = 0.5,
                         subsetAdjust = "average")
xd <- adjustRtime(xd, param = param, msLevel = c(1L))

## regroup post retention time correction
param <- PeakDensityParam(sampleGroups = sampleData(xd)$sample.type,
                          minFraction = 0.1, binSize = 0.01, ppm = 10,
                          bw = 1.8)
xd <- groupChromPeaks(xd, param = param, msLevel = c(1L))
xd <- groupChromPeaks(xd, param = param, msLevel = c(2L), add = TRUE)


# fill peaks
xd <- fillChromPeaks(xd, param = ChromPeakAreaParam(), chunkSize = n.threads)
xd <- fillChromPeaks(xd, param = ChromPeakAreaParam(), chunkSize = n.threads, msLevel = c(2L))

save(xd, file = paste0(out.dir, "/xcms.object.Rdata"))  ## load(paste0(out.dir, "/xcms.object.Rdata"))

### end XCMS
#####################################



#####################################
### RAMClustR
### RAMClustR package takes XCMS feature output, and clusters features into compounds/spectra
### to reduce the complexity of the dataset and increase the interpretability.  

### bring XCMS data into ramclustR ready format
RC <- rc.get.xcms.data2(xcmsObj = xd)

### cluster features
RC <- RAMClustR::rc.ramclustr(ramclustObj = RC)

### annotate ions in each cluster
RC <- doFindmain2(ramclustObj = RC)

### assign fragment ions from MSe DIA data
RC <- assign.xcms.mse.peaks(ramclustObj = RC)

### export spectra for annotation using Sirius
export.ms.formatted.spectra(ramclustObj = RC, output.dir = out.dir)

save(RC, file = paste0(out.dir, "/ramclustr.object.Rdata"))  ## load(paste0(out.dir, "/ramclustr.object.Rdata"))

### end RAMClustR
#####################################

