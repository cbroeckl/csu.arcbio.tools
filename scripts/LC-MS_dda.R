
## this script is derived in part from the RforMassSpectrometry 'Metabonoaut' efforts
## modified to make it compatible with LC or GC Q-TOF data using Waters MSe based DIA MS/MS.  

### install and load all packages
# install.packages(c('readxl', 'BiocManager', 'remotes', 'pander', 'pheatmap', 'vioplot', 'ggfortify', 'gridExtra', 'ggVennDiagram', 'UpSetR', 'fastmatch'), update = FALSE)
# BiocManager::install(c('xcms', "CAMERA", "BiocStyle", "alabaster.se", "MAI",
#                        "RforMassSpectrometry/RforMassSpectrometry",
#                        "RforMassSpectrometry/MsIO", "msdata",
#                        "RforMassSpectrometry/MsBackendMetaboLights",
#                        "RforMassSpectrometry/MsBackendMgf",
#                        "AnnotationHub", "CompoundDb", "MetaboAnnotation"),
#                      update = FALSE)
# 
# remotes::install_github(c('cbroeckl/RAMClustR', 'cbroeckl/csu.pmf.tools', 'cbroeckl/csu.arcbio.tools'))

## load libraries
## Data Import and handling
library(readxl)
library(MsExperiment)
# library(MsIO)
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
library(plotly)

#####
# VARIABLES
#####

## how many parallel processes.  Budget about 6-8 GB RAM per process
## i.e. if you have 48 threads available, but only 64 GB RAM, 
## do not use n.threads <- 48, but rather n.threads <- 8 or fewer
## more threads consume more RAM.  More threads will run faster
## than fewer until you run out of RAM, then it slows down terribly.
n.threads <- 4

## files to process
## your mzML files and a file called sequence.csv should be in the 'in.dir' location.  
## output files will be written to the 'out.dir' location. 
in.dir  <- 'R:/RSTOR-PMF/Projects/20250606-JF-6046/lc-tof/mzml2/'  #R:\RSTOR-PMF\Projects\20250628-MH-6116\qTof
out.dir <- 'R:/RSTOR-PMF/Projects/20250606-JF-6046/lc-tof/mzml2/out/'


color.bpc.by <- "sample.type"
filename.column.name <- "filename"




## set up parallel processing: 
## for linux/unix, use MulticoreParam, and for everything else SnowParam
if (.Platform$OS.type == "unix") {  
  register(MulticoreParam(n.threads))
} else {
  register(SnowParam(n.threads))
}

if(!dir.exists(out.dir <- 'R:/RSTOR-PMF/Projects/20250606-JF-6046/lc-tof/mzml2/out/')) {
  dir.create(out.dir <- 'R:/RSTOR-PMF/Projects/20250606-JF-6046/lc-tof/mzml2/out/')
}



# 
# files <- list.files(in.dir, pattern = ".mzML", full.names = TRUE)
# files <- sort(files[grep("_Org_", files)])
# basename.files <- basename(files)
# cat(basename.files, sep = '\n')

## read in sequence with sample metadata/phenotype
sd <- DataFrame(readxl::read_xlsx(paste0(in.dir, 'sequence.xlsx')))
sd <- sd[order(sd$run.order),]
filenames <- paste0(sd[,'filename'], ".mzML", sep = '')
filepaths <- paste0(in.dir, "/", filenames, sep = '')
spec <- Spectra::Spectra(filepaths, backend = MsBackendMzR())

## set up MsExperiment object called 'mse'
msefiles <- MsExperimentFiles(mzML_files = paste0(sd[,filename.column.name], ".mzML", sep = ''))
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
sample.types <- unique(sampleData(mse)[,color.bpc.by])
col_phenotype <- brewer.pal(length(sample.types), name = "Set1")
names(col_phenotype) <- sample.types
col_sample <- col_phenotype[sampleData(mse)[,color.bpc.by]]


## get base peak ion chromatograms
bpc <- chromatogram(mse, aggregationFun = "max")
save(bpc, file = paste0(out.dir, "/xcms.bpc.Rdata"))  ## load(paste0(out.dir, "/xcms.object.Rdata"))


plot(bpc, col = paste0(col_sample, 80), main = "BPC", lwd = 1.5, xlim = c(30, 750))
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

rownames(bpcmap) <- colnames(bpcmap) <- paste0(sampleData(mse)$sample.type, ".", sampleData(mse)$sample.id)
ann <- data.frame(phenotype = sampleData(mse)$sample.type)
rownames(ann) <- rownames(bpcmap)

#' Plot heatmap
pheatmap(bpcmap, annotation_col = ann,
         annotation_colors = list(phenotype = col_phenotype))

## looks like QC1, sample.46, and sample.12 are outliers
rm <- c(which(sampleData(mse)$sample.id == "QC1"),
        which(sampleData(mse)$sample.id == "46"),
        which(sampleData(mse)$sample.id == "12")
)

# Subset the MsExperiment object
mse <- mse[-rm]



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
param <- CentWaveParam(peakwidth = c(1.5, 12), ppm = 8, integrate = 1, prefilter = c(3,100),
                       extendLengthMSW = TRUE, snthresh = 3,verboseBetaColumns = TRUE)

xd <- findChromPeaks(mse, param = param, chunkSize = 2)
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

## group peaks with PeakDensity
## minFraction is calcluated by sample.type.  you can change this to any other variable you wish
## just keep in mind that you want it to be a variable that is replicated

## pick up here.  determine best column(s) to group by
head(sampleData(xd))
param <- PeakDensityParam(sampleGroups = sampleData(xd)$treatment,
                          minFraction = 0.4,
                          binSize = 0.02, ppm = 8,
                          bw = 2)
xd <- groupChromPeaks(xd, param = param, msLevel = 1L)

## Rt alignment
param <- PeakGroupsParam(minFraction = 0.4, extraPeaks = 50, span = 0.5,
                         subsetAdjust = "average")
xd <- adjustRtime(xd, param = param, msLevel = c(1L))

## regroup post retention time correction
param <- PeakDensityParam(sampleGroups = sampleData(xd)$sample.type,
                          minFraction = 0.6, binSize = 0.01, ppm = 8,
                          bw = 1.5)
xd <- groupChromPeaks(xd, param = param, msLevel = c(1L))


# fill peaks
xd <- fillChromPeaks(xd, param = ChromPeakAreaParam(), chunkSize = n.threads)

save(xd, file = paste0(out.dir, "/xcms.object.Rdata"))  ## load(paste0(out.dir, "/xcms.object.Rdata"))


xd.msms <- xcms::featureSpectra(xd, ppm = 100, skipFilled = TRUE, return.type = "Spectra")
save(xd.msms, file = paste0(out.dir, "/xcms.msms.object.Rdata"))  ## load(paste0(out.dir, "/xcms.object.Rdata"))

### end XCMS
#####################################



#####################################
### RAMClustR
### RAMClustR package takes XCMS feature output, and clusters features into compounds/spectra
### to reduce the complexity of the dataset and increase the interpretability.  

### bring XCMS data into ramclustR ready format
RC <- csu.arcbio.tools::rc.get.xcms.data(xcmsObj = xd, fill.na = FALSE)

### impute missing values using 
RC <- RAMClustR::rc.feature.replace.na.mai(ramclustObj = RC)

RC <- RAMClustR::rc.feature.filter.blanks(ramclustObj = RC, blank.tag = c("blank", "sample.type"),
                                          qc.tag = c("qc", "sample.type"), sn = 5)

### normalize to QC samples
RC <- RAMClustR::rc.feature.normalize.qc(ramclustObj = RC, order = RC$phenoData$run.order, 
                                         qc = RC$phenoData$sample.type == "qc")



### cluster features
RC <- RAMClustR::rc.ramclustr(ramclustObj = RC)

### annotate ions in each cluster
RC <- csu.arcbio.tools::doFindmain(ramclustObj = RC, mode = "positive")

### assign fragment ions from DDA data
RC <- csu.arcbio.tools::assign.xcms.dda.peaks(ramclustObj = RC, xcms.feature.msmsObj = xd.msms)

### export spectra for annotation using Sirius
csu.arcbio.tools::export.ms.formatted.spectra(ramclustObj = RC, output.dir = out.dir)

save(RC, file = paste0(out.dir, "/ramclustr.object.Rdata"))  ## load(paste0(out.dir, "/ramclustr.object.Rdata"))

### end RAMClustR
#####################################

### QC Report
#####################################
setwd(out.dir)
RC <- RAMClustR::rc.qc(ramclustObj = RC, qc.tag = c("qc", "sample.type"), remove.qc = TRUE)

### STATS
#####################################
RC <- csu.pmf.tools::pmfpca(ramclustObj = RC, which.factors = 'treatment')
RC <- csu.pmf.tools::pmfanova(ramclustObj = RC, anova.name = "treatment", anova.call = "treatment")


### Spectra Matching
#####################################
library(Spectra)
load("R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.pos.Rdata")
## set chunksize - smaller chunks seems to process faster
nist.hr.pos@processingChunkSize <- 100

rc.spec <- as.list(rep(NA, length(RC$dda.spectra)))
low_int <- function(x, ...) {
  x > max(x, na.rm = TRUE) * 0.02
}
divide_intensities <- function(x, y, ...) {
  x[, 2] <- 1000* (x[, 2] / max(x[,2]))
  x
}

for(i in 1:length(rc.spec)) {
  # RC$dda.spectra
  len <- length(RC$dda.spectra[[i]])
  if(len == 0) next
  spd <- DataFrame(msLevel = rep(2L, len), rtime = rep(RC$clrt[i], len), 
                   name = rep(RC$cmpd[[i]], len), index = 1:len)
  spd$mz <- IRanges::NumericList(0)
  spd$intensity <- IRanges::NumericList(0)
  for(j in 1:len) {
    spd$mz[[j]] <- RC$dda.spectra[[i]][[j]]$spectrum$mz
    spd$intensity[[j]] <- RC$dda.spectra[[i]][[j]]$spectrum$int 
  }
  # spd$adduct <- sapply(1:length(RC$dda.spectra[[i]]), FUN = function(x) RC$dda.spectra[[i]][x][[1]]$adduct.type)
  rc.spec[[i]] <- Spectra(spd, backend = MsBackendDataFrame())
  rc.spec[[i]] <- filterIntensity(rc.spec[[i]], intensity = low_int)
  rc.spec[[i]] <- Spectra::addProcessing(rc.spec[[i]], divide_intensities)
  rc.spec[[i]] <- applyProcessing(rc.spec[[i]])
}

rc.spec2 <- (rc.spec[[1]])
rc.spec2 <- applyProcessing(rc.spec2)
for(i in 2:length(rc.spec)) {
  if(is.vector(rc.spec[[i]][1])) next
  rc.spec2 <- c(rc.spec2, applyProcessing(rc.spec[[i]]))
  applyProcessing(rc.spec2)
}

rc.spec2@processingChunkSize <- 100
nist.hr.pos@processingChunkSize <- 100

res <- csu.pmf.tools::lib.search(
  prod.tol.ppm = 20, prec.tol.ppm = 20, 
  exp.spectra = rc.spec2, lib.spectra = nist.hr.pos)
save(res, file = paste0(out.dir, "/dda.search.results.crude.Rdata"))

tar.cmpd <- "C2082"
rc.spec2.index <- grep(tar.cmpd, rc.spec2$name)
res.index <- which(res[,'experimental.index'] %in% rc.spec2.index)
res[res.index,]
use <- which.max(res[res.index, "spectral.similarity"])
res[res.index[use],]
nist.hr.pos[res[res.index[use], "library.index"]]
plotSpectraMirror(rc.spec2[rc.spec2.index], nist.hr.pos[res[res.index[use], "library.index"]], 
                  main = nist.hr.pos[res[res.index[use], "library.index"]]$name)

names <- rc.spec2$name[res[,"experimental.index"]]
adducts <- sapply(1:length(names), 
                  FUN = function(x) {
                    cmpd <- as.integer(as.numeric(gsub("C", "", rc.spec2[as.vector(res[x,"experimental.index"])]$name)))
                    ind <- rc.spec2[as.vector(res[x,"experimental.index"])]$index
                    RC$dda.spectra[[cmpd]][[ind]]$adduct.type
                  }
)
precursor.mz <- sapply(1:length(names), 
                       FUN = function(x) {
                         cmpd <- as.integer(as.numeric(gsub("C", "", rc.spec2[as.vector(res[x,"experimental.index"])]$name)))
                         ind <- rc.spec2[as.vector(res[x,"experimental.index"])]$index
                         RC$dda.spectra[[cmpd]][[ind]]$precursor.mz
                       }
)
precursor.rt <- sapply(1:length(names), 
                       FUN = function(x) {
                         cmpd <- as.integer(as.numeric(gsub("C", "", rc.spec2[as.vector(res[x,"experimental.index"])]$name)))
                         ind <- rc.spec2[as.vector(res[x,"experimental.index"])]$index
                         RC$dda.spectra[[cmpd]][[ind]]$retention.time
                       }
)
precursor.xcms.id <- sapply(1:length(names), 
                            FUN = function(x) {
                              cmpd <- as.integer(as.numeric(gsub("C", "", rc.spec2[as.vector(res[x,"experimental.index"])]$name)))
                              ind <- rc.spec2[as.vector(res[x,"experimental.index"])]$index
                              RC$dda.spectra[[cmpd]][[ind]]$xcms.name
                            }
)

match.name <- nist.hr.pos$name[res[,"library.index"]]
match.adduct <- nist.hr.pos$adduct[res[,"library.index"]]
match.cid <- nist.hr.pos$cid[res[,"library.index"]]
match.exactmass <- nist.hr.pos$exactmass[res[,"library.index"]]
match.formula <- nist.hr.pos$formula[res[,"library.index"]]
match.inchikey <- nist.hr.pos$inchikey[res[,"library.index"]]
match.precursorMz <- nist.hr.pos$precursorMz[res[,"library.index"]]

res2 <- data.frame(
  cmpd.name = names,
  precursor.adduct = adducts,
  precursor.mz,
  precursor.rt,
  precursor.xcms.id,
  res,
  match.name,
  match.precursorMz,
  match.adduct,
  match.cid,
  match.exactmass,
  match.formula,
  match.inchikey
)

save(res2, file = paste0(out.dir, "/dda.search.results.Rdata"))

is <- read.csv("is.list.csv", header = TRUE, check.names = FALSE)
is$feature.index <- rep(NA, nrow(is))
is$cmpd.index <- rep(NA, nrow(is))
for(i in 1:nrow(is)) {
  ft <- which(abs(1000000*((RC$fmz - is$`Isotope Peaks m/z`[i])/is$`Isotope Peaks m/z`[i])) <= 10)
  if(length(ft) == 1) {
    is$feature.index[i] <- ft
    is$cmpd.index[i] <- RC$featclus[ft]
  }
  if(length(ft) > 1) {
    is$feature.index[i] <- 10000
    is$cmpd.index[i] <- 10000
  }
}

s <- Spectra::Spectra("R:/RSTOR-PMF/Projects/20250606-JF-6046/lc-tof/mzml2/JF_Dog_Plasma_Oct_2025_redo_005.mzML")


## manual
load(paste0(out.dir, "/xcms.object.Rdata"))
load(paste0(out.dir, "/xcms.msms.object.Rdata"))
load(paste0(out.dir, "/ramclustr.object.Rdata"))
load(paste0(out.dir, "rc.spec2.Rdata"))
load(paste0(out.dir, "rc.spec.Rdata"))
load(paste0(out.dir, "dda.search.results.Rdata"))


cmpd <- 1736
plot_ly(x = RC$fmz[which(RC$featclus == cmpd)], 
        y = RC$msint[which(RC$featclus == cmpd)],
        type = "bar") %>% 
  layout(title = RC$cmpd[cmpd])  
RC$ms1.isotopes[[cmpd]]
spec <- 1
plot_ly(x = RC$dda.spectra[[cmpd]][[spec]]$spectrum$mz, 
        y = RC$dda.spectra[[cmpd]][[spec]]$spectrum$int,
        type = "bar") %>% 
  layout(title = RC$dda.spectra[[cmpd]][[spec]]$adduct.type)  
tmp <- res2[which(res2$cmpd == paste0("C", formatC(cmpd, flag = 0, width = 4))),]
tmp <- tmp[order(tmp$spectral.similarity, decreasing = TRUE),]
head(tmp)


# plot_ly(x = RC$fmz[which(RC$featclus == cmpd)],
#   y = RC$msint[which(RC$featclus == cmpd)], 
#   type = "bar", width = 0.01)
# 

great.matches <- which(
  (res2$precursor.adduct == res2$match.adduct) &
    ((res2$precursor.mz - res2$match.precursorMz) < 0.01) &
    res2$spectral.similarity >= 0.7
)

spec.matches <- res2[great.matches,]
spec.matches <- spec.matches[order(spec.matches$spectral.similarity, decreasing = TRUE),]
keep <- duplicated(spec.matches[,1:2])
spec.matches <- (spec.matches[!keep,])
unique(spec.matches$match.name)
writexl::write_xlsx(spec.matches, path = 'spectra/nist.msmslibrary.matches.xlsx')

dir.create("spectra/ms.significant")
sig <- RC$cmpd[which(RC$anova.pval_treatment[1,] < 0.05)]
ms.files <- list.files("spectra/ms")
keep <- vector(mode = 'integer')
for(i in 1:length(sig)) {
  tmp <- grep(sig[i], ms.files) 
  if(length(tmp) > 0) {
    keep <- c(keep, tmp) 
  }
}
file.copy(from = paste0(out.dir, "spectra/ms/", ms.files[keep]), to = paste0(out.dir, "spectra/ms.significant/", ms.files[keep]))










### custom stats/feature filtering
cmpds <- c(77, 238, 280, 488, 646, 686, 804, 1053, 1153, 1362, 1387, 1407, 1417, 1586, 1607, 1736, 1823, 1914, 2082, 2146, 2207)
out <- data.frame(
  "cmpd" = rep(0,0),
  "feat" = rep("A", 0),
  "feat.mz" = rep(0,0),
  "feat.rt" = rep(0,0),
  "t.test.pval" = rep(0,0),
  "fold.change" = rep(0,0),
  "chi.sq.pval" = rep(0,0),
  "CD.detected.proportion" = rep(0,0),
  "Control.detected.proportion" = rep(0,0)
)

for(i in 1:length(cmpds)) {
  cmpd <- cmpds[i]
  use.features <- RC$featnames[which(RC$featclus == cmpd)]
  use.index <- which(RC2$featnames %in% use.features)
  use.features.2 <- RC2$featnames[use.index]
  d <- data.frame(
    trt = RC2$phenoData$treatment,
    RC2$MSdata.unfilled[,use.index]
  )
  d <- d[d$trt %in% c("CD", "Control"),]
  for(j in 2:ncol(d)) {
    d.sub <- d[!is.na(d[,j]),c(1, j)]
    if(min(table(d.sub$trt)) > 1 & length(table(d.sub$trt))>1) {
      tp <- t.test(d[,j]~d[,1])$p.value
      fc <- t.test(d[,j]~d[,1])$estimate[1]/t.test(d[,j]~d[,1])$estimate[2]
    } else {
      tp <- NA
      fc <- NA
    }
    d.sub <- d[,c(1, j)]
    d.sub[is.na(d.sub[,2]),2] <- 0
    d.sub[d.sub[,2]>0 ,2] <- 1
    qp <- chisq.test(x=d.sub[,2], y=as.factor(d.sub[,1]))$p.value
    CD.detected = length(which(d.sub$trt == "CD" & d.sub[,2] == 1)) / length(which(d.sub$trt == "CD"))
    Control.detected = length(which(d.sub$trt == "Control" & d.sub[,2] == 1)) / length(which(d.sub$trt == "Control"))
    new.out.row <- data.frame(
      "cmpd" = cmpd,
      "feat" = RC2$featnames[use.index[j-1]],
      "feat.mz" = RC2$fmz[use.index[j-1]],
      "feat.rt" = RC2$frt[use.index[j-1]],
      "t.test.pval" = tp,
      "fold.change" = fc,
      "chi.sq.pval" = qp,
      "CD.detected.proportion" = CD.detected,
      "Control.detected.proportion" = Control.detected
    )
    out <- rbind(out, new.out.row)
  }
  
}

boxplot(out$t.test.pval~out$cmpd)
boxplot(out$chi.sq.pval ~out$cmpd)

row.names(out) <- NULL
writexl::write_xlsx(out, paste0(out.dir, "/feature.detection.frequency.stats.xlsx"))


## now do for all compounds
cmpds <- 1:ncol(RC$SpecAbund)
out <- data.frame(
  "cmpd" = rep(0,0),
  "feat" = rep("A", 0),
  "feat.mz" = rep(0,0),
  "feat.rt" = rep(0,0),
  "t.test.pval" = rep(0,0),
  "fold.change" = rep(0,0),
  "chi.sq.pval" = rep(0,0),
  "CD.detected.proportion" = rep(0,0),
  "Control.detected.proportion" = rep(0,0)
)

for(i in 1:length(cmpds)) {
  cmpd <- cmpds[i]
  use.features <- RC$featnames[which(RC$featclus == cmpd)]
  use.index <- which(RC2$featnames %in% use.features)
  use.features.2 <- RC2$featnames[use.index]
  d <- data.frame(
    trt = RC2$phenoData$treatment,
    RC2$MSdata.unfilled[,use.index]
  )
  d <- d[d$trt %in% c("CD", "Control"),]
  for(j in 2:ncol(d)) {
    d.sub <- d[!is.na(d[,j]),c(1, j)]
    if(min(table(d.sub$trt)) > 1 & length(table(d.sub$trt))>1) {
      tp <- t.test(d[,j]~d[,1])$p.value
      fc <- t.test(d[,j]~d[,1])$estimate[1]/t.test(d[,j]~d[,1])$estimate[2]
    } else {
      tp <- NA
      fc <- NA
    }
    d.sub <- d[,c(1, j)]
    d.sub[is.na(d.sub[,2]),2] <- 0
    d.sub[d.sub[,2]>0 ,2] <- 1
    if(length(table(d.sub[,2])) > 1) {
      qp <- chisq.test(x=d.sub[,2], y=as.factor(d.sub[,1]))$p.value
      CD.detected = length(which(d.sub$trt == "CD" & d.sub[,2] == 1)) / length(which(d.sub$trt == "CD"))
      Control.detected = length(which(d.sub$trt == "Control" & d.sub[,2] == 1)) / length(which(d.sub$trt == "Control")) 
    } else {
      qp <- NA
      CD.detected = NA
      Control.detected = NA
    }
    new.out.row <- data.frame(
      "cmpd" = cmpd,
      "feat" = RC2$featnames[use.index[j-1]],
      "feat.mz" = RC2$fmz[use.index[j-1]],
      "feat.rt" = RC2$frt[use.index[j-1]],
      "t.test.pval" = tp,
      "fold.change" = fc,
      "chi.sq.pval" = qp,
      "CD.detected.proportion" = CD.detected,
      "Control.detected.proportion" = Control.detected
    )
    out <- rbind(out, new.out.row)
  }
  
}

boxplot(out$t.test.pval~out$cmpd)
boxplot(out$chi.sq.pval ~out$cmpd)

row.names(out) <- NULL
writexl::write_xlsx(out, paste0(out.dir, "/feature.detection.frequency.stats.all.features.xlsx"))
keep <- which((out$t.test.pval < 0.05 | out$chi.sq.pval < 0.05) & out$fold.change > 4)
out[keep,]



### MS1 viewer
cmpd <- 1736
mzs <- RC$findmain[[cmpd]]$details[[1]]$mz
ints <- RC$findmain[[cmpd]]$details[[1]]$int
adducts <- RC$findmain[[cmpd]]$details[[1]]$adduct
hover <- sapply(1:length(adducts), FUN = function(x) {
  paste0(
    'featname: ', RC$featname[RC$findmain[[cmpd]]$details[[1]]$feature.index[x]],
    '\nisogroup: ', RC$findmain[[cmpd]]$details[[1]]$isogr[x], 
        '\nisotope: ', RC$findmain[[cmpd]]$details[[1]]$iso[x], 
        '\ncharge: ', RC$findmain[[cmpd]]$details[[1]]$charge[x],
        '\nmz: ', round(RC$findmain[[cmpd]]$details[[1]]$mz[x], 5),
        '\nintensity: ', round(RC$findmain[[cmpd]]$details[[1]]$int[x], 2)
  )
}
)

tmp <- plot_ly(x = mzs, y = ints, 
               type="scatter", 
               showlegend = FALSE,
               hoverinfo = "text", 
               size = 100, 
               mode="none", 
               text = hover) %>%
  layout(title = paste(RC$cmpd[cmpd], ":", round(RC$clrt[cmpd], 1), "sec") ,
    xaxis = list(zeroline = F,
                      showgrid = F,
                      showticklabels = F),
         yaxis = list(zeroline = F,
                      showgrid = F,
                      showticklabels = F)
  )
for(i in 1:length(mzs)) {
  tmp <- tmp %>%
    add_segments(x = mzs[i], xend = mzs[i],hoverinfo = "text", text = hover[i],
                 y = 0, yend = ints[i],
                 showlegend = FALSE, line = list(width = 1, color = "forestgreen"))
}

tmp <- tmp %>% add_trace(
  x = mzs, y = ints + (max(ints)*.02), 
  type="scatter",  showlegend = TRUE, 
  text = adducts, textfont = list(size = 18), 
  mode="text", fill = 'white', hoverinfo = "text",
  text = hover)
tmp
dda.feat.names <- xd.msms$feature_id
cmpd.feat.names <- RC$featname[RC$findmain[[cmpd]]$details[[1]]$feature.index]
dda.feat.names.sub <- dda.feat.names[dda.feat.names %in% cmpd.feat.names]
dda.feat.names.table <- table(dda.feat.names.sub)
length(RC$dda.spectra[[cmpd]])



