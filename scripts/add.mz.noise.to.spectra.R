library(Spectra)

## load positive ion mode fiora predicted spectra
## from pubchembio 20250828
## so we can add mass error noise
load("C:/Temp/20250828/fiora/fiora.pos.Rdata")

## get all structures from library, sample 0.1% of all randomly.  we will 
## use these structures for both positive and negative ionization
## mode.
n.inchikeys <- length(unique(fiora.pos$inchikey))
set.seed(123456)
inchikeys <- sample(fiora.pos$inchikey, round(n.inchikeys/1000))

## then also add a selection of structural isomers - same first block of inchikey
isomers <- sort(table(fiora.pos$inchikey.first.block), decreasing = TRUE)
sum(isomers[1:30])
isomer.classes <- c(
  trisaccharides = "FBJQEBRMDXPWNX",
  neurofurans = "UCBGJOSAVMKADE",
  sterols = "HCXVJBMSMIARIN",
  caffeolquinates = "CWVRJTMFETXNAD"
)
sum(isomers[isomer.classes])
## adds 381 compounds in four isomer groups
isomer.inchikeys <- fiora.pos$inchikey[fiora.pos$inchikey.first.block %in% isomer.classes]
## add to target inchikey list
inchikeys <- unique(c(inchikeys, isomer.inchikeys))

## create a subset library containing only those inchikeys
keep <- which(fiora.pos$inchikey %in% inchikeys)
fiora.pos.sub <- fiora.pos[keep]  ##3627 spectra, which includes three collision energies for each compound

## we will be adding mass error based on the RMS mass error concept
## for a given mass, rms mass error can generated as below
mz <- 1000
n <- 1000
ppm <- 0.5
actual <- rep(mz, n)
measured <- rnorm(n, mean = mz, sd = ppm*mz/1000000)
hist(actual - measured)
rmse <- sqrt(mean((actual  - measured)^2))
rmse.ppm <- 1000000*(rmse)/mz
rmse.ppm

## now to do this for all precursor and fragment ion masses
## this is an example.  we will loop through various ppm levels

## for precursors
actual <- fiora.pos.sub$precursorMz
measured <- rnorm(mean = actual, sd = ppm*actual/1000000, n = length(actual))
## double check that PPM centers on zero with rmse roughly the same as the target ppm
ppm.dev <- 1000000*(actual  - measured)/actual
rmse.ppm <- sqrt(mean((ppm.dev)^2))
hist(ppm.dev, main = paste0("measured.ppm = ", round(rmse.ppm, 2)))

## for product ions, one spectrum at a time
for(i in 1:1) {
  actual <- fiora.pos.sub@backend@peaksData[[i]][,"mz"]
  measured <- rnorm(mean = actual, sd = ppm*actual/1000000, n = length(actual))
  fiora.pos.sub@backend@peaksData[[i]][,"mz"] <- measured
}

## also create sirius command line for each condition.  base command is
## this line has worked in the past.  need to sub in directory/name and ppm error variables
sirius.cmd <- 'sirius -i "C:\Temp\fiora.sirius\prec_0.5_prod_0.5" -o "C:\Temp\fiora.sirius\output.summaries\prec_0pt5_prod_0pt5.sirius" --cores=8 config --IsotopeSettings.filter=true --CandidateFormulas=, --FormulaSettings.enforced=H,C,N,O,F,P,I --Timeout.secondsPerInstance=600 --Timeout.secondsPerTree=600 --AlgorithmProfile=qtof --AdductSettings.ignoreDetectedAdducts=false --AdductSettings.prioritizeInputFileAdducts=true --UseHeuristic.useHeuristicAboveMz=300 --IsotopeMs2Settings=IGNORE --MS1MassDeviation.allowedMassDeviation=0.5ppm --MS2MassDeviation.allowedMassDeviation=0.5ppm --FormulaSearchSettings.performDeNovoBelowMz=0 --FormulaSearchSettings.applyFormulaConstraintsToDatabaseCandidates=false --EnforceElGordoFormula=true --NumberOfCandidatesPerIonization=20 --AdductSettings.fallback=[[M+H]+,[M+Na]+,[M+K]+,[M+H3N+H]+] --FormulaSearchSettings.performBottomUpAboveMz=Infinity --FormulaSearchSettings.applyFormulaConstraintsToBottomUp=false --UseHeuristic.useOnlyHeuristicAboveMz=650 --FormulaSearchDB=PUBCHEM,METACYC,BloodExposome,CHEBI,COCONUT,DSSTox,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,PUBMED,SUPERNATURAL,TeroMol,YMDB --AdductSettings.enforced=, --FormulaSettings.detectable=B,S,Cl,Se,Br --NumberOfCandidates=20 --FormulaResultThreshold=true --ExpansiveSearchConfidenceMode.confidenceScoreSimilarityMode=OFF --StructureSearchDB=PUBCHEM,METACYC,BloodExposome,CHEBI,COCONUT,DSSTox,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,PUBMED,SUPERNATURAL,TeroMol,YMDB --RecomputeResults=true formulas fingerprints classes structures summaries --top-k-summary=200'

## set precursor and product errors
## prec = product for most, 
## but also explore precursor and product independently
## to explore value of product ion accuracy, specifically
ppm.conditions <- data.frame(
  "precursor" = c(5, 2, 1, 0.5, 0.1, 5,    2,   1,   1),
  "product" =   c(5, 2, 1, 0.5, 0.1, 0.5,  0.5, 0.5, 0.1)
)

out.bat <- vector(mode = "character", length = 0)

for(x in 1:nrow(ppm.conditions)) {
  ## get precursor and product masses
  prec <- ppm.conditions$precursor[x]
  prod <- ppm.conditions$product[x]
  
  ## configure sirius call for bat output
  base.dir <- "C://Temp//fiora.sirius//"
  if(!dir.exists(base.dir)) dir.create(base.dir)
  on.name <- paste0("prec_", prec, "_prod_", prod)
  in.dir <-  paste0(base.dir, on.name)
  if(!dir.exists(in.dir)) dir.create(in.dir)
  out.dir <- paste0(base.dir, "output.summaries//")
  if(!dir.exists(out.dir)) dir.create(out.dir)
  out.name <- paste0(out.dir, gsub(".", "pt", on.name, fixed = TRUE), ".sirius")
  sirius.template <- paste0(
    'sirius -i ', in.dir, ' -o ', out.name, 
    ' config --AlgorithmProfile=qtof ', 
    ' --MS1MassDeviation.allowedMassDeviation=', 3*prec, 'ppm --MS2MassDeviation.allowedMassDeviation=', 3*prod, 
    'ppm --SpectralSearchDB=PUBCHEM,METACYC,BloodExposome,CHEBI,COCONUT,DSSTox,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,PUBMED,SUPERNATURAL,TeroMol,YMDB --AdductSettings.fallback=[[M+H]+] --FormulaSettings.enforced=H,C,N,O,P --IdentitySearchSettings.precursorDeviation=20.0ppm --FormulaSearchSettings.performBottomUpAboveMz=0 --FormulaSearchDB=, --StructureSearchDB=METACYC,BloodExposome,CHEBI,COCONUT,DSSTox,FooDB,GNPS,HMDB,HSDB,KEGG,KNAPSACK,LOTUS,LIPIDMAPS,MACONDA,MESH,MiMeDB,NORMAN,PLANTCYC,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONFOOD,PUBCHEMANNOTATIONSAFETYANDTOXIC,PUBMED,SUPERNATURAL,TeroMol,YMDB formulas fingerprints classes structures summaries --top-k-summary=200'
  )
  out.bat <- c(out.bat, sirius.template, '\n')
  
  ## create copy of sub library
  tmp.lib <- fiora.pos.sub
  
  ## add noise to precursor ions
  actual <- tmp.lib$precursorMz
  measured <- rnorm(mean = actual, sd = ppm*actual/1000000, n = length(actual))
  tmp.lib$precursorMz <- measured
  
  ## add noise to product ions
  for(i in 1:length(tmp.lib)) {
    actual <- tmp.lib@backend@peaksData[[i]][,"mz"]
    measured <- rnorm(mean = actual, sd = ppm*actual/1000000, n = length(actual))
    tmp.lib@backend@peaksData[[i]][,"mz"] <- measured
  }
  
  ## create output spectrum
  for(i in 1:length(tmp.lib)) {
    out.ms <- paste0(
      ">compound ", tmp.lib$title[i], '\n',
      ">parentmass ", tmp.lib$precursorMz[i], '\n',
      ">ionization ", tmp.lib$adduct[i], '\n',
      ">charge ", tmp.lib$precursorCharge[i], '\n', '\n',
      # ">ms1 ", '\n',
      # ">parentmass ", tmp.lib$precursorMz[i], " ", 1000, '\n', '\n',
      ">collision ", tmp.lib$collisionEnergy[i], '\n'
    )
    sp <- cbind(
      mz = unlist(mz(tmp.lib[i])),
      intensity = unlist(intensity(tmp.lib[i]))
    )
    for(j in 1:nrow(sp)){
      out.ms <- paste0(out.ms, round(sp[j,1], 5), " ", sp[j,2], '\n')
    }
    
    sink(paste0(in.dir, "//", tmp.lib$title[i], ".ms"))
    cat(out.ms)
    sink()
    
  }
}

sink("C://Temp//fiora.sirius//sirius.bat")
cat(out.bat)
sink()


