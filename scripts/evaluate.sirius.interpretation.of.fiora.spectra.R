library(Spectra)
library(enviPat)
data(isotopes)

## load positive ion mode fiora predicted spectra
## from pubchembio 20250828
## so we can add mass error noise
# load("C:/Temp/20250828/fiora/fiora.pos.Rdata")

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
# keep <- which(fiora.pos$inchikey %in% inchikeys)
# fiora.pos.sub <- fiora.pos[keep]  ##3627 spectra, which includes three collision energies for each compound
# save(fiora.pos.sub, file = "C:/Temp/fiora.sirius/output.summaries/fiora.pos.sub.Rdata")
load("C:/Temp/20250828/fiora/fiora.pos.sub.Rdata")


names <- fiora.pos.sub$title
formulas <- fiora.pos.sub$formula
inchikey.first.block <- fiora.pos.sub$inchikey.first.block
inchikey <- fiora.pos.sub$inchikey
precursorMz <- fiora.pos.sub$precursorMz
mimw <- precursorMz - isotopes[which(isotopes$isotope == "1H"), "mass"]

res.files <- list.files('C:/Temp/fiora.sirius/output.summaries/', pattern = "cmd2.sirius")
prec.rms <- as.numeric(gsub("pt", ".", sapply(1:length(res.files), FUN = function(x) unlist(strsplit(res.files[x], "_"))[2])))
prod.rms <- as.numeric(gsub("pt", ".", sapply(1:length(res.files), FUN = function(x) unlist(strsplit(res.files[x], "_"))[4])))

out <- data.frame(
  prec.rms = vector(mode = 'numeric'),
  prod.rms = vector(mode = 'numeric'), 
  name = vector(mode = 'character'), 
  inchikey = vector(mode = 'character'),
  inchikey.first.block = vector(mode = 'character'),
  precursorMz = vector(mode = 'numeric'), 
  mimw = vector(mode = 'numeric'),
  formula.rank = vector(mode = 'numeric'),
  formula.total = vector(mode = 'numeric'),
  formula.score = vector(mode = 'numeric'),
  formula.score.best = vector(mode = 'numeric'),
  formula.score.worst = vector(mode = 'numeric')
)

for(i in 1:length(res.files)) {
  dir <- paste0('C:/Temp/fiora.sirius/output.summaries/', gsub(".sirius", "", res.files[i]))
  
  ## read in formula and structure results for each configuration
  ## for each source spectrum, determine its ranked position in sirius results
  form <- read.delim(paste0(dir, "/formula_identifications_top-200.tsv"))
  form.names <- unique(sapply(1:nrow(form), FUN = function(x) unlist(strsplit(form$mappingFeatureId[x], "_HCD_"))[2]))
  #  struc <- read.delim(paste0(dir, "/structure_identifications_top-200.tsv"). )
  for(j in 1:length(names)) {
    form.sub <- form[grepl(names[j], form$mappingFeatureId, fixed = TRUE),]
    struc.sub <- struc[grepl(names[j], struc$mappingFeatureId, fixed = TRUE),]
  }
  
}

