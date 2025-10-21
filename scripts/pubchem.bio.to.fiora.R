library(pubchem.bio)
library(csu.arcbio.tools)
packageVersion('pubchem.bio')
local.pubchem.directory <- "C:/Temp/20250828"  ## i suggest this naming scheme, but feel free to chose your own.
load(paste0(local.pubchem.directory, "/pc.bio.Rdata"))
pc.bio <- build.element.count(pubchem.bio.object = pc.bio)
pc.bio.fiora <- pubchem.bio.to.fiora(
  pubchem.bio.object = pc.bio,
  param.table = data.frame(
    Precursor_type = c(rep("[M+H]+", 3), rep("[M-H]-", 3)), 
    CE = c(10, 20, 40, 10, 20, 40), 
    Instrument_type = rep("HCD", 6)),
  out.dir = 'C:/Temp/for.fiora/'
)




## at this point you need to run Fiora on Linux
## once complete, you can import files
## into 'Spectra' object to enable
## spectral searching


## we then import the spectra into a 'Spectra' object
mgf.dir <- "C:/Temp/20250828/fiora/"
fiora.lib <- fiora.mgf.to.spectra(mgf.dir = mgf.dir)

## we also want to bring in metadata from the pubchem.bio object
## should probably put this into the fiora.mgf.to.spectra function
load(paste0(local.pubchem.directory, "/pc.bio.Rdata"))

## load pos, then neg fiora data
load(paste0(mgf.dir, "/fiora.pos.Rdata"))
f.sm <- fiora.pos$precursor.smiles
f.sm.mtch <- match(f.sm, pc.bio$smiles)
fiora.pos$inchikey <- pc.bio$inchikey[f.sm.mtch]
fiora.pos$inchikey.first.block <- pc.bio$inchikey.first.block[f.sm.mtch]
fiora.pos$XLogP <- pc.bio$XLogP[f.sm.mtch]
fiora.pos$nAcid <- pc.bio$nAcid[f.sm.mtch]
fiora.pos$nBase <- pc.bio$nBase[f.sm.mtch]
fiora.pos$TopoPSA <- pc.bio$TopoPSA[f.sm.mtch]
save(fiora.pos, file = paste0(mgf.dir, "/fiora.pos.Rdata"))

## load pos, then neg fiora data
load(paste0(mgf.dir, "/fiora.neg.Rdata"))
f.sm <- fiora.neg$precursor.smiles
f.sm.mtch <- match(f.sm, pc.bio$smiles)
fiora.neg$inchikey <- pc.bio$inchikey[f.sm.mtch]
fiora.neg$inchikey.first.block <- pc.bio$inchikey.first.block[f.sm.mtch]
fiora.neg$XLogP <- pc.bio$XLogP[f.sm.mtch]
fiora.neg$nAcid <- pc.bio$nAcid[f.sm.mtch]
fiora.neg$nBase <- pc.bio$nBase[f.sm.mtch]
fiora.neg$TopoPSA <- pc.bio$TopoPSA[f.sm.mtch]
save(fiora.neg, file = paste0(mgf.dir, "/fiora.neg.Rdata"))
