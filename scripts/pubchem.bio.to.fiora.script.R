library(csu.arcbio.tools)
library(pubchem.bio)
local.pubchem.directory <- 'C:/Temp/20251217/'
load(paste0(local.pubchem.directory, "/pc.bio.Rdata"))

pubchem.bio.to.fiora(pubchem.bio.object = pc.bio, max.cmpds.per.file = 1000, 
                                       param.table = data.frame(
                                         Precursor_type = c(rep("[M+H]+", 3), rep("[M-H]-", 3)), 
                                         CE = c(10, 20, 40, 10, 20, 40), 
                                         Instrument_type = "HCD"), min.elements = 8,
                                       min_prob = 0.05, annotation = TRUE, 
                                       rm.charged = TRUE, mw.range = c(50, 2000), out.dir = "C:/Temp/20251217_fiora/")

## still not removing 260_[M-H]-_40_HCD	Br

