library('Spectra')
library('MsBackendMsp')

msp.files <- paste0("R:/RSTOR-PMF/Software/current_NIST_software/NIST17/MSSEARCH/", c("nist_msms.MSP", "nist_msms2.MSP")) 
nist17 <- Spectra::backendInitialize(MsBackendMsp::MsBackendMsp(), msp.files)
save(nist17, file = "R:/RSTOR-PMF/Software/current_NIST_software/NIST17/MSSEARCH/nist17-msms.Rdata")

spectraVariables(nist17)
length(unique(nist17$Name))

load("R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.all.Rdata")
nist17.in.nist23 <- nist17$Name %in% nist.hr.all$name
nist.17.cas <- nist17$CASNO
nist.23.cas <- nist.hr.all$CAS.
nist.23.cas <- sapply(1:length(nist.23.cas), FUN = function(x) {unlist(strsplit(nist.23.cas[x], ";",  fixed = TRUE))[1]})
nist.23.cas <- gsub("-", "", nist.23.cas)

nist.17.no <- nist17$NISTNO
nist.23.no <- nist.hr.all$CAS.
nist.23.no <- sapply(1:length(nist.23.no), FUN = function(x) {unlist(strsplit(nist.23.no[x], "NIST#: ",  fixed = TRUE))[2]})

table(nist.17.cas %in% nist.23.cas)
table(nist.17.no %in% nist.23.no)
table(nist.23.cas %in% nist.17.cas)
table(nist.23.no %in% nist.17.no)

## use only spectra from nist23 which are not part of the nist17 cas set. 
## using NIST NO would exclude spectra, but potentially allow for using the same compound
## name is just messy
nist.23.cas <- nist.hr.all$CAS.
nist.23.cas <- sapply(1:length(nist.23.cas), FUN = function(x) {unlist(strsplit(nist.23.cas[x], ";",  fixed = TRUE))[1]})
nist.23.cas.dedashed <- gsub("-", "", nist.23.cas)
nist.17.cas <- nist17$CASNO
use.cas.23 <- unique(nist.23.cas[which(nist.23.cas.dedashed %in% nist.17.cas)])
sink("R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/cas.numbers.unique.to.23.vs.17.txt")
cat(paste0(use.cas.23, collapse='\n'))
sink()

nist.23.cas <- nist.hr.all$CAS.
nist.23.cas <- sapply(1:length(nist.23.cas), FUN = function(x) {unlist(strsplit(nist.23.cas[x], ";",  fixed = TRUE))[1]})


table(nist17.in.nist23)
nist.hr.all$CAS.
nist17.in.nist23.cas <- nist17$CASNO %in% nist.hr.all$CAS.
table(nist17.in.nist23.cas)
nist23.in.nist17 <- nist.hr.all$name %in% nist17$Name
table(nist23.in.nist17)

## exported a test msp for one spectrum from NIST23 MSSearch
test.msp <- paste0("R:/RSTOR-PMF/Software/current_NIST_software/NIST17/MSSEARCH/", c("test.MSPEC")) 
nistTest <- Spectra::backendInitialize(MsBackendMsp::MsBackendMsp(), test.msp)
spectraVariables(nistTest)

## still no smiles or inchikey,  have to use CAS and/or Name for mapping between nist17 and nist23.  kinda ridiculous, actually. no point in using the click bot export.