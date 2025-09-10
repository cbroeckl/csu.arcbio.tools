#' fiora.mgf.to.spectra
#'
#' imports fiora .mgf files and formats, with naming scheme from 'pubchem.bio.to.fiora', into R as a 'Spectra' object, saves as new .Rdata object on disk
#' @details takes a pubchem.bio object (data.table) and writes out input files for fiora. optionally splits into multiple files for RAM management, writes out text file of commands suitable for running batch processing. 
#' @param mgf.dir character file path.  directory containing fiora output .mgf files.  
#' @param split.by.polarity logical.  if true, spectra are split into positive ('1') and negative ('0') polarities before saving. default = TRUE.
#' @param out.dir character. full file path to save 'Spectra' object(s) to.  
#' @return spectra object.  File(s) also sent to output directory 'out.dir'.    
#' @author Corey Broeckling
#' 
#' @export
#' 
#'
pubchem.bio.to.fiora <- function(
    mgf.dir = NULL,
    split.by.polarity = TRUE,
    out.dir = NULL
    ) {

  requireNamespace('Spectra')
  requireNamespace('MsBackendMgf')
  
  if(!dir.exists(mgf.dir)) stop("mgf.dir does not exist: ", mgf.dir, '\n')
  mgf.files <- list.files(mgf.dir, pattern = "mgf", full.names = TRUE)
  
  if(is.null(out.dir)) {
    out.dir <- mgf.dir
  } else {
    if(!dir.exists(out.dir)) {
      dir.create(out.dir)
    }
  }
  
  normalize.intensity <- function(x, ...) {
    maxint <- max(x[, 2], na.rm = TRUE)
    x[, 2] <- as.integer(100 * x[, 2] / maxint)
    x
  }
  
  mapping = c(
    title = "TITLE",
    precursorMz = "PRECURSOR_MZ",
    precursor.smiles = "SMILES",
    adduct = "PRECURSORTYPE",
    collisionEnergy = "COLLISIONENERGY", 
    activation.type = "INSTRUMENTTYPE",
    formula = "FORMULA",
    comment = "COMMENT"
  )

  for(i in 1:length(mgf.files)) {
    cat(i, " ")
    ## read in mgf file
    
    s <- backendInitialize(MsBackendAnnotatedMgf(), mgf.files[i], mapping = mapping)
    add.type <- s$adduct
    polarity <- rep(1, length(add.type))
    is.neg <- grep("[M-H]-", add.type, fixed = TRUE)
    if(any(is.neg)){
      polarity[is.neg] <- 0
    }
    polarity(s) <- polarity
    plotSpectra(s[1])
    nm <- s$title
    cid <- as.numeric(sapply(1:length(nm), FUN = function(x) unlist(strsplit(nm[x], "_"))[1]))
    s$cid <- as.integer(cid)
    s$precursorCharge <- 1L
    s <- Spectra::Spectra(s)
    
    ## clean spectra
    s <- filterMzRange(s, c(50,2000))
    s <- addProcessing(s, normalize.intensity)
    
    ## apply processing
    s <- applyProcessing(s)
    
    ## bind spectra to prior spectra
    if(i == 1) {
      fiora.library <- s
    } else {
      fiora.library <- concatenateSpectra(fiora.library, s)
    }
    rm(s)
    gc()
    
  }
  
  fiora.library <- Spectra::setBackend(fiora.library, backend = MsBackendMemory())
  save(fiora.library, file = paste0(mgf.dir, "fiora.pos.Rdata"))
  
  if(split.by.polarity) {
    fiora.pos <- Spectra::filterPolarity(fiora.library, 1)
    fiora.neg <- Spectra::filterPolarity(fiora.library, 0)
    if(length(fiora.pos) > 0) save(fiora.pos, file = paste0(mgf.dir, "fiora.pos.Rdata"))
    if(length(fiora.neg) > 0) save(fiora.neg, file = paste0(mgf.dir, "fiora.pos.Rdata"))
    rm(fiora.library)
    fiora.library <- as.list(rep(NA, 2))
    names(fiora.library) <- c("pos", "neg")
    fiora.library[[1]] <- fiora.pos
    fiora.library[[2]] <- fiora.neg
  }
  
  return(fiora.library)
  
}
