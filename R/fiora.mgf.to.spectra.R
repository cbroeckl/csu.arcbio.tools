#' fiora.mgf.to.spectra
#'
#' imports fiora .mgf files and formats, with naming scheme from 'pubchem.bio.to.fiora', into R as a 'Spectra' object, saves as new .Rdata object on disk
#' @details takes a pubchem.bio object (data.table) and writes out input files for fiora. optionally splits into multiple files for RAM management, writes out text file of commands suitable for running batch processing. 
#' @param mgf.dir character file path.  directory containing fiora output .mgf files.  
#' @param split.by.polarity logical.  if true, spectra are split into positive ('1') and negative ('0') polarities before saving. default = TRUE.
#' @param out.dir character. full file path to save 'Spectra' object(s) to.  
#' @param return.library logical.  If true, return library to R console.  
#' @return spectra object.  File(s) also sent to output directory 'out.dir'.    
#' @author Corey Broeckling
#' 
#' @export
#' 
#'
fiora.mgf.to.spectra <- function(
    mgf.dir = NULL,
    split.by.polarity = TRUE,
    out.dir = NULL, 
    return.library = TRUE
) {
  
  requireNamespace('Spectra')
  requireNamespace('MsBackendMgf')
  
  
  if(!dir.exists(mgf.dir)) stop("mgf.dir does not exist: ", mgf.dir, '\n')
  mgf.files <- list.files(mgf.dir, pattern = "mgf", full.names = TRUE)
  # mgf.files <- mgf.files[1:1000]
  
  
  if(is.null(out.dir)) {
    out.dir <- mgf.dir
  } else {
    if(!dir.exists(out.dir)) {
      dir.create(out.dir)
    }
  }
  
  normalize.intensity <- function(x, ...) {
    maxint <- max(x[, 2], na.rm = TRUE)
    x[, 2] <- as.integer(1000 * x[, 2] / maxint)
    x
  }
  
  # ## read in all individual .mgf files and export to single file called 'all.mgf'
  # 
  # if(aggregate.files) {
  #   tmp <- as.list(rep(NA, length(mgf.files)))
  #   for(i in 1:length(mgf.files)) {
  #     tmp[[i]] <- paste(readLines(mgf.files[i]), collapse = '\n')
  #     # tmp <- c(tmp, tm)
  #     # cat(i, " ")
  #   }
  #   
  #   out <- sapply(tmp, paste, collapse= '\n')
  #   sink(paste0(mgf.dir, "/all.MGF"))
  #   cat(out)
  #   sink()
  #   rm(out)
  #   rm(tmp)
  #   gc()
  # }
  # 
  
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
  
  # if(aggregate.files) {
  #   fiora.library <- Spectra::backendInitialize(MsBackendMgf::MsBackendAnnotatedMgf(), paste0(mgf.dir, "/all.MGF"), mapping = mapping)
  # } else {
    fiora.library <- Spectra::backendInitialize(MsBackendMgf::MsBackendAnnotatedMgf(), mgf.files, mapping = mapping)
  # }
  # print(mgf.files)
  add.type <- fiora.library$adduct
  polarity <- rep(1L, length(add.type))
  is.neg <- grep("[M-H]-", add.type, fixed = TRUE)
  if(any(is.neg)){
    polarity[is.neg] <- 0L
  }
  fiora.library$polarity <- polarity
  nm <- fiora.library$title
  cid <- as.numeric(sapply(1:length(nm), FUN = function(x) unlist(strsplit(nm[x], "_"))[1]))
  fiora.library$cid <- as.integer(cid)
  fiora.library$precursorCharge <- 1L
  fiora.library <- Spectra::Spectra(fiora.library)
  
  ## clean spectra
  fiora.library <- Spectra::filterMzRange(fiora.library, c(50,2000))
  fiora.library <- Spectra::addProcessing(fiora.library, normalize.intensity)
  
  ## apply processing
  fiora.library <- Spectra::applyProcessing(fiora.library)
  
  gc()
  
  fiora.library <- Spectra::setBackend(fiora.library, backend = Spectra::MsBackendMemory())
  save(fiora.library, file = paste0(mgf.dir, "fiora.pos.Rdata"))
  
  if(split.by.polarity) {
    fiora.pos <- Spectra::filterPolarity(fiora.library, 1L)
    fiora.neg <- Spectra::filterPolarity(fiora.library, 0L)
    if(length(fiora.pos) > 0) save(fiora.pos, file = paste0(mgf.dir, "fiora.pos.Rdata"))
    if(length(fiora.neg) > 0) save(fiora.neg, file = paste0(mgf.dir, "fiora.neg.Rdata"))
    rm(fiora.library)
    fiora.library <- as.list(rep(NA, 2))
    names(fiora.library) <- c("pos", "neg")
    fiora.library[[1]] <- fiora.pos
    fiora.library[[2]] <- fiora.neg
  }
  
  if(return.library) {
    return(fiora.library)
  }
  
}
