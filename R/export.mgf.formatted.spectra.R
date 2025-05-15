#' export.mgf.formatted.spectra
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param output.dir character: path to write MGF formatted spectra to a new directory called 'spectra' is created (if it does not exist), and within it a directory called 'mgf'.  all spectra files are written to the 'spectra/ms' directory within the specified directory, or the working directory, if not specified. 
#' @param out.file logical: if TRUE, all spectra are exported in a single mgf file.  Else each spectrum is in its own file. default = TRUE
#' @param out.file.name character: file name of mgf file when 'one.file' = TRUE. default = 'spectra.mgf'
#' @details This function exports spectra in .mgf format, suitable for importing into GNPS, NIST MSsearch, etc.
#' @return  nothing - files written to disk
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept xcms
#' @author Corey Broeckling
#' @export 
#' 
#' 
export.mgf.formatted.spectra <- function(
    ramclustObj = NULL,
    output.dir = NULL, 
    one.file = TRUE,
    one.file.name = "spectra.mgf"
) {
  
  if(is.null(ramclustObj$mse.spectra) & is.null(ramclustObj$dda.spectra)){
    stop("neither mse nor dda spectra have been assigned. Please use one of 'assign.xcms.dda.peaks' or 'assign.xcms.mse.peaks' first.", '\n')
  }
  
  
  if(is.null(output.dir)) {
    output.dir <- getwd()
  }
  if(!dir.exists(output.dir)) {
    dir.create(output.dir)
  }
  if (!dir.exists(paste0(output.dir, "/spectra"))) {
    dir.create(paste0(output.dir, "/spectra"))
  }
  if (!dir.exists(paste0(output.dir, "/spectra/mgf"))) {
    dir.create(paste0(output.dir, "/spectra/mgf"))
  }
  
  if(!is.null(ramclustObj$dda.spectra)) {
    spec <- ramclustObj$dda.spectra
  }
  if(!is.null(ramclustObj$mse.spectra)) {
    if(any(ls() == 'spec')) {
      spec <- c(spec, ramclustObj$mse.spectra)
    } else {
      spec <- ramclustObj$mse.spectra
    }
    
  }
  
  if(one.file) {
    out.list <- as.list(rep(NA, 0))
  }
  
  for (i in 1:length(spec)) {
    if(!is.null(spec[[i]])) {
      if(length(spec[[i]]) == 0) next
      for(j in 1:length(spec[[i]])) {
        sp.name <- paste(
          spec[[i]][[j]]$ramclustr.cmpd, 
          spec[[i]][[j]]$xcms.name, 
          spec[[i]][[j]]$adduct.type,
          spec[[i]][[j]]$spectrum.type, sep = ".")
        out <- paste("BEGIN IONS", '\n',
                     paste0("TITLE=", sp.name), '\n',
                     paste0("PEPMASS=", round(spec[[i]][[j]]$precursor.mz, 5)), '\n',
                     paste0("CHARGE=", spec[[i]][[j]]$charge), '\n',
                     paste0("MSLEVEL=2"), '\n', sep = "")
        
        # out <- paste(out, '\n', 
        #              ">collision ", spec[[i]][[j]]$collision.energy, '\n',
        #              sep = "")
        for (k in 1:nrow(spec[[i]][[j]]$spectrum)) {
          out <- paste(out, round(spec[[i]][[j]]$spectrum$mz[k], 5), " ", round(spec[[i]][[j]]$spectrum$int[k], 2), "\n",
                       sep = "")
        }
        out <- paste0(out, "END IONS", '\n')
        if(one.file) {
          out.list[[length(out.list)+1]] <- out
        } else {
          write(out, file = paste0(output.dir, "/spectra/mgf/", sp.name,
                                 ".mgf"))
          }
        
        rm(out)
      }
    }
  }
  if(one.file) {
    out <- unlist(out.list)
    out <- paste(out, collapse = '\n')
    write(out, file = paste0(output.dir, "/spectra/mgf/spectra.mgf"))
  }
  
}
