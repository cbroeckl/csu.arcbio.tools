#' pubchem.bio.to.fiora
#'
#' takes a pubchem.bio object (data.table) and writes out input files for fiora. 
#' @details takes a pubchem.bio object (data.table) and writes out input files for fiora. optionally splits into multiple files for RAM management, writes out text file of commands suitable for running batch processing. 
#' @param pubchem.bio.object R data.table, generally produced by build.pubchem.bio.
#' @param max.cmpds.per.file integer. default = 1000.  use more if you have lots of RAM, less if you have little. 
#' @param param.table data frame. should contain headers 'Precursor_type', 'CE', and 'Instrument_type'.  Each row becomes a parameter set to predict. See Fiora documentation for valid options.
#' @param min_prob numeric. fiora option for minimum probability for a fragment ion to be reported. default = 0.001 in Fiora, 0.01 in this function.
#' @param annotation logical. should fragment ion structures be reported in output mgf format? default = TRUE.
#' @param rm.charged logical. should charged molecules be removed?  default = FALSE.  
#' @param mw.range numeric vector, length two.  c(50, 1200) by default. only compounds in this range are exported.  
#' @param out.dir path.  full path to directory files are to be written to. Directory must already exist.
#' @return data.frame of Fiora input parameters.  Files also sent to output directory 'out.dir'.    
#' @author Corey Broeckling
#' 
#' @export
#' 
#'
pubchem.bio.to.fiora <- function(
    pubchem.bio.object = NULL,
    max.cmpds.per.file = 1000,
    param.table = data.frame(Precursor_type = "[M+H]+", CE = 20, Instrument_type = "HCD"), 
    min_prob = 0.01,
    annotation = TRUE,
    rm.charged = FALSE,
    mw.range = c(50,1200),
    out.dir = NULL
    ) {
  
  if(!dir.exists(out.dir)) {
    stop('out.dir must point to a valid directory', '\n')
  }
  
  pubchem.bio.object <- pubchem.bio.object[
    which(pubchem.bio.object$monoisotopic.mass > mw.range[1] & pubchem.bio.object$monoisotopic.mass < mw.range[2]),
    ]
  
  if(rm.charged) {
    charged <- grepl("+", pubchem.bio.object$formula, fixed = TRUE) | grepl("-", pubchem.bio.object$formula, fixed = TRUE)
    pubchem.bio.object <- pubchem.bio.object[!charged,]
  }
  
  out <- data.frame(
    Name = rep("", 0),
    SMILES = rep("", 0),
    Precursor_type = rep("", 0),
    CE = rep(0, 0),
    Instrument_type = rep("", 0)
  )
  
  for(i in 1:nrow(param.table)) {
    tmp <- data.frame(
      Name = paste(pubchem.bio.object$cid, param.table$Precursor_type[i], param.table$CE[i], param.table$Instrument_type[i], sep = "_"),
      SMILES = pubchem.bio.object$smiles,
      Precursor_type = rep(param.table$Precursor_type[i], nrow(pubchem.bio.object)),
      CE = rep(param.table$CE[i], nrow(pubchem.bio.object)),
      Instrument_type = rep(param.table$Instrument_type[i], nrow(pubchem.bio.object))
    )
    
    out <- rbind(tmp, out)
  }
  
  suppressWarnings(chnks <- split(1:nrow(out), ceiling(seq_along(1:nrow(out))/max.cmpds.per.file)))
  bat.out <- rep("", 0)
  
  in.names <- paste0(
    "fiora/data/fiora.input_", 
    formatC(1:length(chnks), digits = (nchar(length(chnks))-1), flag = 0),
    ".csv"
    )
  out.names <-  paste0(
    "fiora/data/fiora.output_", 
    formatC(1:length(chnks), digits = (nchar(length(chnks))-1), flag = 0),
    ".mgf"
  )
  
  
  for(i in 1:length(chnks)) {
  # for(i in 1:3) {
    do <- chnks[[i]]
    out.sub <- out[do,]
    write.csv(out.sub, file = paste0(out.dir, "/", basename(in.names[i])), row.names = FALSE) 
    bat.out <- c(
      bat.out,
      paste("fiora-predict", "-i", in.names[i], "-o", out.names[i], 
            if(annotation) {"--annotation"}, 
            paste0("--min_prob=", min_prob))
    )
  }
  
  sink(paste0(out.dir, "/", "fiora.script.sh"))
  cat(paste(bat.out, collapse = '\n'))
  sink()
  
  message("After copying file to linux server, ensure  the '.sh' file is saved in 'unix' format (MobaXTerm converts to DOS upon transfer). ")
  message("You will then need to make the file executable using 'chmod +x data/fiora.script.sh', for example.")
  return(out)
}
