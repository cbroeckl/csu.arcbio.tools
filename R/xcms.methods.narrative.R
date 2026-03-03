
xcms.processing.narrative <- function(
    xcmsObj = xd,
    out.file = "xcms-ms.methods.txt") {
  
  out.tmp <- paste(
    "Data processing was performed in R", paste0(R.Version()$major, ".", R.Version()$minor), "running XCMS version", paste0(package.version('xcms'), "."),
    xcmsObj@processHistory[[1]]@type,
    "was performed using the following", as.character(class(xcmsObj@processHistory[[1]]@param)),  "parameters: ")
  param.names <- slotNames(xcmsObj@processHistory[[1]]@param)
  for(j in 1:length(param.names)) {
    out.tmp <- paste0(" ", out.tmp, param.names[j], "=", paste0(slot(xcmsObj@processHistory[[1]]@param, param.names[j]), collapse = ","))
    if(j < length(param.names)) {
      out.tmp <- paste0(out.tmp, "; ")
    } else {
      out.tmp <- paste0(out.tmp, ".")
    }
  }
  out.tmp <- trimws(out.tmp)
  
  for(i in 2:length(xcmsObj@processHistory)) {
    out.tmp <- paste(out.tmp, xcmsObj@processHistory[[i]]@type,
    "was performed using the following", as.character(class(xcmsObj@processHistory[[i]]@param)),  "parameters: ")
    
    param.names <- slotNames(xcmsObj@processHistory[[i]]@param)
    for(j in 1:length(param.names)) {
      tmp.val <- slot(xcmsObj@processHistory[[i]]@param, param.names[j])
      if(is.function(tmp.val)) {
        tmp.val <- capture.output(print(tmp.val))[2]
      } else {
        if(length(tmp.val) > 2) {next}
      }
      out.tmp <- paste0(" ", out.tmp, param.names[j], "=", paste0(tmp.val, collapse = ","))
      if(j < length(param.names)) {
        out.tmp <- paste0(out.tmp, "; ")
      } else {
        out.tmp <- paste0(out.tmp, ".")
      }
    }
  }

  
  out.tmp <- trimws(out.tmp)
  
  sink(out.file)
  cat(out.tmp)
  sink()
}