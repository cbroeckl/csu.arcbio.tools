rc.annotate.pubchem.bio <- function(
    ramclustObj = NULL,
    pubchem.bio.object = NULL,
    ppm.tolerance = 5,
    library.ppm.tolerance = 10,
    resolution = 20000,
    rt.tol = 15, 
    ri.tol = 20,
    library.level1 = NULL,
    library.level2a = NULL,
    library.level2b = NULL,
    n.threads = 3,
    processingChunkSize = 1000
) {
  
  if(!requireNamespace("InterpretMSSpectrum", quietly = TRUE)) {
    stop("The use of this function requires package 'InterpretMSSpectrum'.")
  }
  if(!requireNamespace("enviPat", quietly = TRUE)) {
    stop("The use of this function requires package 'enviPat'.")
  }
  
  if(is.null(ramclustObj$mse.spectra) & is.null(ramclustObj$dda.spectra)) {
    warning("No MSMS spectra are contained in this ramclustObj. Annotation will be based only on accurate mass in the MS1 data.", '\n')
  } 
  
  # # load object and rename 'pubchem.bio.object' for simplicity
  # objects.1 <- ls()
  # load(pubchem.bio.file)
  # objects.2 <- ls()
  # new.object <- objects.2[!(objects.2 %in% objects.1)]
  # new.object <- new.object[-which(new.object == "objects.1")]
  # # envir <- environment()
  # # cat("environment =", environmentName(envir), '\n')
  # value <- get(new.object, envir=environment())
  # assign(x='pubchem.bio.object', value=value, envir=environment())
  # rm(new.object)
  # # cat(paste(ls()))
  
  ## make sure we have inchikey and first block in libraries 2a and 2b. use CID as reference
  lib.cid <- library.level2a$cid.parent
  mtch <- match(lib.cid, pubchem.bio.object$cid)
  lib.inchikey <- pubchem.bio.object$inchikey[mtch]
  lib.first.block <- pubchem.bio.object$inchikey.first.block[mtch]
  library.level2a$inchikey <- lib.inchikey
  library.level2a$inchikey.first.block <- lib.first.block
  
  lib.cid <- library.level2b$cid
  mtch <- match(lib.cid, pubchem.bio.object$cid)
  lib.inchikey <- pubchem.bio.object$inchikey[mtch]
  lib.first.block <- pubchem.bio.object$inchikey.first.block[mtch]
  library.level2b$inchikey <- lib.inchikey
  library.level2b$inchikey.first.block <- lib.first.block
  
  
  ## for each compound, assign putative formula/structure
  pubchem.bio.annotations <- as.list(rep(NA, length(ramclustObj$clrt)))
  for(i in 1:length(ramclustObj$clrt)) {
    
    pubchem.bio.annotations[[i]] <- list("annotation.table" = NA, "isotope.spectra" = NA)
    
    ## find neutral masses in pubchem.bio object consistent with assigned neutral masses from find.main
    fm <- ramclustObj$findmain[[i]]$summary
    fmm <- fm$neutral_mass
    ppm.err <- outer(fmm, pubchem.bio.object$monoisotopic.mass, FUN = "-")
    ppm.err <- 1000000*ppm.err/pubchem.bio.object$monoisotopic.mass
    keep <- which(abs(ppm.err) <= ppm.tolerance, arr.ind = TRUE)
    if(nrow(keep))
      annotations <- data.frame(
        fm[keep[,1],], 
        fm.hypothesis = keep[,1],
        pubchem.bio.object[keep[,2],], 
        ppm.error = round(ppm.err[keep],3)
      )
    annotations$ppm.score = round(exp(-(annotations$ppm.error^2)/(2*(ppm.tolerance ^ 2))), digits = 3)
    # annotations$ion.formula <- adductFormula()
    
    
    # use MS/MS fragment ions to help prioritize correct formula assignments
    # msms <- ramclustObj$ms
    
    pubchem.bio.annotations[[i]]$annotation.table <- annotations
  }
  
  ## from interpretMSSpectrum
  ## modified to return data for MetaboCoreUtils::adducts compatibility
  ## mass_multi = 1/nmol
  ## mass_add = massdiff
  ## formula_add new
  ## formula_sub new
  ## positive = logical charge > 0
  getRuleFromIonSymbol <- function(ions="[M+H]+") {
    checkSymbol <- function(ion) {
      regexpr("\\[[0-9]{0,2}M.*\\][0-9]{0,2}[\\+\\-]{1,2}", ion) != -1
    }
    shortCuts <- cbind(
      c("M+H", "M+Na", "M+K", "M+NH4", "M+", "M", "M-H","M+Cl-", "M-"),
      c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M]+", "[M-H]-","[M+Cl]-", "[M]-")
    )
    em <- 0.0005485799
    chemical_elements <- NULL
    utils::data(chemical_elements, envir=environment(), package="InterpretMSSpectrum")
    on.exit(rm(chemical_elements))
    out <- lapply(ions, function(ion) {
      if(ion %in% shortCuts[,1]) ion <- shortCuts[,2][ which(shortCuts[,1] == ion) ]
      if(!checkSymbol(ion)) stop("invalid ion")
      nmol <- sub(".*[^0-9M]([0-9]?M).*", "\\1", ion)
      nmol <- sub("M", "", nmol)
      nmol <- as.numeric(ifelse(nmol=="", 1, nmol))
      ch <- sub(".*[^0-9]([0-9]{0,2}[\\+\\-])$", "\\1", ion)
      sgn <- sub("[^\\+\\-]", "", ch)
      sgn <- ifelse(sgn=="+", 1, -1)
      ch <- sub("[\\+\\-]", "", ch)
      ch[ch==""] <- "1"
      ch <- as.numeric(ch)
      ch <- ch * sgn
      x <- ion
      x <- sub("^.*\\[", "", x)
      x <- sub("\\].*", "", x)
      x <- sub("[0-9]?M", "", x)
      starts <- gregexpr("[\\+\\-]", x)[[1]]
      ends <- c(starts[-1]-1, nchar(x))
      n <- length(starts)
      spl <- lapply(1:n, function(i) substr(x, starts[i], ends[i]))
      massdiff <- lapply(spl, function(y) {
        sgn <- sub("^([\\+\\-]).*", "\\1", y)
        sgn <- ifelse(sgn=="+", 1, -1)
        el <- sub("^[\\+\\-]", "", y)
        if (regexpr("^[0-9]+[A-Za-z]+", el) != -1) el <- gsub("([0-9]+)([A-Za-z]+)", "\\2\\1", el)
        #browser()
        #el <- tabulateElements(el)
        el <- InterpretMSSpectrum::CountChemicalElements(x=el)
        masses <- sapply(names(el), function(a) {
          chemical_elements[,2][ which(chemical_elements[,1] == a)[1] ]
        })
        return(
          (sum(masses * el) * sgn)
        )
      })
      addsub <- sapply(spl, function(y) {
        sgn <- sub("^([\\+\\-]).*", "\\1", y)
        sgn <- ifelse(sgn=="+", 1, -1)
        el <- sub("^[\\+\\-]", "", y)
        if (regexpr("^[0-9]+[A-Za-z]+", el) != -1) el <- gsub("([0-9]+)([A-Za-z]+)", "\\2\\1", el)
        #browser()
        #el <- tabulateElements(el)
        el <- InterpretMSSpectrum::CountChemicalElements(x=el)
        el[which(el == 1)] <- ""
        f <- paste0(paste0(names(el), el, collapase = ""), collapse = "")
        return(
          as.vector(f)
        )
      })
      
      
      spl <- unlist(spl)
      splp <- addsub[grepl("+", spl, fixed = TRUE)]
      if(length(splp) > 0) {
        formula_add <- splp[1]
        if(length(splp) > 1) {
          for(j in 2:length(splp)){
            formula_add <- MetaboCoreUtils::addElements(formula_add, splp[j])
          }
        }
      } else {
        formula_add <- "C0"
      }
      
      spln <- addsub[grepl("-", spl, fixed = TRUE)]
      if(length(spln) > 0) {
        formula_sub <- spln[1]
        if(length(spln) > 1) {
          for(j in 2:length(spln)){
            formula_sub <- MetaboCoreUtils::addElements(formula_sub, spln[j])
          }
        }
      } else {
        formula_sub <- "C0"
      }
      formula_sub <- gsub("-", "", formula_sub, fixed = TRUE)
      
      massdiff <- sum(unlist(massdiff), na.rm=TRUE) + ch * -em
      add.out <- data.frame(
        name=ion, mass_multi=1/nmol, mass_add=massdiff, 
        formula_add=formula_add, formula_sub=formula_sub, 
        charge=ch, positive = ch>0, stringsAsFactors = FALSE)
      return(add.out)
    })
    return(do.call("rbind", out))
  }
  
  all.adducts <- unique(unlist(sapply(1:length(pubchem.bio.annotations), FUN = function(x){unique(pubchem.bio.annotations[[x]]$annotation.table$adducthyp)})))
  ims.table <- getRuleFromIonSymbol(ions = all.adducts)
  row.names(ims.table) <- ims.table$name
  
  
  ## isotope pattern for precursor ion(s) score
  utils::data(isotopes, envir=environment(), package="enviPat")
  for(i in 1:length(pubchem.bio.annotations)) {
    ann <- pubchem.bio.annotations[[i]]$annotation.table
    rm <- c(grep("+", ann$formula, fixed = TRUE), grep("-", ann$formula, fixed = TRUE))
    if(length(rm) > 0) {
      ann <- ann[-rm,]
    }
    # cat(i, nrow(ann), '\n')
    if(nrow(ann) == 0) {
      pubchem.bio.annotations[[i]]$annotation.table$isotope.score <- NA
      next
    }
    add.forms.all <- paste(ann$formula, ann$adducthyp)
    add.forms <- unique(add.forms.all)
    forms <- sapply(1:length(add.forms), FUN = function(x) {unlist(strsplit(add.forms[x], " [", fixed = TRUE))[1]})
    adds <- sapply(1:length(add.forms), FUN = function(x) {paste0("[", unlist(strsplit(add.forms[x], " [", fixed = TRUE))[2])})
    ion.forms <- sapply(1:length(adds), FUN = function(x) {
      ## make tryCatch, NA if fail
      ## cat(x, " "); 
      
      tryCatch(
        expr = {
          MetaboCoreUtils::adductFormula(formulas = forms[x], adduct = ims.table[adds[x],])
        },
        warning = function(w) {
        },
        error = function(e) {
          NA
        },
        finally = {
          
        }
      )
    })
    
    ch <- sub(".*[^0-9]([0-9]{0,2}[\\+\\-])$", "\\1", ion.forms)
    ion.forms <- sub("^.*\\[", "", ion.forms)
    ion.forms <- sub("\\].*", "", ion.forms)
    if(any(is.na(ion.forms))) {ion.forms <- ion.forms[-which(is.na(ion.forms))]}
    if(any(ion.forms == "NULL")) {ion.forms <- ion.forms[-(which(ion.forms == "NULL")) ]}
    if(any(ion.forms == "NA")) {ion.forms <- ion.forms[-(which(ion.forms == "NA")) ]}
    if(length(ion.forms) == 0) {
      pubchem.bio.annotations[[i]]$annotation.table$isotope.score <- NA
      next()
    }
    suppressMessages({
      envs <- enviPat::envelope(enviPat::isopattern(emass = 0.00054857990924, isotopes = isotopes, chemforms = ion.forms, plotit = FALSE, verbose = FALSE), resolution = resolution, plotit = FALSE, verbose = FALSE)
    })
    cents <- enviPat::vdetect(envs, plotit = FALSE, verbose = FALSE)
    exp.spec <- ramclustObj$ms1.spectrum[[i]]
    add.form.sim <- rep(NA, length(ion.forms))
    spectra.pairs <- as.list(rep(NA, length(cents)))
    for(j in 1:length(cents)) {
      cents[[j]] <- cents[[j]][which(cents[[j]][,"abundance"] > 1),, drop = FALSE]
      mz.range <- range(cents[[j]][,"m/z"]) + c(-0.5, 0.5)
      sub.exp.spec <- exp.spec[((exp.spec$mz > mz.range[1]) & (exp.spec$mz < mz.range[2])),]
      sub.exp.spec$int <- 100*sub.exp.spec$int/max(sub.exp.spec$int)
      spd <- DataFrame(id = c("experimental", "theoretical"))
      spd$mz <- list(
        sub.exp.spec$mz, 
        cents[[j]][,1]
      )
      spd$intensity <- list(
        sub.exp.spec$int, 
        cents[[j]][,2]
      )
      sps <- Spectra::Spectra(spd)
      spectra.pairs[[j]] <- sps
      names(spectra.pairs)[j] <- ion.forms[j]
      # plotSpectraMirror(sps[1], sps[2], main = j, ppm = 2*ppm.tolerance)
      # Sys.sleep(0.5)
      add.form.sim[j] <- Spectra::compareSpectra(sps[1], sps[2], ppm = 2*ppm.tolerance)
      # cat(j, add.form.sim[j], '\n')
    }
    ann <- pubchem.bio.annotations[[i]]$annotation.table
    add.forms.all <- paste(ann$formula, ann$adducthyp)
    isotope.annotation.score <- rep(NA, length(add.forms.all))
    adduct.formula <- rep(NA, length(add.forms.all))
    for(j in 1:length(add.forms)) {
      isotope.annotation.score[which(add.forms.all == add.forms[j])] <- add.form.sim[j]
      adduct.formula[which(add.forms.all == add.forms[j])] <- ion.forms[j]
    }
    pubchem.bio.annotations[[i]]$annotation.table$isotope.score <- isotope.annotation.score
    pubchem.bio.annotations[[i]]$annotation.table$adduct.formula <- adduct.formula
    pubchem.bio.annotations[[i]]$isotope.spectra <- spectra.pairs
  }
  
  ## next task is to add a fragment ion subformula scoring.
  ## this could be either direct spectral matching, or by a more qualitative
  ## metric, asking what proportion of fragment ions are consistent with a
  ## any subformula of the assigned formula(s)
  ## or we just allow Sirius to do this.  less validation that way.  need to figure 
  ## out how to do this one at a time, constraining formula and structures for each search
  
  ## set up parallel processing
  unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  
  
  ## cleanup functions
  constrain.entropy <- function(x, max.entropy, min.int, min.peaks, ...) {
    # x <- t(matrix(c(unlist(mz(rd.cotinine[10])), unlist(intensity(rd.cotinine[10]))), nrow = 2, byrow = TRUE)); dimnames(x)[[2]] <- c("mz", "intensity")
    # min.int <- 0.01; max.entropy <- 3
    # min.peaks <- 10
    maxint <- max(x[, 2], na.rm = TRUE)
    x[, 2] <- (100 * x[, 2] / maxint)
    x <- x[which(x[,2] > 100*min.int), , drop = FALSE]
    x <- x[order(x[,2], decreasing = TRUE),]
    x <- matrix(x, ncol = 2)
    # nrow(x)
    ent <- msentropy::calculate_spectral_entropy(x)
    if(nrow(x) > min.peaks) {
      if(ent > max.entropy) {
        use <- nrow(x)
        while(ent > max.entropy & use > min.peaks) {
          ent <- msentropy::calculate_spectral_entropy(matrix(x[1:use, , drop = FALSE], ncol = 2))
          use <- use - 1
        }
      } else {
        use <- min.peaks
      }
      
      x <- x[1:(use),]
      x <- x[order(x[,1]),]
    }
    dimnames(x)[[2]] <- c("mz", "intensity")
    x
  }
  normalize.intensity <- function(x, ...) {
    maxint <- max(x[, 2], na.rm = TRUE)
    x[, 2] <- 100 * x[, 2] / maxint
    x
  }
  reorder.by.mz <- function(x, ...) {
    x <- x[order(x[,1]),, drop = FALSE]
  }
  
  if(any(ls() == 'search.spectra')) rm(search.spectra)
  
  if(!is.null(ramclustObj$mse.spectra)) {
    
    for(i in 1:length(ramclustObj$mse.spectra)) {
      cmpd <- ramclustObj$mse.spectra$cmpd[i]
      cmpd.ind <- as.numeric(gsub("C", "", cmpd))
      sp <- ramclustObj$mse.spectra[i]
      prec.mzs <- unique(pubchem.bio.annotations[[cmpd.ind]]$annotation.table$adductmz)
      for(j in 1:length(prec.mzs)) {
        sp.tmp <- ramclustObj$mse.spectra[i]
        sp.tmp$precursorMz <- prec.mzs[j]
        sp.tmp <- filterPrecursorPeaks(sp.tmp)
        sp.tmp <- applyProcessing(sp.tmp)
        if(!any(ls() == ("search.spectra"))) {
          search.spectra <- sp.tmp
        } else {
          search.spectra <- c(search.spectra, sp.tmp)
        }
      }
    }
  }
  
  if(!is.null(ramclustObj$dda.spectra)) {
    ### add dda spectra to search.spectra list, if available. 
    ### if search.spectra does not exist, create it.  
    for(i in 1:length(ramclustObj$dda.spectra)) {
      cmpd <- ramclustObj$cmpd[i]
      sp <- ramclustObj$dda.spectra[[i]]
      if(length(sp) == 0) next
      feat.names <- sapply(1:length(sp), FUN = function(x) sp[[x]]$xcms.name)
      u.feat.names <- unique(feat.names)
      prec.mzs <- sapply(1:length(sp), FUN = function(x) sp[[x]]$precursor.mz)
      for(j in 1:length(u.feat.names)) {
        do <- which (feat.names == u.feat.names[j])
        sp.sub <- sp[do]
        spd <- DataFrame(id = feat.names[do], precursorMz = prec.mzs[do], cmpd = rep(cmpd, length(do)))
        spd$mz <- lapply(1:length(sp.sub), FUN = function(x) sp.sub[[x]]$spectrum$mz)
        spd$intensity <- lapply(1:length(sp.sub), FUN = function(x) unlist(sp.sub[[x]]$spectrum$int))
        sps <- Spectra::Spectra(spd)
        
        if(length(sps) > 1) {
          stop()
          sps <- combinePeaks(
            sps,
            intensityFun = median,
            mzFun = median,
            weighted = FALSE,
            tolerance = 0,
            ppm = 10,
            timeDomain = FALSE,
            peaks = "intersect",
            f = sps$feature_id,
            minProp = 0.5,
          )
          
          sps <- applyProcessing(sps)
        }
        if(!any(ls() == ("search.spectra"))) {
          search.spectra <- sps
        } else {
          search.spectra <- c(search.spectra, sps)
        }
      }
    }
  }
  
  search.spectra <- Spectra::addProcessing(search.spectra, FUN = constrain.entropy, min.int = 0.02, max.entropy = 3, min.peaks = 10)
  search.spectra <- Spectra::addProcessing(search.spectra, FUN = normalize.intensity)
  search.spectra <- Spectra::addProcessing(search.spectra, FUN = reorder.by.mz)
  
  processingChunkSize(search.spectra) <- processingChunkSize
  if(!is.null(library.level1)) {
    #### get search results against in-house library
  }
  
  if(!is.null(library.level2a)) {
    unregister_dopar()
    snowparam <- SnowParam(workers = n.threads, type = "SOCK")
    register(snowparam)
    
    processingChunkSize(library.level2a) <- processingChunkSize
    
    library.level2a <- Spectra::addProcessing(library.level2a, FUN = constrain.entropy, min.int = 0.02, max.entropy = 3, min.peaks = 10)
    library.level2a <- Spectra::addProcessing(library.level2a, FUN = normalize.intensity)
    library.level2a <- Spectra::addProcessing(library.level2a, FUN = reorder.by.mz)
    
    ent.prm <- CompareSpectraParam(MAPFUN = joinPeaksNone, FUN = msentropy::msentropy_similarity,
                                   ms2_tolerance_in_ppm = 5, 
                                   ms2_tolerance_in_da = -1, THRESHFUN = function(x) which(x >= 0.5))
    mtch.ent <- matchSpectra(search.spectra, library.level2a, param = ent.prm, BPPARAM = snowparam)
    mtch.data <- matchedData(mtch.ent)
    q.cmpds <- mtch.data$cmpd
    for(i in 1:length(ramclustObj$cmpd)) {
      pubchem.bio.annotations[[i]]$annotation.table$library.level2a.score <- NA
      cmpd <- ramclustObj$cmpd[i]
      if(!(cmpd %in% q.cmpds)) next
      # stop()
      # lib.ind <- q.inds[which(q.cmpds == cmpd)]
      mtch.sub <- data.frame(mtch.data[which(q.cmpds == cmpd),])
      mtch.sub <- data.frame(mtch.sub[which(!is.na(mtch.sub$score)),])
      if(nrow(mtch.sub) == 0) next
      # stop()
      lib.cid <- mtch.sub$target_cid
      lib.name <-  mtch.sub$target_name
      lib.adducts <- mtch.sub$target_adduct
      query.precursorMz <- mtch.sub$precursorMz
      target.precursorMz <- mtch.sub$target_precursorMz
      target.adduct <- mtch.sub$target_adduct
      # mtch.sub$query_precursorMz <- query.precursorMz
      # mtch.sub$target_precursorMz<- target.precursorMz
      # mtch.sub$target_cid <- lib.cid
      # mtch.sub$target_adduct <- target.adduct
      pc.bio.cids <- pubchem.bio.annotations[[i]]$annotation.table$cid
      pc.bio.short.inchi <- pubchem.bio.annotations[[i]]$annotation.table$inchikey.first.block
      pc.bio.full.inchi <- pubchem.bio.annotations[[i]]$annotation.table$inchikey
      unique.lib.cids <- unique(mtch.sub$target_cid)
      unique.short.inchikeys <- unique(mtch.sub$target_inchikey.first.block)
      cid.match <- any(pubchem.bio.annotations[[i]]$annotation.table$cid %in% unique.lib.cids)
      inchi.match <- any(pubchem.bio.annotations[[1]]$annotation.table$inchikey.first.block %in% unique.short.inchikeys)
      if(cid.match | inchi.match) {
        # stop()
        for(j in 1:length(unique.short.inchikeys)) {
          ## inchi first block matches get 0.95x match score. CID match gets full match score
          mtch.sub2 <- mtch.sub[which(mtch.sub$target_inchikey.first.block == unique.short.inchikeys[j]),]
          if(nrow(mtch.sub2) == 0) next
          mtch.sub2 <- mtch.sub2[which.max(mtch.sub2$score),]
          add.msms.score <- which(
            ((abs(pubchem.bio.annotations[[i]]$annotation.table$adductmz - mtch.sub2$precursorMz[1])-pubchem.bio.annotations[[i]]$annotation.table$adductmz) < 0.1) & 
              (pubchem.bio.annotations[[i]]$annotation.table$inchikey.first.block == mtch.sub2$target_inchikey.first.block)
          )
          if(length(add.msms.score) > 0) {
            pubchem.bio.annotations[[i]]$annotation.table[add.msms.score, "library.level2a.score"] <- 0.95*mtch.sub2$score
          }
        }
        for(j in 1:length(unique.lib.cids)) {
          mtch.sub2 <- mtch.sub[which(mtch.sub$target_cid == unique.lib.cids[j]),]
          if(nrow(mtch.sub2) == 0) next
          mtch.sub2 <- mtch.sub2[which.max(mtch.sub2$score),]
          add.msms.score <- which(
            (pubchem.bio.annotations[[i]]$annotation.table$adductmz == mtch.sub2$precursorMz[1]) & 
              (pubchem.bio.annotations[[i]]$annotation.table$cid == mtch.sub2$target_cid)
          )
          if(length(add.msms.score) > 0) {
            pubchem.bio.annotations[[i]]$annotation.table[add.msms.score, "library.level2a.score"] <- mtch.sub2$score
          }
        }
      }
    }
    
  }
  
  if(!is.null(library.level2b)) {
    
    unregister_dopar()
    snowparam <- SnowParam(workers = n.threads, type = "SOCK")
    register(snowparam)
    
    processingChunkSize(library.level2b) <- processingChunkSize
    
    library.level2b <- Spectra::addProcessing(library.level2b, FUN = constrain.entropy, min.int = 0.02, max.entropy = 3, min.peaks = 10)
    library.level2b <- Spectra::addProcessing(library.level2b, FUN = normalize.intensity)
    library.level2b <- Spectra::addProcessing(library.level2b, FUN = reorder.by.mz)
    
    ent.prm <- CompareSpectraParam(MAPFUN = joinPeaksNone, FUN = msentropy::msentropy_similarity,
                                   ms2_tolerance_in_ppm = 5, 
                                   ms2_tolerance_in_da = -1, THRESHFUN = function(x) which(x >= 0.1))
    mtch.ent <- matchSpectra(search.spectra, library.level2b, param = ent.prm, BPPARAM = snowparam)
    mtch.data <- matchedData(mtch.ent)
    q.cmpds <- mtch.data$cmpd
    for(i in 1:length(ramclustObj$cmpd)) {
      pubchem.bio.annotations[[i]]$annotation.table$library.level2b.score <- NA
      cmpd <- ramclustObj$cmpd[i]
      if(!(cmpd %in% q.cmpds)) next
      # stop()
      # lib.ind <- q.inds[which(q.cmpds == cmpd)]
      mtch.sub <- data.frame(mtch.data[which(q.cmpds == cmpd),])
      mtch.sub <- data.frame(mtch.sub[which(!is.na(mtch.sub$score)),])
      if(nrow(mtch.sub) == 0) next
      # stop()
      lib.cid <- mtch.sub$target_cid
      lib.name <-  mtch.sub$target_name
      lib.adducts <- mtch.sub$target_adduct
      query.precursorMz <- mtch.sub$precursorMz
      target.precursorMz <- mtch.sub$target_precursorMz
      target.adduct <- mtch.sub$target_adduct
      # mtch.sub$query_precursorMz <- query.precursorMz
      # mtch.sub$target_precursorMz<- target.precursorMz
      # mtch.sub$target_cid <- lib.cid
      # mtch.sub$target_adduct <- target.adduct
      pc.bio.cids <- pubchem.bio.annotations[[i]]$annotation.table$cid
      pc.bio.short.inchi <- pubchem.bio.annotations[[i]]$annotation.table$inchikey.first.block
      pc.bio.full.inchi <- pubchem.bio.annotations[[i]]$annotation.table$inchikey
      unique.lib.cids <- unique(mtch.sub$target_cid)
      unique.short.inchikeys <- unique(mtch.sub$target_inchikey.first.block)
      cid.match <- any(pubchem.bio.annotations[[i]]$annotation.table$cid %in% unique.lib.cids)
      inchi.match <- any(pubchem.bio.annotations[[1]]$annotation.table$inchikey.first.block %in% unique.short.inchikeys)
      if(cid.match | inchi.match) {
        # stop()
        for(j in 1:length(unique.short.inchikeys)) {
          ## inchi first block matches get 0.95x match score. CID match gets full match score
          mtch.sub2 <- mtch.sub[which(mtch.sub$target_inchikey.first.block == unique.short.inchikeys[j]),]
          if(nrow(mtch.sub2) == 0) next
          mtch.sub2 <- mtch.sub2[which.max(mtch.sub2$score),]
          add.msms.score <- which(
            ((abs(pubchem.bio.annotations[[i]]$annotation.table$adductmz - mtch.sub2$precursorMz[1])-pubchem.bio.annotations[[i]]$annotation.table$adductmz) < 0.1) & 
              (pubchem.bio.annotations[[i]]$annotation.table$inchikey.first.block == mtch.sub2$target_inchikey.first.block)
          )
          if(length(add.msms.score) > 0) {
            pubchem.bio.annotations[[i]]$annotation.table[add.msms.score, "library.level2b.score"] <- 0.95*mtch.sub2$score
          }
        }
        for(j in 1:length(unique.lib.cids)) {
          mtch.sub2 <- mtch.sub[which(mtch.sub$target_cid == unique.lib.cids[j]),]
          if(nrow(mtch.sub2) == 0) next
          mtch.sub2 <- mtch.sub2[which.max(mtch.sub2$score),]
          add.msms.score <- which(
            (pubchem.bio.annotations[[i]]$annotation.table$adductmz == mtch.sub2$precursorMz[1]) & 
              (pubchem.bio.annotations[[i]]$annotation.table$cid == mtch.sub2$target_cid)
          )
          if(length(add.msms.score) > 0) {
            pubchem.bio.annotations[[i]]$annotation.table[add.msms.score, "library.level2b.score"] <- mtch.sub2$score
          }
        }
      }
    }
    
  }
  
  annotation.summary <- list()

  
  for(i in 1:length(pubchem.bio.annotations)) {
    pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.3 <- {
      round((pubchem.bio.annotations[[i]]$annotation.table$total_score * pubchem.bio.annotations[[i]]$annotation.table$ppm.score * pubchem.bio.annotations[[i]]$annotation.table$isotope.score)^(1/3), 4)
    }
    pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2b <- {
      round((pubchem.bio.annotations[[i]]$annotation.table$total_score * pubchem.bio.annotations[[i]]$annotation.table$ppm.score * pubchem.bio.annotations[[i]]$annotation.table$isotope.score
       * pubchem.bio.annotations[[i]]$annotation.table$library.level2b.score )^(1/4), 4)
    }
    pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2a <- {
      round((pubchem.bio.annotations[[i]]$annotation.table$total_score * pubchem.bio.annotations[[i]]$annotation.table$ppm.score * pubchem.bio.annotations[[i]]$annotation.table$isotope.score
       * pubchem.bio.annotations[[i]]$annotation.table$library.level2a.score )^(1/4), 4)
    }
    # pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.3.taxon <- {
    #   round((pubchem.bio.annotations[[i]]$annotation.table$total_score * pubchem.bio.annotations[[i]]$annotation.table$ppm.score * pubchem.bio.annotations[[i]]$annotation.table$isotope.score
    #    * pubchem.bio.annotations[[i]]$annotation.table$taxonomy.lca.similarity.aggregate)^(1/4), 4)
    # }
    # pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2b.taxon <- {
    #   round((pubchem.bio.annotations[[i]]$annotation.table$total_score * pubchem.bio.annotations[[i]]$annotation.table$ppm.score * pubchem.bio.annotations[[i]]$annotation.table$isotope.score
    #    * pubchem.bio.annotations[[i]]$annotation.table$library.level2b.score * pubchem.bio.annotations[[i]]$annotation.table$taxonomy.lca.similarity.aggregate)^(1/5), 4)
    # }
    # pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2a.taxon <- {
    #   round((pubchem.bio.annotations[[i]]$annotation.table$total_score * pubchem.bio.annotations[[i]]$annotation.table$ppm.score * pubchem.bio.annotations[[i]]$annotation.table$isotope.score
    #    * pubchem.bio.annotations[[i]]$annotation.table$library.level2a.score * pubchem.bio.annotations[[i]]$annotation.table$taxonomy.lca.similarity.aggregate)^(1/5), 4)
    # }

    # pubchem.bio.annotations[[i]]$annotation.table$annotation.score <- NULL
    # pubchem.bio.annotations[[i]]$annotation.table$taxonomy.annotation.score <- NULL
    
    # msms.detected <- !all(is.na(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2b)) | !all(is.na(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2a))
    # 
    # if(msms.detected) stop()
    
    pubchem.bio.annotations[[i]]$annotation.table$taxonomy.score <- pubchem.bio.annotations[[i]]$annotation.table$taxonomy.lca.similarity.aggregate
    
    pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$pmid.ct, decreasing = TRUE),]
    pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.3, decreasing = TRUE),] 
    pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2b, decreasing = TRUE),]
    pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2a, decreasing = TRUE),]
    pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$taxonomy.score, decreasing = TRUE),]
    
    # pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.3.taxon, decreasing = TRUE),] 
    # pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2b.taxon, decreasing = TRUE),]
    # pubchem.bio.annotations[[i]]$annotation.table <- pubchem.bio.annotations[[i]]$annotation.table[order(pubchem.bio.annotations[[i]]$annotation.table$ann.score.lev.2a.taxon, decreasing = TRUE),]
    
    use <- 1
    if(!all(is.na(pubchem.bio.annotations[[i]]$annotation.table$library.level2a.score))) {
      if(any(pubchem.bio.annotations[[i]]$annotation.table$library.level2a.score > 0.75, na.rm = TRUE)) {
        use <- which.max(pubchem.bio.annotations[[i]]$annotation.table$library.level2a.score)
      }
    }

    
    pubchem.bio.annotations[[i]]$annotation.table$assigned <- rep(FALSE, nrow(pubchem.bio.annotations[[i]]$annotation.table))
    pubchem.bio.annotations[[i]]$annotation.table$assigned[use] <- TRUE

    keep.columns <- c('adductmz',    'adducthyp', 'neutral_mass', 'total_score',
                      'cid', 'name', 'formula', 'monoisotopic.mass', 'inchikey',
                      'lca', 'pmid.ct', 'pathway.ct', 'taxonomy.ct', 'XLogP',
                      'ppm.score', 'isotope.score',
                      'library.level2a.score', 'library.level2b.score', 'ann.score.lev.3', 
                      'ann.score.lev.2b', 'ann.score.lev.2a', 'taxonomy.score')
    
    annotation.summary[[i]] <- data.frame(
      cmpd = ramclustObj$cmpd[i],
      rt = ramclustObj$clrt[i], 
      med.int = median(ramclustObj$SpecAbund[,i]),
      pubchem.bio.annotations[[i]]$annotation.table[use,keep.columns],
      link = paste0('https://pubchem.ncbi.nlm.nih.gov/compound/', pubchem.bio.annotations[[i]]$annotation.table[use, 'cid'])
    )
  }
  

  
  
  annotation.summary.df <- do.call(rbind.data.frame, annotation.summary)
  ## assign annotation confidence levels and threshold
  annotation.summary.df$annotation.confidence <- '5'
  use.3 <- which(annotation.summary.df$ann.score.lev.3 > 0.7)
  annotation.summary.df$annotation.confidence[use.3] <- '3'
  use.2b <- which(annotation.summary.df$ann.score.lev.2b > 0.5)
  annotation.summary.df$annotation.confidence[use.2b] <- '2b'
  use.2a <- which(annotation.summary.df$ann.score.lev.2a > 0.75)
  annotation.summary.df$annotation.confidence[use.2b] <- '2a'
  
  ramclustObj$pubchem.bio.annotations <- pubchem.bio.annotations
  ramclustObj$pubchem.bio.annotation.summary <- annotation.summary.df
  return(ramclustObj)
}

# use <- which(annotation.summary.df$rt > 40)
# rt <- annotation.summary.df$rt[use]
# logp <- annotation.summary.df$XLogP[use]
# plot(rt, logp, cex = annotation.summary.df$ann.score.lev.3.taxon, pch = 19)

# use <- which(annotation.summary.df$library.level2b.score > 0.30)
# annotation.summary.df[use,]
# plot(annotation.summary.df[use, "library.level2a.score"], annotation.summary.df[use, "library.level2b.score"])
