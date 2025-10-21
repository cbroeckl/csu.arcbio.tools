
## add CID to NIST msp-derived spectral library
library(MetaboAnnotation)
library(Spectra)

## Load spectra (on my windows desktop)
load("C:/Temp/20250828/fiora/fiora.pos.Rdata")  ## fiora predicted pos mode spectra for pubchembio compounds from 20250828. 
load("R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.pos.Rdata") ## NIST 23 pos mode spectra 

## generate a random subset for benchmarking
set.seed(123456)
nist.sub.ind <- sample(1:length(nist.hr.pos), 10000)
nist.sub <- nist.hr.pos[nist.sub.ind]

fiora.sub.ind <- sample(1:length(fiora.pos), 10000)
fiora.sub <- fiora.pos[fiora.sub.ind]

## MatchForwardReverseParam returns both forward and reverse match score
## this sets up parameters.  
csp.fr <- MatchForwardReverseParam(
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 25,
  FUN = MsCoreUtils::ndotproduct,
  requirePrecursor = TRUE,
  requirePrecursorPeak = FALSE,
  THRESHFUN = function(x) which(x >= 0.2),  ## all matches above 0.2 returned
  THRESHFUN_REVERSE = function(x) which(x >= 0.2), ## all matches above 0.2 returned
  toleranceRt = Inf,
  percentRt = 0
)

p1 <- BiocParallel::SnowParam(1)
p2 <- BiocParallel::SnowParam(2)
p3 <- BiocParallel::SnowParam(3)
p4 <- BiocParallel::SnowParam(4)
library(microbenchmark)
microbenchmark(
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p1),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p2),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p3),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p4),
  times = 3
)
## when all fit in one chunk, parallelizing slows things down
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p1) 15.28431 15.37370 15.50163 15.46310 15.61029 15.75747     3 a  
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p2) 24.26542 24.35954 26.76186 24.45366 28.01009 31.56652     3  b 
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p3) 34.00575 34.64650 34.95765 35.28725 35.43360 35.57994     3  bc
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p4) 33.77469 34.96326 38.41339 36.15183 40.73274 45.31366     3   c

processingChunkSize(nist.sub) <- 1000
library(microbenchmark)
microbenchmark(
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p1),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p2),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p3),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p4),
  times = 3
)
# with smaller chunks, there may be speed to be had with more threads
# Unit: seconds
# expr       min       lq     mean   median       uq      max neval cld
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p1)  9.867816 12.69523 13.74772 15.52265 15.68767 15.85270     3  a 
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p2) 21.319669 21.59279 22.80442 21.86590 23.54680 25.22770     3  ab
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p3) 28.109204 31.21368 32.58852 34.31816 34.82818 35.33820     3   b
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p4) 27.764921 28.54710 33.66759 29.32928 36.61892 43.90855     3   b

#  try with bigger fiora library
fiora.sub.ind <- sample(1:length(fiora.pos), 100000)
fiora.sub <- fiora.pos[fiora.sub.ind]
processingChunkSize(nist.sub) <- 2000
library(microbenchmark)
microbenchmark(
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p1),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p2),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p3),
  matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p4),
  times = 3
)

# trend similar with a bigger library to search against, 1 still faster. 
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p1) 18.91276 26.12816 32.96830 33.34355 39.99607 46.64859     3   a
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p2) 34.56181 47.14283 51.45775 59.72385 59.90572 60.08759     3   a
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p3) 41.65327 41.89187 50.93195 42.13046 55.57129 69.01212     3   a
# matchSpectra(nist.sub, fiora.sub, param = csp.fr, BPPARAM = p4) 53.24308 55.23580 55.90669 57.22851 57.23849 57.24847     3   a



register(BiocParallel::SnowParam(2))  ## tried also with 6
nist.mh <- which(nist.hr.pos$adduct == "[M+H]+")
set.seed(123456)
nist.mh.samp <- sample(nist.mh, round(length(nist.mh)/5))  ## 117306 random pos mode M+H spectra, searched against  ## fiora predicted M+H spectra
a <- Sys.time()
matches.fr <- matchSpectra(nist.hr.pos[nist.mh.samp], fiora.pos, param = csp.fr, BPPARAM = p1)
b <- Sys.time()

processingChunkSize(nist.hr.pos)
processingChunkSize(fiora.pos)



length(table(matches.fr@matches$query_idx))
tar.mtch <- 66
nist.hr.pos$inchikey[nist.sample[tar.mtch]]
nist.hr.pos$name[nist.sample[tar.mtch]]
nist.hr.pos$cid[nist.sample[tar.mtch]]
nist.hr.pos$cid.parent[nist.sample[tar.mtch]]
nist.hr.pos$precursorMz[nist.sample[tar.mtch]]
nist.hr.pos$formula[nist.sample[tar.mtch]]
fiora.pos$title[matches.fr@matches$target_idx[which(matches.fr@matches$query_idx == tar.mtch)]]
fiora.pos$formula[matches.fr@matches$target_idx[which(matches.fr@matches$query_idx == tar.mtch)]]
mtch.ind <- which(matches.fr@matches$query_idx == tar.mtch)
mtch <- data.frame(matches.fr@matches[mtch.ind,], 'cid'=fiora.pos$cid[matches.fr@matches[mtch.ind,"target_idx"]])
mtch[order(mtch$score, decreasing = TRUE),]

fragment.smiles <- unlist(fiora.pos$V1)
length(unique(fragment.smiles))

plot(enviPat::envelope(enviPat::isopattern(isotopes, "C39H76NO7P"), resolution = 80000)[[1]], type = "l", xlim = c(701.5, 701.57))
points(enviPat::envelope(enviPat::isopattern(isotopes, "C39H74NO7P"), resolution = 80000)[[1]], type = "l", col = 2)
