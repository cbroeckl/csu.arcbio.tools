#' setup.arcbio.R
#'
#' just a convenience function to set up all packages needed for routine metabolomics work at CSU ARC-BIO
#' 
#' When setting up on Windows for the first time, you will need to manually install
#'  - 64 Bit Java
#'  - OpenBabel
#' 
#' When setting up on linux
#' using package manager:
#'  - libcurl-devel 
#'  - libcurl-openssl-devel
#'  - librsvg2-devel
#'  - fontconfig-devel
#'  - fribidi
#'  - fribidi-devel
#'  - freetype-devel 
#'  - libpng-devel 
#'  - libtiff-devel 
#'  - libjpeg-devel 
#'  - libwebp-devel
#'  - openbabel
#'  - openbabel-devel
#' From Source:
#'  - netcdf
#'  - zlib (required by netcdf)
#'  - hdf5 (optional requirement by netcdf)
#' 
#' 
#' @details can use "source('https://raw.githubusercontent.com/cbroeckl/csu.arcbio.tools/master/R/setup.arcbio.R.R'); setup.arcbio.R()" to run
#' @return nothing  
#' @author Corey Broeckling
#' 
#' @export 
#' 

setup.arcbio.R <- function() {
  options(install.packages.check.source = "no")
  options(install.packages.compile.from.source = "no")
  
  installed <- installed.packages()[,c(1)]
  
  
  cran.packages <- c(
    'devtools', 'BiocManager', 'rstudioapi', 'ggplot2', 'ggfortify', 'curl',
    'emmeans', 'effects', 'lme4', 'lmerTest', 'pbkrtest', 'ClassDiscovery', 
    'PCDimension', 'CHNOSZ', 'remotes', 'pubchem.bio', 'fioRa'  
  )
  cran.packages <- cran.packages[!(cran.packages %in% installed)]
  
  if(length(cran.packages) > 0) {
    for(i in 1:length(cran.packages)) {
      install.packages(cran.packages[i], dependencies = TRUE, update = TRUE, ask = FALSE, repos='https://mirror.las.iastate.edu/CRAN/', verbose = FALSE)
    }
  }
  
  bioc.packages <- c(
    'Biobase', 'BiocGenerics',
    'MsCoreUtils', 'MetaboCoreUtils', 'Spectra', 'MetaboAnnotation', 'CompoundDb', 'MetaboAnnotation',
    'MsBackendMgf', 'MsBackendMassbank', 'MsBackendMsp', 'mzR', 'xcms', 'ChemmineR', 'ChemmineOB'
  )
  
  bioc.packages <- bioc.packages[!(bioc.packages %in% installed)]
  if(length(bioc.packages) > 0) {
    for(i in 1:length(bioc.packages)) {
      BiocManager::install(bioc.packages[i], dependencies = TRUE, update = TRUE, ask = FALSE, verbose = FALSE)
    }
  }
  
  ## from GITHUB
  remotes::install_github('cbroeckl/RAMClustR', dependencies = TRUE, upgrade = 'always')
  remotes::install_github('cbroeckl/csu.pmf.tools', dependencies = TRUE, upgrade = 'always')
  remotes::install_github('cbroeckl/csu.arcbio.tools', dependencies = TRUE, upgrade = 'always')

}


