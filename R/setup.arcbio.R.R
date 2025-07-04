#' setup.arcbio.R
#'
#' just a convenience function to set up all packages needed for routine metabolomics work at CSU ARC-BIO
#' @details can use "source('https://github.com/cbroeckl/csu.arcbio.tools/tree/main/R/setup.arcbio.R.R'); setup.arcbio.R()" to run
#' @return nothing.  all data are saved to disk for later loading
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
    'PCDimension', 'CHNOSZ', 'remotes'  
  )
  cran.packages <- cran.packages[!(cran.packages %in% installed)]
  
  if(length(cran.packages) > 0) {
    for(i in 1:length(cran.packages)) {
      install.packages(cran.packages[i], dependencies = TRUE, update = TRUE, ask = FALSE, repos='https://mirror.las.iastate.edu/CRAN/', verbose = FALSE)
    }
  }
  
  bioc.packages <- c(
    'biobase', 'biocgenerics',
    'MsCoreUtils', 'MetaboCoreUtils', 'Spectra', 'MetaboAnnotation', 'CompoundDb',
    'MsBackendMgf', 'MsBackendMassbank', 'MsBackendMsp', 'mzR', 'xcms', 'ChemmineR', 'ChemmineOB'
  )
  
  bioc.packages <- bioc.packages[!(bioc.packages %in% installed)]
  if(length(bioc.packages) > 0) {
    for(i in 1:length(bioc.packages)) {
      BiocManager::install(bioc.packages[i], dependencies = TRUE, update = TRUE, ask = FALSE, verbose = FALSE)
    }
  }
  
  ## from GITHUB
  devtools::install_github('cbroeckl/RAMClustR', dependencies = TRUE, upgrade = 'always')
  devtools::install_github('cbroeckl/csu.pmf.tools', dependencies = TRUE, upgrade = 'always')
  devtools::install_github('cbroeckl/csu.arcbio.tools', dependencies = TRUE, upgrade = 'always')
  devtools::install_github('cbroeckl/pubchem.bio', dependencies = TRUE, upgrade = 'always')
  
}


