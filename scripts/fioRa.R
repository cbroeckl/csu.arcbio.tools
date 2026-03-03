remotes::install_github("janlisec/fioRa")
fioRa::run_app()  ## didn't work, prompted to install reticulate package
install.packages("reticulate")  
fioRa::run_app() ## didn't work, No valid 'default_path' provided and no valid reticulate::miniconda_path
