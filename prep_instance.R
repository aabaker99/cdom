local({
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.r-project.org" 
  options(repos=r)
})
install.packages('BiocManager')
requireNamespace('BiocManager')

# TODO is this needed?
BiocManager::install('biomaRt')

install.packages('RSQLite')
install.packages('argparse')
