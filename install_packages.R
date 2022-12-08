inst_pack <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg,  repo = 'https://cran.rstudio.com', dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('BiocManager',  'rmarkdown', 'devtools', 'curl', 'httr', 'Seurat', 'SingleCellExperiment', 'xml2', 'tidyverse', 'dplyr', 'Matrix', 'scales', 'RCurl', 'svglite')
inst_pack(packages)

devtools::install_github("mojaveazure/seurat-disk", upgrade = "always")