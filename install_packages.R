inst_pack <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg,  repo = 'https://cran.rstudio.com', dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('curl', 'httr', 'Seurat', 'anndata', 'xml2', 'tidyverse', 'dplyr', 'Matrix', 'scales', 'RCurl', 'svglite', 'patchwork', 'readr', 'ggpubr', 'DT')
inst_pack(packages)

BiocManager::install(c("SingleCellExperiment", "zellkonverter", "limma", "UCell", "scuttle", "SingleR", "celldex", "DT"))
devtools::install_github("mojaveazure/seurat-disk", upgrade = "always")