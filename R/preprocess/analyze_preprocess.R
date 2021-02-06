
## Run the preprocessing pipeline
setwd("/csc/mustjoki/aml_tim3/") ## FIMM
setwd("/Users/hru/Dropbox/aml_tim3/") ## Local

source("src/R/main.R")
source("src/R/preprocess/run_preprocessRNA.R")
source("src/R/preprocess/run_preprocessTCRab.R")
source("src/R/preprocess/run_preprocessRNATCRab.R")

message("Fin")
