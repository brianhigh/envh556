# Install all of the packages used in this repository.

mylib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(mylib)) dir.create(mylib, recursive = T, showWarnings = F)

if (!suppressWarnings(require("pacman"))) install.packages("pacman")

pacman::p_load(tinytex)
if (!dir.exists(tinytex_root())) tinytex::install_tinytex()

pacman::p_load(BiocManager)
if (!suppressWarnings(require("limma"))) BiocManager::install("limma")

pacman::p_load(plyr, reshape2)
pacman::p_load(knitr, readr, dplyr, tidyr, ggplot2, formatR)
pacman::p_load(broom, tidyverse, rgdal, sp, foreign, downloader)
pacman::p_load(stringr, forcats, Hmisc, EnvStats, codetools)
pacman::p_load(multcomp, modelr, car, lme4, VCA, parallel)  
pacman::p_load(geoR, maps, scatterplot3d, ggmap, funModeling, scales, akima)
pacman::p_load(MKmisc)
pacman::p_load(stats4, boot, gsl)

