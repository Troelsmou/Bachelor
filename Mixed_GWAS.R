library(tidyverse)

ProFABA_blue <- read_csv()

norfab_blue <- read_csv()

norfab_blue %>%
  mutate(Line = str_remove(Line, pattern = "Name")) -> norfab_blue

ProFABA_blue %>%
  mutate(Line = str_remove(Line, pattern = "Name")) -> ProFABA_blue

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
setwd("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Mixed/ProFABA")
myY <- ProFABA_blue %>%
  as.data.frame()
myGM <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Mixed/Data/ProFABA_GM.csv")
myGD <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Mixed/Data/ProFABA_GD.csv")

myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  SNP.MAF = 0.05,
  model = c("FarmCPU")
)

setwd("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Mixed/norfab")
myY <- norfab_blue %>%
  as.data.frame()

myGM <- read_tsv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Mixed/Data/myGM.txt")
myGD <- read_tsv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Mixed/Data/myGD.txt")

myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  SNP.MAF = 0.05,
  model = c("FarmCPU")
)

