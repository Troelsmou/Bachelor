library(tidyverse)

filelist <- list.files(pattern="MET_GWAS_results")

exp <- read_delim(file = filelist[1], delim = " ")

traits_trials = str_split(filelist, pattern = "[.]")

for (i in 2:length(filelist)) {
  
  temp <- read_delim(file = filelist[i], delim = " ")
  print(i)
  print(length(filelist))
  exp <- bind_rows(exp, temp)
}

write_csv(x = exp, file = "ProFABA_AVG_cat_Results.csv")

