library(tidyverse)

myGM <- read_tsv("./Data/myGM.txt") %>%
  rename(Taxa = Name)
myGD <- read_tsv("./Data/myGD.txt")

df_transpose <- function(df) {
  
  first_name <- colnames(df)[1]
  
  temp <-
    df %>% 
    tidyr::pivot_longer(-1) %>%
    tidyr::pivot_wider(names_from = 1, values_from = value)
  
  colnames(temp)[1] <- first_name
  temp
}

df_transpose(myGD) %>%
  left_join(y = myGM, by = "Taxa") -> temp

temp[c(197, 198, 2:196)] %>%
  write_csv(file = "EMMAnorfab_gen.csv")

norfab_emma <- read_csv("./Data/norfab_EMMA.csv")
norfab_emma %>%
  dplyr::select(-any_of(c("Violet_11", "Black_26"))) -> norfab_emma

write_csv(norfab_emma, file = "./Data/norfab_EMMA.csv")
