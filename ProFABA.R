library(tidyverse)
ProFABA <- read_csv(file = "./Data/ProFABA.csv")

replace_dates <- function(df) {
  ProFABA %>%
    filter(DescriptionOfTrait == "End of flowering" & !(is.na(Score)) & TRID %in% c("AGVG001", "AGVG002")) %>%
    mutate(Score = str_replace_all(Score, "-", "/")) -> endflo_df
  ProFABA %>%
    filter(!(DescriptionOfTrait == "End of flowering" & !(is.na(Score)) & TRID %in% c("AGVG001", "AGVG002"))) %>%
    mutate(Score = as.numeric(Score)) -> ret
  for (i in c("AGVG001", "AGVG002")) {
    endflo_df %>%
      filter(TRID == i) %>%
      mutate(Score = as.POSIXct(Score, format = "%d/%m/%Y")) %>%
      mutate(Score = as.numeric(Score)/86400) %>%
      mutate(Score = round(Score)) %>%
      mutate(Score = Score-min(Score)) %>%
      bind_rows(ret) -> ret
  }
  return(ret)
}

reformat_df <- function(df) {
  df %>%
    unite(Trait_and_Trial, c(DescriptionOfTrait, TRID)) %>%
    dplyr::select(Score, Trait_and_Trial, Name) %>% 
    na.omit() %>%
    group_by(Name, Trait_and_Trial) %>%
    summarise(Score = mean(Score, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(id_cols = Name, names_from = Trait_and_Trial, values_from = Score) %>%
    return()
}

remove_no_variance_col <- function(df) {
  df %>%
    dplyr::select(where(is.numeric)) %>%
    colnames() -> iter
  
  keep <- c()
  for (i in iter) {
    df[[i]] %>%
      var(na.rm = T) -> variance
    if (variance > 0) {
      keep <- c(keep, i)
    }
  }
  df[c("Name",keep)] %>%
    return()
}

names_translation <- read_csv(file = "./Data/Translation.csv", col_names = c("Name", "GPID"))

vec <- names_translation$Name

names(vec) <- names_translation$GPID

ProFABA <- replace_dates(ProFABA)

ProFABA %>%
  mutate(Name = vec[GPID]) -> ProFABA

myPlot <- function(df, x, y) {
  full_traits <- df$DescriptionOfTrait %>% unique()
  
  df %>%
    filter(DescriptionOfTrait == full_traits[x]) %>%
    group_by(TRID) %>%
    summarise(mean = mean(na.omit(Score))) -> means
  
  df %>%
    filter(DescriptionOfTrait == full_traits[x]) %>%
    ggplot(aes(x = Score)) +
    geom_histogram(bins = y, fill = "red", color = "black") +
    # xlim(c(0,200)) +
    # ylim(c(0, 100)) +
    facet_wrap(~TRID) +
    geom_vline(aes(xintercept = mean), data = means) +
    labs(title = full_traits[x])
}

traits <- c("End of flowering", "Flowering time", "Flower number per node","Flowering period duration","Branch number per plant",
            "Aerial tillers number per plant (tillering higher than basal)", "Plant height at end of flowering", "Plant number",
            "Plant height at maturity", "First pod position", "Pod number per node", "Pod length", "Seed number per pod",
            "Seed yield","Hundred seed weight", "Number of seeds", "Seed health status", "Rust disease score (plant reaction-based scale)",
            "Rust disease score (percentage of covered leaf area)", "Chocolate spot disease score", "Aphid resistance index",
            "Broomrape infestation index", "Ascochyta blight severity score (leaves)", "Ascochyta blight severity score (pods)",
            "Herbicide damage", "Downy mildew (Peronospora vicia)", "Virus", "Sitona resistance index",
            "Percent of bruchid-infested seeds", "Field emergence", "Delayed emergence", "Lodging incidence")

Weird_traits <- c("Leaflet shape", "Stem anthocyanin coloration", "Standard petal colour", "Wing petal colour", "Pod colour",
                  "Seed shape", "Seed testa ground colour", "Hilum color", "Germination")

ProFABA %>%
  filter(!(DescriptionOfTrait == "Hundred seed weight" & Score > 400)) %>%
  filter(!(DescriptionOfTrait == "Seed yield" & Score > 3000)) %>%
  filter(!(is.na(Score))) -> ProFABA

ProFABA %>%
  filter(DescriptionOfTrait %in% traits) %>%
  reformat_df() -> n_df



#n_df %>%
#  group_by(DescriptionOfTrait, TRID) %>%
#  summarise(u = unique(Date)) %>%
#  group_by(DescriptionOfTrait, TRID) %>%
#  summarise(n = n()) %>%
#  arrange(-n) %>%
#  filter(DescriptionOfTrait %in% traits[c(19, 21, 23, 25, 27, 18, 20, 22, 24, 26, 26)])

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[1]) %>%
  mutate(Score = ifelse(Score < 1.6, "Narrow elongated", Score)) %>%
  mutate(Score = ifelse(Score >= 1.6 & Score < 2.5, "Sub elliptic", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Sub orbicular", Score)) %>%
  mutate(Score = ifelse(Score == 4, "Mixed", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) -> store_df

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[2] & Score > 0) %>%
  mutate(Score = ifelse(Score == 1, 0, Score)) %>%
  mutate(Score = ifelse(Score == 3, 1, Score)) %>%
  reformat_df() %>%
  merge(y = n_df, all = T) -> n_df

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[3] & Score %in% c(1, 2, 5, 6)) %>%
  mutate(Score = ifelse(Score == 1, "White", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Violet", Score)) %>%
  mutate(Score = ifelse(Score == 5, "Pink", Score)) %>%
  mutate(Score = ifelse(Score == 6, "Red", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[4] & Score %in% 1:3) %>%
  mutate(Score = ifelse(Score == 1, "White", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Coloured", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Spotted", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

#ProFABA %>%
#  filter(DescriptionOfTrait == Weird_traits[5]) %>%
#  mutate(Score = ifelse(Score < 1.5, "Light yellow", Score)) %>%
#  mutate(Score = ifelse(Score < 2.5, "Dark", Score)) %>%
#  mutate(Score = ifelse(Score >= 2.5 & Score < 4, "Mixed", Score)) %>%
#  unite(Score, c(Score, DescriptionOfTrait)) %>%
#  bind_rows(store_df) -> store_df

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[6] & Score %in% 1:3) %>%
  mutate(Score = ifelse(Score == 1, "Flat", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Angular", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Round", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[7] & Score %in% 1:9) %>%
  mutate(Score = ifelse(Score == 1, "Black", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Light brown", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Dark brown", Score)) %>%
  mutate(Score = ifelse(Score == 4, "Light green", Score)) %>%
  mutate(Score = ifelse(Score == 5, "Dark green", Score)) %>%
  mutate(Score = ifelse(Score == 6, "Red", Score)) %>%
  mutate(Score = ifelse(Score == 7, "Violet", Score)) %>%
  mutate(Score = ifelse(Score == 8, "Yellow", Score)) %>%
  mutate(Score = ifelse(Score == 9, "White", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[8] & Score %in% 1:4) %>%
  mutate(Score = ifelse(Score == 1, "Black", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Grey", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Colorless", Score)) %>%
  mutate(Score = ifelse(Score == 4, "Mixed", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

ProFABA %>%
  group_by(Name, TRID) %>%
  summarise(DescriptionOfTrait = "Germination", Score = 50) %>%
  filter(!(TRID == "GOET003" & Name %in% c("Sel97-1", "75573", "ViciapaucijugaVF635", "Kassa", "Farah", "Victus")))-> temp

ProFABA %>%
  filter(DescriptionOfTrait == Weird_traits[9]) %>%
  dplyr::select(Name, TRID, DescriptionOfTrait, Score) %>%
  merge(y = temp, all = TRUE) %>%
  filter(TRID == "GOET003") %>%
  reformat_df() %>%
  merge(y = n_df, all = T) -> n_df

remove(temp)

store_df %>%
  dplyr::select(Name, Score, TRID) %>%
  na.omit() %>%
  group_by(Name, TRID, Score) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Name, TRID), names_from = Score, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) -> store_df

pivot_longer(store_df, cols = colnames(store_df)[-c(1,2)]) %>%
  rename(DescriptionOfTrait = name) %>%
  rename(Score = value) %>%
  reformat_df() %>%
  merge(y = n_df, all = T) %>%
  remove_no_variance_col() -> myY

myGD <- read_csv(file = "./Data/ProFABA_GD.csv")
myGM <- read_csv(file = "./Data/ProFABA_GM.csv")

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  SNP.MAF = 0.05,
  model = c("FarmCPU")
)
