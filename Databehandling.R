library(tidyverse)
SNPs <- read_csv("SNPs.csv")
color <- read_csv("Imfaba_color.csv")


SNPs %>%
  mutate(Chromosome = as.factor(Chromosome)) %>%
  mutate(Trait = as.factor(Trait)) %>%
  mutate(Trial = as.factor(Trial)) -> SNPs

SNPs %>%
  group_by(Trait) %>%
  summarise(Blinks = sum(Blink), FarmCPUs = sum(FarmCPU), Total = sum(Blinks+FarmCPUs)) %>%
  ggplot(mapping = aes(y = Total, x = Trait)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

SNPs %>%
  group_by(SNP, Trait, Trial) %>%
  summarise(n = sum(Blink+FarmCPU)) %>%
  ggplot(mapping = aes(x = n)) +
  geom_histogram(bins = 4, color = "Black", fill = "Firebrick") +
  labs(title = "Distribution of times a SNP is significant", subtitle = "Sum of times a SNP is significant per Trait and Trial")

SNPs %>%
  group_by(SNP, Trait) %>%
  summarise(Blink = sum(Blink), FarmCPU = sum(FarmCPU)) %>%
  filter(Blink > 0 & FarmCPU > 0)

tibble(Values = c(312, 85, 65), Names = c("Total unique",
                                          "Significant with both methods",
                                          "Unique SNPs significant in multiple models")) %>%
  ggplot(mapping = aes(y = Values, x = Names)) +
  geom_col(fill = "Firebrick", color = "Black") +
  labs(title = "Times when a SNP was significant in both models", subtitle = "Per trait for all trials")

SNPs %>%
  group_by(SNP, Trait) %>%
  summarise(Blink = sum(Blink), FarmCPU = sum(FarmCPU)) %>%
  filter(Blink > 0 & FarmCPU > 0) %$%
  unique(SNP) -> doubl


SNPs %>%
  filter(SNP %in% doubl) %>%
  mutate(id = paste(SNP, Trait, Trial)) %>%
  dplyr::select(id, P.value, Blink, FarmCPU, Trait) %>%
  filter(id != "AX-416758221 Botrytis Fabae (Per) 22") %>%
  pivot_wider(id_cols = c(id, Trait), values_from = P.value, names_from = c(Blink, FarmCPU)) %>%
  rename("Blink_P_value" = '1_0') %>%
  rename("FarmCPU_P_value" = '0_1') %>%
  ggplot(mapping = aes(x = Blink_P_value, y = FarmCPU_P_value, fill = Trait)) +
  geom_point(aes(color = Trait, size = 2)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Correlation between FarmCPU and Blink P values", subtitle = "For SNPs significant in both models for one trait and trial, scales are log transformed")      
  
SNPs %>%
  group_by(SNP, Trait, Trial) %>%
  summarise() %>%
  ungroup() %>%
  group_by(SNP, Trait) %>%
  summarise(n = n()) %>%
  ggplot(mapping = aes(x = n, fill = Trait)) +
  geom_histogram(bins = 5, color = "Black") +
  labs(title = "Times a SNP was significant for one trait")

#===============================================================================
#Magic population
#===============================================================================
read_csv("Magic_snps.csv") -> Magic

SNPs %>%
  group_by(SNP) %>%
  summarise(Total = sum(Blink + FarmCPU)) %>%
  filter(Total > 1) %$%
  SNP -> double_snps

SNPs$SNP %>% 
  unique() -> total_snps

Magic %>%
  filter(SNP %in% double_snps)

chromosomes <- c("Chr1S", "Chr1L", "Chr2",  "Chr3",  "Chr4",  "Chr5",  "Chr6")

SNPs %>%
  mutate(Chromosome = chromosomes[Chromosome]) -> SNPs

color %>%
  mutate(Chromosome = chromosomes[Chromosome]) -> Color

findDistance <- function(df1, df2) {
  df1 %>%
    group_by(SNP, Chromosome, Position) %>%
    summarise(n = n()) -> search_df
  
  snps1 <- c()
  snps2 <- c()
  distance <- c()
  for(i in 1:length(search_df$SNP)) {
    chromosome <- search_df$Chromosome[i]
    position <- search_df$Position[i]
    find_df <- df2 %>%
      filter(Chromosome == chromosome)
    for (j in 1:length(find_df$SNP)) {
      snps1 <- c(snps1, search_df$SNP[i])
      snps2 <- c(snps2, find_df$SNP[j])
      distance <- c(distance, abs(position - find_df$Position[j]))
    }
  }
  ret_df <- tibble(SNP1 = snps1,
                   SNP2 = snps2,
                   distance = distance)
  return(ret_df)
}

findDistance(SNP, Magic) -> Distance

Distance %>%
  ggplot(mapping = aes(x = distance, fill = "firebrick", color = "black")) +
  geom_histogram(bins = 100) +
#  scale_x_log10()
  scale_x_continuous(limits = c(0, 10^6))

Distance %>%
  filter(distance <= 10^7) -> Distance

SNPs %>%
  filter(SNP %in% Distance$SNP1) %>%
  group_by(SNP, Chromosome, Position, Trait) %>%
  summarise(total = sum(Blink + FarmCPU))

find_traits <- function(df1, df2) {
  findDistance(df1, df2) %>%
    filter(distance <= 10^7) -> Distance
  Trait1 <- c()
  Trait2 <- c()
  for (i in 1:nrow(Distance)) {
    df1 %>%
      filter(SNP == Distance$SNP1[i]) %$%
      unique(Trait) -> trait_vec1
      Trait1 <- c(Trait1, str_c(trait_vec1, collapse = " | "))
    df2 %>%
      filter(SNP == Distance$SNP2[i]) %$%
      unique(Trait) -> trait_vec2
    Trait2 <- c(Trait2, str_c(trait_vec2, collapse = " | "))
  }
  tttbl <- tibble(Trait1 = Trait1, Trait2 = Trait2)
  return(tttbl)
}
write_distance <- function(df1, df2, name) {
  findDistance(df1, df2) %>%
    filter(distance <= 10^7) -> Distance
  
    find_traits(df1, df2) %>%
    bind_cols(Distance) %>%
    write_csv(file = name)
}

SNPs %>%
  mutate(Method = ifelse(Blink == 1, "Blink", "FarmCPU")) -> coded

condenseSNPs <- function(df) {

  df %>%
    group_by(SNP, Chromosome, Position) %>%
    summarise() -> out_df
  
  nrow(out_df) -> n
  
  maf_store <- rep(NA, n)
  trait_store <- rep(NA, n)
  count_store <- rep(NA, n)
  p_value_store <- rep(NA, n)
  trial_store <- rep(NA, n)
  method_store <- rep(NA, n)
  
  for(i in 1:n) {
    
    indices <- df$SNP == out_df$SNP[i]
    
    df$maf[indices] %>%
      str_flatten(collapse = " | ") -> maf_store[i]
    df$P.value[indices] %>%
      str_flatten(collapse = " | ") -> p_value_store[i]
    df$Trial[indices] %>%
      str_flatten(collapse = " | ") -> trial_store[i]
    df$Method[indices] %>%
      str_flatten(collapse = " | ") -> method_store[i]
    
    df$Trait[indices] -> traits
    traits %>% 
      unique() -> u_traits
    u_traits %>%
      str_flatten(collapse = " | ") -> trait_store[i]
    
    count <- c()
    for (j in 1:length(u_traits)) {
      traits == u_traits[j] -> temp
      count <- c(count, sum(temp))
    }
    
    count %>%
      str_flatten(collapse = " | ") -> count_store[i]
  }
  
  out_df$maf <- maf_store
  out_df$Trait <- trait_store
  out_df$Count <- count_store
  out_df$P.value <- p_value_store
  out_df$Trial <- trial_store
  out_df$Method <- method_store
  return(out_df)
}
#===============================================================================
#Effect plots
#===============================================================================

myGD <- read.table("./Data/myGD.txt", head = TRUE) %>%
  as.tibble()

interesting <- c(temp$SNP1, double_snps) %>%
  unique() %>%
  str_replace_all("AX-", "AX.")

x <- 101

SNPs %>%
  filter(SNP == interesting[x]) %$%
  unique(Trait)

trait <- "Seed Width"

helper <- myGD[[interesting[x]]] %>%
  setNames(nm = myGD$Taxa)

faba_df_filtered %>%
  dplyr::select(Score, DescriptionOfTrait, Name, TRID) %>%
  filter(DescriptionOfTrait == trait) %>%
  group_by(Name, TRID) %>%
  summarise(meanscore.y = mean(Score)) %>%
  mutate(genotype = helper[Name]) %>%
  na.omit() %>%
  ggplot(aes(x=factor(genotype), y=meanscore.y, fill=factor(TRID)) ) + 
  geom_jitter(position=position_jitter(width=0.3), aes(colour=factor(TRID)), alpha=0.2) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  facet_wrap(~TRID, scales = "free_y") +
  theme(strip.text.x = element_text(size=9, color="black", face="bold"), legend.position = "none") +
  ylab(trait) +
  theme_classic() +
  labs(title = paste(trait, interesting[x], sep = " | "))

#ggsave(paste('EffectPlotMeansTrial',SNP,trait,Sys.Date(),'.png',sep="_"), plot = p2, width = 15, height = 10, unit = 'cm')

#===============================================================================
#ProFABA prelimenary data wrangling
#===============================================================================

ProFABA <- read_csv("SNPs_ProFABA.csv")

#Takes the first dataframe with name1 vector in it, finds the matches in df2
#in name2 and returns the corresponding values from values vector in df2 in order.
get_corresponding_values <- function(df1, df2, name1, name2, values) {
  ret <- c()
  searchkey <- df1[[name1]]
  search <- df2[[name2]]
  searchreturn <- df2[[values]]
  
  for (i in 1:length(searchkey)) {
    x <- searchreturn[which(search == searchkey[i])]
    
    ret <- c(ret, x)
  }
  
  return(ret)
}

ProFABA %>%
  mutate(Trial = ifelse(Trial == "AGVG01", "AGVG001", Trial)) %>%
  mutate(Trial = ifelse(Trial == "AGVG02", "AGVG002", Trial)) -> ProFABA

ProFABA[1:8] -> ProFABA

ProFABA %>%
  mutate(Chromosome = get_corresponding_values(ProFABA, myGM, "SNP", "Name", "Chromosome")) %>%
  mutate(Chromosome = chromosomes[Chromosome]) -> ProFABA

Magic %>%
  mutate(Population = "Magic") -> Magic

SNPs %>%
  bind_rows(color) %>%
  mutate(population = "Norfab") -> SNPs


