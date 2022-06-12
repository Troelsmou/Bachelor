norfab <- read_csv("./Data/norfab.csv")
ProFABA <- read_csv("./Data/ProFABA.csv")
SNPs <- read_csv("Interesting_SNPs_summarized.csv")
myGM <- read_tsv("./Data/norfab_GM.tsv")
norfab_pheno <- read_csv("./Data/norfab_pheno_DMU.csv")
ProFABA_pheno <- read_csv("./Data/ProFABA_pheno_DMU.csv")
norfab_geno <- as.matrix(read.delim("./Data/norfab_GD.tsv", row.names = 1, h = T)) %>% 
  magrittr::set_rownames(plyr::mapvalues(rownames(.), from = norfab_pheno$Name, to = norfab_pheno$GPID, warn_missing = F))
ProFABA_geno <- as.matrix(read.csv("./Data/ProFABA_GD.csv", row.names = 1, h = T)) %>% 
  magrittr::set_rownames(plyr::mapvalues(rownames(.), from = ids$old_name, to = ids$DMU_id))
ids <- read_csv("./Data/Translation.csv", col_names = F,
                show_col_types = FALSE) %>% 
  mutate(DMU_id = substr(X2, 6, nchar(X2)) %>% as.numeric()) %>% 
  magrittr::set_colnames(c('old_name', 'GPID', 'DMU_id'))
merged_pheno <- rbind(norfab_pheno, ProFABA_pheno)
rbind(ProFABA_geno, norfab_geno) -> merged_geno

myEffectPlot <- function(pheno, geno, SNP) {
  merge <- geno[,SNP]
  data.frame(SNP = as.factor(merge), DMU_id = as.numeric(names(merge))) -> merge_df
  
  pheno %>% 
   # filter(DescriptionOfTrait == trait) %>% 
    mutate(Score = as.numeric(Score)) %>% 
    dplyr::select(DMU_id, Score) %>% 
    inner_join(merge_df, by = 'DMU_id') -> gg_df#%>% 
    #mutate(TRID = as.numeric(as.factor(TRID))) -> gg_df
  
  
  gg_df %>%
    group_by(SNP) %>%
    summarise(mean = mean(Score)) -> means
  print(means)
  myFit <- lm(Score ~ SNP, data = gg_df)
  R_squared <- broom::glance(myFit)[[1]] %>%
    round(digits = 5)
  
  gg_df %>%
    ggplot(aes(x = SNP, y = Score)) +
    #geom_boxplot(fill = "grey") +
    geom_jitter(height = 0) +
    geom_point(mapping = aes(y = mean, x = SNP), data = means, shape = "_", color = "red", size = 20) +
    theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 20), legend.position = "none") +
    geom_text(mapping = aes(x = -Inf, y = Inf, hjust = -2, vjust = 1.5, label = paste0("R2 = ",R_squared), color = "red")) +
    xlab("Genotype") +
    ylab("Hundred seed weight (g)") +
    labs(title = "Hundred seed weight effect plot") -> p1
  
  return(p1)  
}

myMixedModel <- function(df) {
  df %>%
    group_by(DescriptionOfTrait) %>%
    summarise(n = length(unique(TRID))) %>%
    filter(n > 1) -> real
  
  
  genetic <- c()
  enviromental <- c()
  genetic_enviromental_interaction <- c()
  variation <- c()
  residual <- c()
  for (i in real$DescriptionOfTrait) {
    print(i)
    
    df %>%
      filter(DescriptionOfTrait == i) %>%
      mutate(Name = as.factor(Name), DescriptionOfTrait = as.factor(DescriptionOfTrait), TRID = as.factor(TRID),
             replicate = as.factor(replicate)) -> y
    
    myFit <- lme4::lmer(Score ~ (1|Name) + (1|TRID) + (1|Name:TRID) + (1|replicate:TRID), data = y)
    as.data.frame(lme4::VarCorr(myFit))$sdcor -> x
    
    genetic_enviromental_interaction <- c(genetic_enviromental_interaction, x[1])
    genetic <- c(genetic, x[2])
    variation <- c(variation, x[3])
    enviromental <- c(enviromental, x[4])
    residual <- c(residual, x[5])
  }
  tibble(Trait = real$DescriptionOfTrait,
         Genetic = genetic,
         Enviromental = enviromental,
         'Genetic by Enviromental \n Interaction' = genetic_enviromental_interaction,
         Replicate = variation,
         Residual = residual) %>%
    return()
}

norfab %>%
  filter(DescriptionOfTrait == "Plant Height") %>%
  group_by(TRID) %>%
  summarise(mean = mean(na.omit(Score))) -> means
  
norfab %>%
  filter(DescriptionOfTrait == "Plant Height", TRID == 25) %>%
  ggplot(aes(x = Score)) +
  geom_histogram(bins = 50, fill = "dark blue", color = "black") +
  xlab("Plant Height (cm)") +
  #facet_wrap(~TRID) +
  #geom_vline(aes(xintercept = mean), data = means) +
  theme(plot.title = element_text(face = "bold.italic", color = "#000000", hjust = 0.5, size = 17)) +
  labs(title = "Plant Height Outlier")

myGM %>%
  filter(SNP %in% SNPs$SNP) %>%
  ggplot(mapping = aes(x = Position, y = Chromosome, color = as.factor(Chromosome))) +
  geom_point() +
#  geom_jitter(height = 0.3, width = NULL) +
  guides(color = "none") +
  scale_y_continuous(breaks = 1:7, minor_breaks = NULL) +
  labs(title = "SNPs chromosomal distribution") +
  theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 20))

# x = c(2, 17, 88)
x = 17
myEffectPlot(norfab_pheno, norfab_geno, SNPs$SNP[x], "TGW")

mf <- myMixedModel(norfab_pheno)
mf <- reshape2::melt(mf, id.vars = "Trait")

colnames(mf)[2] <- "Effect"

ggplot(mf, mapping = aes(x = Trait, y = value, fill = Effect)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Variance explaned (%)") +
  labs(title = "Norfab mixed model variance breakdown per trait") +
  theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 20), axis.text.x = element_text(angle = 90))

mf <- myMixedModel(ProFABA_pheno)
mf %>%
  mutate(Trait = ifelse(Trait == "Aerial tillers number per plant (tillering higher than basal)", "Aerial tillers", Trait)) %>%
  mutate(Trait = ifelse(Trait == "Rust disease score (percentage of covered leaf area)", "Rust disease score", Trait)) -> mf
mf <- reshape2::melt(mf, id.vars = "Trait")

colnames(mf)[2] <- "Effect"

ggplot(mf, mapping = aes(x = Trait, y = value, fill = Effect)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Variance explained (%)") +
  labs(title = "ProFaba mixed model variance breakdown per trait") +
  theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 20), axis.text.x = element_text(angle = 90))

traits <- c("End of flowering", "Flowering period duration", "Flowering time", "Chocolate spot disease score",
            "Rust disease score (percentage of covered leaf area)")

ProFABA_pheno %>%
  filter(DescriptionOfTrait %in% traits) %>%
  mutate(DescriptionOfTrait = ifelse(DescriptionOfTrait == traits[5], "Rust", DescriptionOfTrait)) %>%
  mutate(DescriptionOfTrait = ifelse(DescriptionOfTrait == traits[4], "Chocolate spot", DescriptionOfTrait)) %>%
  group_by(Name, DescriptionOfTrait) %>%
  dplyr::summarise(mean = mean(Score, na.rm = T)) %>%
  pivot_wider(id_cols = Name, names_from = DescriptionOfTrait, values_from = mean) %>%
  corrgram(lower.panel = panel.pts, upper.panel = panel.conf)

traits <- c("Seed Width", "Seed Length", "Seed Area", "TGW", "Number of ovules")

norfab_pheno %>%
  filter(DescriptionOfTrait %in% traits) %>%
  group_by(Name, DescriptionOfTrait) %>%
  dplyr::summarise(mean = mean(Score, na.rm = T)) %>%
  pivot_wider(id_cols = Name, names_from = DescriptionOfTrait, values_from = mean) %>%
  corrgram(lower.panel = panel.pts, upper.panel = panel.conf)

#===============================================================================
# Data behandling
#===============================================================================
results <- read_csv("./Mixed/norfab_Blues_Results.csv")
results %>%
  dplyr::select(-c(Chromosome, Position)) %>%
  left_join(myGM, by = "SNP") -> results

results %>%
  filter(Chromosome == 1) %>%
  mutate(effect = str_split(effect, pattern = ",")) %>%
  mutate(effect = sapply(effect, "[[", 1)) %>%
  mutate(P.value = as.numeric(effect)) -> real_ps

results %>%
  filter(Chromosome != 1) %>%
  bind_rows(real_ps) %>%
  dplyr::select(SNP, P.value, Trait, Chromosome, Position) -> results

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_Gapit_Blink_Results.csv")

temp %>%
  filter(is.na(Chromosome)) %>%
  dplyr::select(Trait, `SNP,Chromosome,Position ,P.value,maf,nobs,Rsquare.of.Model.without.SNP,Rsquare.of.Model.with.SNP,FDR_Adjusted_P-values,effect`) -> broken

colnames(broken)[2] <- "broken"

broken$broken %>%
  str_split(pattern = ",") -> info

broken %>%
  mutate(SNP = sapply(info, "[[", 1)) %>%
  mutate(Chromosome = as.numeric(sapply(info, "[[", 2))) %>%
  mutate(Position = as.numeric(sapply(info, "[[", 3))) %>%
  mutate(P.value = as.numeric(sapply(info, "[[", 4))) %>%
  dplyr::select(-broken) -> real_info

temp %>%
  filter(!(is.na(Chromosome))) %>%
  bind_rows(real_info) %>%
  dplyr::select(SNP, Chromosome, Position, P.value, Trait) -> temp


temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_Gapit_FarmCPU_Results.csv")

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_EMMAX_Results.csv")

temp %>%
  rename(Chomosome = chromosomes, Position = positions, P.value = scores) %>%
  filter(mafs > 0.05) %>%
  dplyr::select(Chomosome, Position, P.value, Trait, Trial) -> temp

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_AVG_cat_Results.csv")

temp %>%
  filter(maf > 0.05) %>%
  dplyr::select(-Position) %>%
  left_join(myGM, by = "SNP") -> temp

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_Gapit_seed_Results.csv")

colnames(temp)[1] <- "broken"
temp$broken %>%
  str_split(pattern = ",") -> info

temp %>%
  mutate(SNP = sapply(info, "[[", 1)) %>%
  mutate(Chromosome = as.numeric(sapply(info, "[[", 2))) %>%
  mutate(Position = as.numeric(sapply(info, "[[", 3))) %>%
  mutate(P.value = as.numeric(sapply(info, "[[", 4))) %>%
  dplyr::select(-broken) -> real_info

real_info$Trait -> traits

traits %>%
  str_split(pattern = "_") -> traits

real_info %>%
  mutate(Trait = sapply(traits, "[[", 1)) %>%
  mutate(Trial = as.numeric(sapply(traits, "[[", 2))) -> real_info

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_AVG_cat_Results.csv")

temp %>%
  dplyr::select(-maf, -nobs) -> temp

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_Blues_Results.csv")

temp <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_EMMAX_Results.csv")

results <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFaba_Gapit_Results.csv")

#===============================================================================
# seed
#===============================================================================

seed <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_seed.csv") %>%
  dplyr::select(-Trial)

avg_seed_ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_avg.csv") %>%
  filter(Trait %like% "testa")

avg_seed_norfab <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_avg.csv") %>%
  dplyr::select(-maf, -nobs)

ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_gapit.csv") %>%
  filter(Trait %like% "testa")

emma_norfab <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_norfab.csv") %>%
  filter(Trait %in% avg_seed_norfab$Trait)

emma_ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_ProFABA.csv") %>%
  filter(Trait %like% "testa") %>%
  dplyr::select(-Trial)

ProFABA <- bind_rows(avg_seed_ProFABA, ProFABA_seed) %>%
  bind_rows(emma_ProFABA)

ProFABA %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) -> ProFABA

norfab <- bind_rows(seed, avg_seed_norfab) %>%
  bind_rows(emma_norfab)

bind_rows(ProFABA, norfab) -> seed

seed %>%
  mutate(Trait = ifelse(Trait == "White_Gray", "White-Gray", Trait)) -> seed

#===============================================================================
# Seed size
#===============================================================================

seed_size_vec <- c("TGW", "Seed Length", "Seed Width", "Seed Area", "Number of ovules")


blink <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_blink.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% seed_size_vec)

FarmCPU <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_FarmCPU.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% seed_size_vec)


ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_gapit.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait == "Hundred seed weight")

emma_norfab <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_norfab.csv") %>%
  filter(Trait %in% seed_size_vec) %>%
  select(-Trial)

emma_ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_ProFABA.csv") %>%
  filter(Trait == "Hundred seed weight") %>%
  dplyr::select(-Trial)

norfab_blue <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_blues.csv") %>%
  filter(Trait %in% seed_size_vec)

ProFABA_blue <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_blues.csv") %>%
  filter(Trait == "Hundred seed weight")

bind_rows(blink, FarmCPU) %>%
  bind_rows(emma_norfab) %>%
  bind_rows(emma_ProFABA) %>%
  bind_rows(ProFABA) %>%
  bind_rows(norfab_blue) %>%
  bind_rows(ProFABA_blue) -> seed_size

#===============================================================================
# Seed yield
#===============================================================================
seed_yield_vec <- c("Number of seeds", "Seed yield", "Seed number per pod")

ProFABA_blue <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_blues.csv") %>%
  filter(Trait %in% seed_yield_vec)

emma_ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_ProFABA.csv") %>%
  filter(Trait %in% seed_yield_vec) %>%
  dplyr::select(-Trial)

ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_gapit.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% seed_yield_vec)

bind_rows(ProFABA_blue, emma_ProFABA) %>%
  bind_rows(ProFABA) -> seed_yield

#===============================================================================
# Norfab second cluster
#===============================================================================
vec_names <- c("Internode_Length", "Internode", "Maturation Date", "End Flowering", "Earliness of Flowering", "Duration of flowering",
               "Branching", "Plant Height", "Downy Mildew")

blink <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_blink.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% vec_names | Trait %like% "Botrytis")

FarmCPU <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_FarmCPU.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% vec_names | Trait %like% "Botrytis")

emma_norfab <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_norfab.csv") %>%
  filter(Trait %in% vec_names | Trait %like% "Botrytis") %>%
  select(-Trial)

norfab_blue <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/norfab_blues.csv") %>%
  filter(Trait %in% vec_names | Trait %like% "Botrytis")

bind_rows(blink, FarmCPU) %>%
  bind_rows(emma_norfab) %>%
  bind_rows(norfab_blue) -> norfab_second

#===============================================================================
# ProFABA plant height
#===============================================================================
ProFABA_height_vec <- c("Flower number per node", "Plant height at end of flowering", "Plant height at maturity",
                        "Branch number per plant")

ProFABA_blue <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_blues.csv") %>%
  filter(Trait %in% ProFABA_height_vec)

emma_ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_ProFABA.csv") %>%
  filter(Trait %in% ProFABA_height_vec) %>%
  dplyr::select(-Trial)

ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_gapit.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% ProFABA_height_vec)

bind_rows(ProFABA_blue, emma_ProFABA) %>%
  bind_rows(ProFABA) -> ProFABA_height

#===============================================================================
# ProFABA disease
#===============================================================================
ProFABA_disease_vec <- c("End of flowering", "Flowering period duration", "Flowering time",
                         "Rust disease score (percentage of covered leaf area)",
                         "Rust disease score (plant reaction-based scale)")

ProFABA_blue <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_blues.csv") %>%
  filter(Trait %in% ProFABA_disease_vec)

emma_ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/EmmaX_ProFABA.csv") %>%
  filter(Trait %in% ProFABA_disease_vec) %>%
  dplyr::select(-Trial)

ProFABA <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/6. Semester/Bachelor/Results/ProFABA_gapit.csv") %>%
  mutate(Trait = str_split(Trait, "_")) %>%
  mutate(Trait = sapply(Trait, "[[", 1)) %>%
  filter(Trait %in% ProFABA_disease_vec)

bind_rows(ProFABA_blue, emma_ProFABA) %>%
  bind_rows(ProFABA) -> ProFABA_disease

snp_database %>%
  filter(Trait %like% "seed testa" | Trait %like% "seed coat") %>%
  mutate(Trait = ifelse(Trait %like% "reen", "Green", Trait)) %>%
 # mutate(Trait = ifelse(Trait %like% "rown", "Beige", Trait)) %>%
  mutate(Trait = ifelse(Trait %like% "hite", "White", Trait)) %>%
  mutate(Trait = ifelse(Trait %like% "eige", "Beige", Trait)) %>%
  mutate(Trait = ifelse(Trait %like% "lack", "Black", Trait)) %>%
  mutate(Trait = ifelse(Trait %like% "iolet", "Violet", Trait)) %>%
  mutate(Trait = ifelse(Trait %like% "Yellow", "Yellow", Trait)) %>%
  mutate(Trait = ifelse(Trait %like% "Red", "Red", Trait)) %>%
  filter(Trait == "Beige") %>%
  myManhattan()


c("AX-416788501", "AX-416814019", "AX-416753479", "AX-416746203", "AX-416727691", "AX-416800905",
  "AX-416782664",
  "AX-416811818",
  "AX-416793049",
  "AX-416762933", "AX-181487700", "AX-416727621", "AX-181487929", "AX-416758714", "AX-416807440", "AX-416765207", "AX-416730290",
  "AX-416758221",
  "AX.416754732", "AX.416755984", "AX.416766685",
  "AX.416722001", "AX.416762521", "AX.416766241",
  "AX-416808136"
)
c(rep("Beige", 6),
  "Green",
  "Violet",
  "Black",
  rep("Seed size", 8),
  "Botrytis",
  rep("Plant height", 3),
  rep("Seed yield", 3),
  "End of flowering"
  )
snp_database %>%
  filter(Trait == "Plant height at maturity" & Chromosome == "Chr5") %>%
  group_by(SNP_ID, Position) %>%
  summarise(n = n())

merged_pheno %>%
  filter(DescriptionOfTrait == "Flowering time") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416755984", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "Botrytis Early(Cat)") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416758221", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "End Flowering") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-181195340", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "Downy Mildew") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416816675", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "Number of seeds") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416788159", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait %like% "Plant height at") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416755984", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "TGW") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416762933", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "TGW") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416727621", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "Hundred seed weight") %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416730290", "[-]", "."))

merged_pheno %>%
  filter(DescriptionOfTrait == "TGW") %>%
#  filter(Score < 10) %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace(save$SNP[23], "[-]", "."))

color_vector <- c("White_Gray", "Beige", "Green", "Violet", "Black")

norfab %>%
  filter(DescriptionOfTrait == "Seed Coat Colour") %>%
  mutate(Color = color_vector[Score]) %>%
  dplyr::select(Name, Color, TRID, GPID) %>%
  na.omit() -> norfab_color
norfab_color %>%
  group_by(Name, Color) %>%
  summarise(n = n()) %>%
  group_by(Name) %>%
  summarise(n = n()) %>%
  filter(n == 1) %>%
  pull(Name) -> homo
 
norfab_color %>%
  filter(Name %in% homo) %>%
  group_by(Name, Color, GPID) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(n = ifelse(n == 0 | is.na(n), 0, 1)) %>%
  pivot_wider(id_cols = c(Name, GPID), names_from = Color, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  pivot_longer(cols = all_of(color_vector)) %>%
  pivot_wider(id_cols = c(Name, GPID), names_from = name, values_from = value) %>%
  rename(Score = Beige, DMU_id = GPID) %>%
  myEffectPlot(geno = merged_geno, SNP = str_replace("AX-416814019", "[-]", "."))

beige <- norfab_color %>%
  filter(Name %in% homo) %>%
  group_by(Name, Color, GPID) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(n = ifelse(n == 0 | is.na(n), 0, 1)) %>%
  pivot_wider(id_cols = c(Name, GPID), names_from = Color, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  pivot_longer(cols = all_of(color_vector)) %>%
  pivot_wider(id_cols = c(Name, GPID), names_from = name, values_from = value) %>%
  rename(Score = Beige, DMU_id = GPID)

R2 <- rep(NA, length(SNPs$SNP_ID))
for(i in 1:length(SNPs$SNP_ID)) {
  merge <- merged_geno[,SNPs$SNP_ID[i]]
  data.frame(SNP = as.factor(merge), DMU_id = as.numeric(names(merge))) -> merge_df
  
  beige %>% 
    # filter(DescriptionOfTrait == trait) %>% 
    mutate(Score = as.numeric(Score)) %>% 
    dplyr::select(DMU_id, Score) %>% 
    inner_join(merge_df, by = 'DMU_id') -> gg_df#%>% 
  #mutate(TRID = as.numeric(as.factor(TRID))) -> gg_df
  
  
  gg_df %>%
    group_by(SNP) %>%
    summarise(mean = mean(Score)) -> means
  print(means)
  myFit <- lm(Score ~ SNP, data = gg_df)
  R_squared <- broom::glance(myFit)[[1]] %>%
    round(digits = 5)
  R2[i] <- R_squared
}

SNPs$R2 <- R2
SNPs %>%
  ggplot(mapping = aes(x = Position, y = R2, color = Chromosome)) +
  geom_point() +
  theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 20)) +
  ylab("Variance explained per SNP") +
  labs(title = "R-squared Beige")

merge <- merged_geno[,save$SNP]
as.data.frame(merge) -> merge_df

merge_df$DMU_id <- str_remove(rownames(merge_df), "X") %>%
  as.numeric()
rownames(merge_df) <- NULL
merged_pheno %>%
  filter(DescriptionOfTrait == "Hundred seed weight") %>%
  mutate(Score = as.numeric(Score)) %>% 
  dplyr::select(DMU_id, Score) %>% 
  inner_join(merge_df, by = 'DMU_id') -> gg_df

myFit <- lm(Score ~ ., data = gg_df)
broom::glance(myFit)[[1]] %>%
  round(digits = 5)
R2[i] <- R_squared