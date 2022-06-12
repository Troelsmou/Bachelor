library(lme4)
library(tidyverse)
library(lmerTest)
library(varhandle)
library(magrittr)
library(scales)
library(nlme)
norfab <- read_csv(file = "./Data/Faba.csv")

norfab <- norfab %>%
  filter(!(DescriptionOfTrait == "Plant Height" & Score > 500)) %>%
  filter(!(DescriptionOfTrait == "Number of ovules" & Score > 9)) %>%
  filter(!(DescriptionOfTrait == "End Flowering" & Score > 55)) %>%
  filter(!(DescriptionOfTrait == "Seed Coat Colour" & Score > 8))

norfab %>%
  dplyr::select(Score, Name, TRID, DescriptionOfTrait) -> norfab

unicize_rows <- function(df) {
  df <- arrange(df, Name, DescriptionOfTrait, TRID) %>%
    filter(!(is.na(Name) | is.na(DescriptionOfTrait) | is.na(TRID) | is.na(Score)))
  
  Name_i <- df$Name[1]
  TRID_i <- df$TRID[1]
  DescriptionOfTrait_i <- df$DescriptionOfTrait[1]
  replicate <- rep(NA, nrow(df))
  replicate[1] <- 1
  for (i in 1:(nrow(df)-1)) {
    
    if (Name_i == df$Name[i+1] & TRID_i == df$TRID[i+1] & DescriptionOfTrait_i == df$DescriptionOfTrait[i+1]) {
      replicate[i+1] <- replicate[i]+1
    }
    else {
      replicate[i+1] <- 1
    }
    
    Name_i <- df$Name[i+1]
    TRID_i <- df$TRID[i+1]
    DescriptionOfTrait_i <- df$DescriptionOfTrait[i+1]
  }
  tibble(Name = df$Name,
         TRID = df$TRID,
         DescriptionOfTrait = df$DescriptionOfTrait,
         Score = df$Score,
         replicate = replicate) %>%
    return()
}

norfab %>%
  unicize_rows() -> norfab

myMixedModel <- function(df, categorical_traits) {
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
  tibble(traits = real$DescriptionOfTrait,
         genetic = genetic,
         enviromental = enviromental,
         genetic_X_enviromental = genetic_enviromental_interaction,
         replicate = variation,
         residual = residual) %>%
    return()
}

norfab_cat <- c("Seed Coat Colour", "Seed Hilum Colour", "Borytis Fabae", "Botrytis Fabae", "Duration of flowering",
                "Internode_length", "Novitron harm", "Stress", "Uniformity")

variances <- myMixedModel(norfab, norfab_cat)

mdf <- reshape2::melt(variances, id.vars = "traits")

ggplot(mdf, mapping = aes(x = traits, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  labs(title = "Norfab mixed model variance breakdown per trait")

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

ProFABA <- replace_dates(ProFABA)

ProFABA %>%
  filter(!(DescriptionOfTrait == "Hundred seed weight" & Score > 400)) %>%
  filter(!(DescriptionOfTrait == "Seed yield" & Score > 3000)) %>%
  filter(!(is.na(Score))) -> ProFABA

ProFABA_cat <- c("Leaflet shape", "Stem anthocyanin coloration", "Standard petal colour", "Wing petal colour", "Pod colour",
                 "Seed shape", "Seed testa ground colour", "Hilum color", "Germination", "Delayed emergence",
                 "Field emergence", "Lodging incidence", "Rust disease score (plant reaction-based scale)",
                 "Sitona resistance index")

ProFABA %>%
  dplyr::select(Score, Name, TRID, DescriptionOfTrait) -> ProFABA

unicize_rows(ProFABA) -> ProFABA

variances <- myMixedModel(ProFABA, ProFABA_cat)

mdf <- reshape2::melt(variances, id.vars = "traits")
ggplot(mdf, mapping = aes(x = traits, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "ProFABA mixed model variance breakdown per trait")

#===============================================================================
#Categorical traits
#===============================================================================

norfab %>%
  filter(DescriptionOfTrait == "Seed Coat Colour") %>%
  mutate(Color = color_vector[Score]) %>%
  dplyr::select(Name, Color, TRID) %>%
  na.omit() %>%
  group_by(Name, TRID, Color) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Name, TRID), names_from = Color, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  pivot_longer(cols = all_of(color_vector)) %>%
  mutate(value = ifelse(value == 0, 0, 1)) -> norfab_color
norfab_color %>%
  rename(DescriptionOfTrait = name, Score = value) -> norfab_color
norfab %>%
  dplyr::select(Name, TRID, DescriptionOfTrait, Score) %>%
  filter(DescriptionOfTrait != "Seed Coat Colour") %>%
  bind_rows(norfab_color) -> norfab_mixed

#===============================================================================
#ANOVA
#===============================================================================
n.obs <- function(df, column) {
  traits <- unique(df$DescriptionOfTrait)
  save <- rep(NA, length(traits))
  
  for (i in 1:length(traits)) {
    df %>%
      filter(DescriptionOfTrait == traits[i]) -> temp
    temp[[column]] %>%
      unique() %>%
      length() -> save[i]
  }
  return(save)
}

myANOVA <- function(df) {
  
  #For storing results
  trait <- unique(df$DescriptionOfTrait)
  mean <- rep(NA, length(trait))
  p_mean <- rep(NA, length(trait))
  VarG <- rep(NA, length(trait))
  p_G <- rep(NA, length(trait))
  n_obs_G <- n.obs(df, "Name")
  VarE <- rep(NA, length(trait))
  p_E <- rep(NA, length(trait))
  n_obs_E <- n.obs(df, "TRID")
  VarGxE <- rep(NA, length(trait))
  p_GxE <- rep(NA, length(trait))
  n_obs_GxE <- rep(NA, length(trait))
  VarRep <- rep(NA, length(trait))
  p_rep <- rep(NA, length(trait))
  n_obs_Rep <- rep(NA, length(trait))
  VarRes <- rep(NA, length(trait))
  n_datapoints <- rep(NA, length(trait))
  H2 <- rep(NA, length(trait))
  
  #Making columns factors
  df %>%
    mutate(Name = as.factor(Name), TRID = as.factor(TRID), replicate = as.factor(replicate)) -> df
  
  #Computing values for each trait
  for (i in 1:length(trait)) {
    df %>%
      filter(DescriptionOfTrait == trait[i]) %>%
      na.omit() -> temp_df
    
    n_datapoints[i] <- nrow(temp_df)
    
    #Is there more than 1 environment tested?
    if (n_obs_E[i] > 1) {
      #Is the trait binary?
      if (all(unique(temp_df$Score) %in% c(0,1))) {
        main_model <- glmer(Score ~ (1|Name) + (1|TRID) + (1|Name:TRID) + (1|replicate:TRID), family='binomial',data = temp_df)
        # the residual variance in a logistic regression is fixed 
        VarRes[i] <- (pi ^ 2) / 3
        
        p_mean[i] <-  coef(summary(main_model))[4]
        
        model_withoutGxE = glmer(Score ~ (1|Name) + (1|TRID) + (1|replicate:TRID), family='binomial', data = temp_df) 
        p_GxE[i] = anova(main_model,model_withoutGxE)$`Pr(>Chisq)`[2]
        
        model_withoutRep = glmer(Score ~ (1|Name) + (1|TRID) + (1|Name:TRID), family='binomial',data = temp_df) 
        p_rep[i] = anova(main_model,model_withoutRep)$`Pr(>Chisq)`[2]
        
        model_withoutG = glmer(Score ~  (1|replicate:TRID) + (1|TRID),family='binomial', data = temp_df) 
        p_G[i] = anova(main_model,model_withoutG)$`Pr(>Chisq)`[2]
        
        model_withoutE = glmer(Score ~ (1|Name),family='binomial', data = temp_df) 
        p_E[i] = anova(main_model,model_withoutE)$`Pr(>Chisq)`[2]
      }
      else {
        main_model <- lmer(Score ~ (1|Name) + (1|TRID) + (1|Name:TRID) + (1|replicate:TRID), data = temp_df)
        
        main_model %>%
          VarCorr() %>%
          as.data.frame() -> estimates
        
        VarRes[i] <- estimates$vcov[5]
        
        p_mean[i] <-  coef(summary(main_model))[5]
        
        test_model <- lmer(Score ~ (1|Name) + (1|TRID) + (1|Name:TRID) + (1|replicate:TRID), data = temp_df, REML = F)
        
        model_withoutGxE = lmer(Score ~ (1|Name) + (1|TRID) + (1|replicate:TRID), data = temp_df, REML = F) 
        p_GxE[i] = anova(test_model,model_withoutGxE)$`Pr(>Chisq)`[2]
        
        model_withoutRep = lmer(Score ~ (1|Name) + (1|TRID) + (1|Name:TRID), data = temp_df, REML = F) 
        p_rep[i] = anova(test_model,model_withoutRep)$`Pr(>Chisq)`[2]
        
        model_withoutG = lmer(Score ~  (1|replicate:TRID) + (1|TRID), data = temp_df, REML = F) 
        p_G[i] = anova(test_model,model_withoutG)$`Pr(>Chisq)`[2]
        
        model_withoutE = lmer(Score ~ (1|Name), data = temp_df, REML = F) 
        p_E[i] = anova(test_model,model_withoutE)$`Pr(>Chisq)`[2]
      }
      
      main_model %>%
        VarCorr() %>%
        as.data.frame() -> estimates
      
      VarGxE[i] <- estimates$vcov[1]
      VarG[i] <- estimates$vcov[2]
      VarRep[i] <- estimates$vcov[3]
      VarE[i] <- estimates$vcov[4]
      
      mean[i] <- coef(summary(main_model))[1]
      n_obs_GxE[i] <- ranef(main_model)$`Name:TRID` %>%
        nrow()
      n_obs_Rep[i] <- ranef(main_model)$`replicate:TRID` %>%
        nrow()
      H2[i] = round(VarG[i]/(VarG[i]+VarGxE[i]/n_obs_E[i]+VarRes[i]/n_obs_Rep[i]),2)
    }
    #Is there more than 1 environment tested?
    else {
      #Is the trait binary?
      if (all(unique(temp_df$Score) %in% c(0,1))) {
        main_model <- glmer(Score ~ (1|Name) + (1|replicate), family='binomial',data = temp_df)
        # the residual variance in a logistic regression is fixed 
        VarRes[i] <- (pi ^ 2) / 3
        
        p_mean[i] <-  coef(summary(main_model))[4]
        
        model_withoutRep = glmer(Score ~ (1|Name), family='binomial',data = temp_df) 
        p_rep[i] = anova(main_model,model_withoutRep)$`Pr(>Chisq)`[2]
        
        model_withoutG = glmer(Score ~ (1|replicate),family='binomial', data = temp_df) 
        p_G[i] = anova(main_model,model_withoutG)$`Pr(>Chisq)`[2]
      }
      else {
        main_model <- lmer(Score ~ (1|Name) + (1|replicate), data = temp_df)
        
        main_model %>%
          VarCorr() %>%
          as.data.frame() -> estimates
        
        VarRes[i] <- estimates$vcov[3]
        
        p_mean[i] <-  coef(summary(main_model))[5]
        
        test_model <- lmer(Score ~ (1|Name) + (1|replicate), data = temp_df, REML = F)
        
        model_withoutRep = lmer(Score ~ (1|Name), data = temp_df, REML = F) 
        p_rep[i] = anova(test_model,model_withoutRep)$`Pr(>Chisq)`[2]
        
        model_withoutG = lmer(Score ~ (1|replicate), data = temp_df, REML = F) 
        p_G[i] = anova(test_model,model_withoutG)$`Pr(>Chisq)`[2]
      }
      
      VarG[i] <- estimates$vcov[1]
      VarRep[i] <- estimates$vcov[2]
      
      mean[i] <- coef(summary(main_model))[1]
      n_obs_Rep[i] <- ranef(main_model)$`replicate` %>%
        nrow()
      
      H2[i] = round(VarG[i]/(VarG[i]+VarRes[i]/n_obs_Rep[i]),2)
    }
  }
  #Construct return tibble
  return_tibble <- tibble(Trait = trait,
                          Mean = mean,
                          P.value_mean = p_mean,
                          Genomic_variance = VarG,
                          P.value_genomic = p_G,
                          n_genomic = n_obs_G,
                          Environmental_variance = VarE,
                          P.value_environmental = p_E,
                          n_environmental = n_obs_E,
                          GenomicxEnvironmental_variance = VarGxE,
                          P.value_GxE = p_GxE,
                          n_GxE = n_obs_GxE,
                          Replicate_variance = VarRep,
                          P.value_replicates = p_rep,
                          n_replicates = n_obs_Rep,
                          Residual_variance = VarRes,
                          n_datapoints = n_datapoints,
                          Broad_sense_heritability = H2)
  return(return_tibble)
}

#===============================================================================
#Blues
#===============================================================================

# BLUEs multienvironmental trials (non-binary) function ----
# Function to refit and extract BLUEs from non-binary traits scored in more environments 
BLUEs_multienvironmental_nonbinary <- function(File,sigfile_multi,trait){
  
  TEMP_DF_obs = File[] %>%
    filter(DescriptionOfTrait == trait)
  
  Relevance_of_Env <- sigfile_multi$P.value_environmental[which(sigfile_multi$Trait==trait)]
  Relevance_of_GxE <- sigfile_multi$P.value_GxE[which(sigfile_multi$Trait==trait)]
  Relevance_of_Rep <- sigfile_multi$P.value_replicates[which(sigfile_multi$Trait==trait)]
  if (is.na(Relevance_of_Env)) {
    Relevance_of_Env <- 1
    Relevance_of_GxE <- 1
  }
  
  TEMP_DF_obs$replicate = as.factor(TEMP_DF_obs$replicate)
  TEMP_DF_obs$Name = as.factor(TEMP_DF_obs$Name)
  TEMP_DF_obs$TRID = as.factor(TEMP_DF_obs$TRID)
  
  # model for relevance of environments, GxE and reps
  if (Relevance_of_Env< 0.05 & Relevance_of_GxE< 0.05 &  Relevance_of_Rep<0.05){
    
    model <- lmer(Score ~ Name + (1|TRID) + (1|Name:TRID) + (1|replicate:TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of environments, GxE
  if (Relevance_of_Env< 0.05 & Relevance_of_GxE< 0.05 &  Relevance_of_Rep>=0.05){
    
    model <- lmer(Score ~ Name + (1|TRID) + (1|Name:TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of environments, reps
  if (Relevance_of_Env< 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep<0.05){
    
    model <- lmer(Score ~ Name + (1|TRID) + (1|replicate:TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of GxE and reps
  if (Relevance_of_Env>= 0.05 & Relevance_of_GxE< 0.05 &  Relevance_of_Rep<0.05){
    
    model <- lmer(Score ~ Name + (1|Name:TRID) + (1|replicate:TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of E
  if (Relevance_of_Env< 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep>=0.05){
    
    model <- lmer(Score ~ Name + (1|TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of GxE
  if (Relevance_of_Env>= 0.05 & Relevance_of_GxE<0.05 &  Relevance_of_Rep>=0.05){
    
    model <- lmer(Score ~ Name + (1|Name:TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of Reps
  if (Relevance_of_Env>= 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep<0.05){
    model <- lmer(Score ~ Name + (1|replicate), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of nothing
  if (Relevance_of_Env>= 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep>=0.05){
    model <- gls(Score ~ Name, data = TEMP_DF_obs) 
    BLUEs <- model$coefficients[2:length(model$coefficients)]
  }
  
  BLUE_df = as.data.frame(BLUEs)
  BLUE_df2 = cbind(rownames(BLUE_df),BLUE_df$BLUEs)
  col2name = trait
  colnames(BLUE_df2)= c("Line",col2name)
  return(BLUE_df2)
  
}

ProFABA_traits <- unique(ProFABA_mixed$DescriptionOfTrait)

for (i in 1:length(ProFABA_traits)) {
  if (i %in% 1:27) {
    next
  }
  if (i == 28) {
    ProFABA_blue <- BLUEs_multienvironmental_nonbinary(ProFABA_mixed, ProFABA_print, ProFABA_traits[i]) %>%
      as_tibble()
  }
  else {
    BLUEs_multienvironmental_nonbinary(ProFABA_mixed, ProFABA_print, ProFABA_traits[i]) %>%
      as_tibble() -> temp
    full_join(ProFABA_blue, temp) -> ProFABA_blue
  }
}

norfab_traits <- unique(norfab_mixed$DescriptionOfTrait)

for (i in 1:length(norfab_traits)) {
  if (i == 33) {
    break
  }
  if (i == 1) {
    norfab_blue <- BLUEs_multienvironmental_nonbinary(norfab_mixed, norfab_print, norfab_traits[i]) %>%
      as_tibble()
  }
  else {
    BLUEs_multienvironmental_nonbinary(norfab_mixed, norfab_print, norfab_traits[i]) %>%
      as_tibble() -> temp
    full_join(norfab_blue, temp) -> norfab_blue
  }
}
norfab_blue %>%
  mutate(Line = str_remove(Line, pattern = "Name")) -> norfab_blue

ProFABA_blue %>%
  mutate(Line = str_remove(Line, pattern = "Name")) -> ProFABA_blue

#===============================================================================
#Binary AVG
#===============================================================================
ProFABA %>%
  filter(DescriptionOfTrait %in% Weird_traits) %>%
  group_by(Name, DescriptionOfTrait, TRID) %>%
  summarise(unique = length(unique(Score))) -> tainted

ProFABA %>%
  filter(DescriptionOfTrait %in% Weird_traits) -> cat_pro
cat_pro %>%
  mutate(unique = get_corresponding_values(cat_pro, tainted, c("Name", "DescriptionOfTrait", "TRID"),
                                           c("Name", "DescriptionOfTrait", "TRID"), "unique")) %>%
  filter(unique == 1) %>%
  select(-c("unique", "TRID")) -> intermediate

intermediate %>%
  filter(DescriptionOfTrait == Weird_traits[1]) %>%
  mutate(Score = ifelse(Score < 1.6, "Narrow elongated", Score)) %>%
  mutate(Score = ifelse(Score >= 1.6 & Score < 2.5, "Sub elliptic", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Sub orbicular", Score)) %>%
  mutate(Score = ifelse(Score == 4, "Mixed", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) -> store_df

intermediate %>%
  filter(DescriptionOfTrait == Weird_traits[3] & Score %in% c(1, 2, 5, 6)) %>%
  mutate(Score = ifelse(Score == 1, "White", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Violet", Score)) %>%
  mutate(Score = ifelse(Score == 5, "Pink", Score)) %>%
  mutate(Score = ifelse(Score == 6, "Red", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

intermediate %>%
  filter(DescriptionOfTrait == Weird_traits[4] & Score %in% 1:3) %>%
  mutate(Score = ifelse(Score == 1, "White", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Coloured", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Spotted", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

intermediate %>%
  filter(DescriptionOfTrait == Weird_traits[6] & Score %in% 1:3) %>%
  mutate(Score = ifelse(Score == 1, "Flat", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Angular", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Round", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df

intermediate %>%
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

intermediate %>%
  filter(DescriptionOfTrait == Weird_traits[8] & Score %in% 1:4) %>%
  mutate(Score = ifelse(Score == 1, "Black", Score)) %>%
  mutate(Score = ifelse(Score == 2, "Grey", Score)) %>%
  mutate(Score = ifelse(Score == 3, "Colorless", Score)) %>%
  mutate(Score = ifelse(Score == 4, "Mixed", Score)) %>%
  unite(Score, c(Score, DescriptionOfTrait)) %>%
  bind_rows(store_df) -> store_df


store_df %>%
  dplyr::select(Name, Score) %>%
  na.omit() %>%
  group_by(Name, Score) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Name), names_from = Score, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) -> store_df

store_df %>%
  mutate(across(.cols = where(is.numeric),.fns = ~replace(., . > 1, 1)))# %>%
write_csv("ProFABA_AVG_cat.csv")



color_vector <- c("White_Gray", "Beige", "Green", "Violet", "Black")

norfab %>%
  filter(DescriptionOfTrait == "Seed Coat Colour") %>%
  group_by(Name, TRID) %>%
  summarise(unique = length(unique(Score))) -> tainted

norfab %>%
  filter(DescriptionOfTrait == "Seed Coat Colour") -> norfab_color
norfab_color %>%
  mutate(unique = get_corresponding_values(norfab_color, tainted, c("Name", "TRID"),
                                           c("Name", "TRID"), "unique")) %>%
  filter(unique == 1) %>%
  select(Name, Score) %>%
  mutate(Score = color_vector[Score]) %>%
  na.omit() %>%
  group_by(Name, Score) %>%
  summarise(n = n()) -> intermediate

intermediate %>%
  ungroup() %>%
  pivot_wider(id_cols = Name, names_from = Score, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(.cols = where(is.numeric),.fns = ~replace(., . > 1, 1))) %>%
  write_csv("norfab_AVG_cat.csv")

#===============================================================================
# GxE try
#===============================================================================
var_norfab <- read_csv("./Data/Variances_norfab.csv")
var_ProFABA <- read_csv("./Data/Variances_ProFABA.csv")

color = c("White_Gray", "Beige", "Green", "Violet", "Black")

var_norfab %>%
  filter(!(Trait %in% color)) %>%
  filter(n_environmental > 1 & P.value_GxE < 0.05 & P.value_environmental < 0.05) -> var_norfab

BLUEs_GxE <- function(File,sigfile_multi,trait){
  
  TEMP_DF_obs = File[] %>%
    filter(DescriptionOfTrait == trait)
  
  Relevance_of_Env <- sigfile_multi$P.value_environmental[which(sigfile_multi$Trait==trait)]
  Relevance_of_GxE <- sigfile_multi$P.value_GxE[which(sigfile_multi$Trait==trait)]
  Relevance_of_Rep <- sigfile_multi$P.value_replicates[which(sigfile_multi$Trait==trait)]
  
  TEMP_DF_obs$replicate = as.factor(TEMP_DF_obs$replicate)
  TEMP_DF_obs$Name = as.factor(TEMP_DF_obs$Name)
  TEMP_DF_obs$TRID = as.factor(TEMP_DF_obs$TRID)
  
  # model for relevance of environments, GxE and reps
  if (Relevance_of_Rep<0.05){
    
    model <- lmer(Score ~ Name + (1|TRID) + (1|Name:TRID) + (1|replicate:TRID), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  # model for relevance of environments, GxE
  if (Relevance_of_Rep>=0.05){
    
    model <- lmer(Score ~ (1|Name) + (1|TRID) + Name:TRID, data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  }
  
  BLUE_df = as.data.frame(BLUEs)
  BLUE_df2 = cbind(rownames(BLUE_df),BLUE_df$BLUEs)
  col2name = trait
  colnames(BLUE_df2)= c("Line",col2name)
  return(BLUE_df2)
  
}

#===============================================================================
# GxE try
#===============================================================================
pheno <- read_csv("./Data/ProFABA_pheno_DMU.csv")

myEGmaker <- function(pheno) {
  myFit <- lme4::lmer(formula = Score ~ (1|Name) + TRID + (1|replicate:TRID) + (1|Name:TRID), data = pheno)
  fixed <- lme4::fixef(myFit)
  
  EG <- as.vector(fixed) %>%
    scale() %>%
    as.vector()
  TRID <- names(fixed[2:length(fixed)]) %>%
    str_remove(pattern = "TRID")
  trials <- unique(pheno$TRID)
  trials <- trials[!(trials %in% TRID)]
  TRID <- c(trials, TRID)
  ret <- tibble(TRID = TRID, EG = EG)
  return(ret)
  
}

myGxEmixedmodelmaker <- function(pheno, Trait) {
  library(data.table)
  pheno %>%
    filter(DescriptionOfTrait == Trait) -> pheno
  
  pheno %>%
    group_by(Name) %>%
    summarise(n = length(unique(TRID))) %>%
    filter(n > 2) -> envs
  
  pheno %>%
    filter(Name %in% envs$Name)
  
  EG <- pheno %>% 
    myEGmaker()
  
  # Organizing fixed and random effects ----
  ph <- pheno %>% 
    mutate(Score = as.numeric(Score)) %>% 
    inner_join(., EG, by = 'TRID')
  
  #myFit <- lme4::lmer(formula = Score ~ Name + (1|EG) + (1|replicate:EG) + Name:EG, data = ph)
  myFit <- lme4::lmer(formula = Score ~ factor(EG) + (EG|Name), data = ph)
  lme4::ranef(myFit) -> random
 # summary(myFit) -> sum_model
#  coef <- sum_model$coefficients
 # coef <- coef[rownames(coef) %like% "EG",]
 # Z_value <- coef[,1]/coef[,2]
  #names <- str_remove(rownames(coef), pattern = "Name") %>%
  #  str_remove(pattern = ":EG")
  
  #tibble(Name = names, trait = Z_value) -> ret
  
  #colnames(ret)[2] <- Trait
  return(random)
}

ProFABA <- read_csv(file = "./Data/ProFABA.csv")

model_tibble <- myGxEmixedmodelmaker(pheno, "Seed yield")

names_translation <- read_csv(file = "./Data/Translation.csv", col_names = c("Name", "GPID"))

model_tibble %>%
  mutate(Name = plyr::mapvalues(model_tibble$Name, from = ProFABA$Name, to = ProFABA$GPID, warn_missing = F)) -> model_tibble

model_tibble$Name <- plyr::mapvalues(model_tibble$Name, from = names_translation$GPID, to = names_translation$Name, warn_missing = F)

model_tibble %>%
  as.data.frame() -> myY

options(rgl.debug = TRUE)

myGD <- read_csv(file = "./Data/ProFABA_GD.csv")
myGM <- read_csv(file = "./Data/ProFABA_GM.csv")

options(rgl.debug = TRUE)
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

myGM %>%
  mutate(Name = str_replace(Name, "-", ".")) -> myGM

myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  SNP.MAF = 0.05,
  model = c("FarmCPU")
)

myFit <- lme4::lmer(formula = Score ~ Name + (1|EG) + (1|replicate:EG) + Name:EG, data = ph)