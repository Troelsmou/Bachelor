Chr_vector <- c("Chr1S", "Chr1L", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6")

ProFABA_Gapit <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/GAPIT/SNPs_ProFABA.csv")
ProFABA_phen <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Data/ProFABA.csv")
norfab_gapit <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/GAPIT/SNPs_norfab.csv")
norfab_phen <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/Data/norfab.csv")
norfab_color <- read_csv("C:/Users/Troels/Desktop/Skole/Uni/Molekylærbiologi/Hestebønner/GAPIT/norfab_color.csv")

#Takes the first dataframe with name1 vector in it, finds the matches in df2
#in name2 and returns the corresponding values from values vector in df2 in order.
get_corresponding_values <- function(df1, df2, name1, name2, values) {
  ret <- c()
  if (length(name1) == 1) {
    searchkey <- df1[[name1]]
    search <- df2[[name2]]
    searchreturn <- df2[[values]]
    
    for (i in 1:length(searchkey)) {
      x <- unique(searchreturn[which(search == searchkey[i])])
      ret <- c(ret, x)
    }
  }
  if (length(name1) == 2) {
    searchkey1 <- df1[[name1[1]]]
    searchkey2 <- df1[[name1[2]]]
    search1 <- df2[[name2[1]]]
    search2 <- df2[[name2[2]]]
    searchreturn <- df2[[values]]
    for (i in 1:length(searchkey1)) {
      match1 <- searchkey1[i] == search1
      match2 <- searchkey2[i] == search2
      match <- match1 & match2
      x <- unique(searchreturn[match])
      ret <- c(ret,x)
    }
  }
  if (length(name1) == 3) {
      searchkey1 <- df1[[name1[1]]]
      searchkey2 <- df1[[name1[2]]]
      searchkey3 <- df1[[name1[3]]]
      search1 <- df2[[name2[1]]]
      search2 <- df2[[name2[2]]]
      search3 <- df2[[name2[3]]]
      searchreturn <- df2[[values]]
      for (i in 1:length(searchkey1)) {
        match1 <- searchkey1[i] == search1
        match2 <- searchkey2[i] == search2
        match3 <- searchkey3[i] == search3
        match <- match1 & match2 & match3
        x <- unique(searchreturn[match])
        ret <- c(ret,x)
      }
  }
  return(ret)
}
#===============================================================================
#ProFABA GAPIT
#===============================================================================

ProFABA_traits <- c("End of flowering", "Flowering time", "Flower number per node","Flowering period duration","Branch number per plant",
            "Aerial tillers number per plant (tillering higher than basal)", "Plant height at end of flowering", "Plant number",
            "Plant height at maturity", "First pod position", "Pod number per node", "Pod length", "Seed number per pod",
            "Seed yield","Hundred seed weight", "Number of seeds", "Seed health status", "Rust disease score (plant reaction-based scale)",
            "Rust disease score (percentage of covered leaf area)", "Chocolate spot disease score", "Aphid resistance index",
            "Broomrape infestation index", "Ascochyta blight severity score (leaves)", "Ascochyta blight severity score (pods)",
            "Herbicide damage", "Downy mildew (Peronospora vicia)", "Virus", "Sitona resistance index",
            "Percent of bruchid-infested seeds", "Field emergence", "Delayed emergence", "Lodging incidence",
            "Stem anthocyanin coloration")

ProFABA_Gapit %>%
  mutate(Trait = ifelse(Trait == "Aerial tillers", ProFABA_traits[6], Trait)) %>%
  mutate(Trait = ifelse(Trait == "Branch number", ProFABA_traits[5], Trait)) %>%
  mutate(Trait = ifelse(Trait == "Downy mildew", ProFABA_traits[26], Trait)) %>%
  mutate(Trait = ifelse(Trait == "Flower number", ProFABA_traits[3], Trait)) %>%
  mutate(Trait = ifelse(Trait == "Rust disease score", ProFABA_traits[19], Trait)) -> ProFABA_Gapit
ProFABA_non_bin <- ProFABA_Gapit %>%
  filter(Trait %in% ProFABA_traits)

ProFABA_non_bin %>%
  mutate(PDID = get_corresponding_values(ProFABA_non_bin, ProFABA_phen, "Trait", "DescriptionOfTrait", "PDID")) -> db

ProFABA_bin <- ProFABA_Gapit %>%
  filter(!(Trait %in% ProFABA_traits))

ProFABA_bin %>%
  mutate(PDID = ifelse(Trait %like% "seed testa ground", 21, NA)) %>%
  mutate(PDID = ifelse(Trait %like% "seed shape", 19, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "wing petal color", 6, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "standard petal color", 5, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "hilum color", 22, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "leaflet shape", 2, PDID)) %>%
  bind_rows(db) -> db
  
ProFABA_phen %>%
  group_by(PDID, TRID) %>%
  summarise(count = length(unique(Name))) -> count_df
  
db %>%
  mutate(Number_of_individuals = get_corresponding_values(db, count_df, c("PDID", "Trial"), c("PDID", "TRID"), "count")) %>%
  rename(SNP_ID = SNP, MAF = maf, GWAS_method = Method) %>%
  mutate(Project = "ProFABA", GWAS_options = "3 PCs", MAC = "", Phenotype_transformation = "",
         Population = "Core_ProFABA", Number_of_SNPs = n) -> db_end
#===============================================================================
#norfab GAPIT
#===============================================================================
norfab_traits <- unique(norfab_phen$DescriptionOfTrait)
  

norfab_gapit %>%
  mutate(GWAS_method = ifelse(Blink == 1, "Blink", "FarmCPU")) %>%
  mutate(Trait = ifelse(Trait == "Novitron Harm", "Novitron harm", Trait)) %>%
  mutate(Trait = ifelse(Trait == "Internode Length", "Internode_length", Trait)) %>%
  mutate(Trait = ifelse(Trait == "Sterile Tillers", "Sterile tillers", Trait)) %>%
  select(-c("Blink", "FarmCPU")) -> norfab_gapit

norfab_gapit %>%
  filter(Trait %in% norfab_traits) -> norfab_easy
norfab_easy %>%
  mutate(PDID = get_corresponding_values(norfab_easy, norfab_phen, c("Trial", "Trait"), c("TRID","DescriptionOfTrait"), "PDID")) -> db

norfab_hard <- norfab_gapit %>%
  filter(!(Trait %in% norfab_traits))

norfab_hard %>%
  mutate(PDID = ifelse(Trait %like% "Seed", 24, NA)) %>%
  mutate(PDID = ifelse(Trait %like% "(Per)", 29, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "(0-9)", 5, PDID)) %>%
  bind_rows(db) -> db

norfab_color %>%
  mutate(PDID = 25, Trait = str_c(Trait, "seed coat colour", sep = " ")) %>%
  rename(GWAS_method = Method) %>%
  bind_rows(db) -> db

norfab_phen %>%
  group_by(PDID, TRID) %>%
  summarise(count = length(unique(Name))) -> count_df

db %>%
  mutate(Number_of_individuals = get_corresponding_values(db, count_df, c("PDID", "Trial"), c("PDID", "TRID"), "count")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Project = "norfab", GWAS_options = "3 PCs", MAC = "", Phenotype_transformation = "",
         Population = "Core_norfab", Number_of_SNPs = n) %>%
  rename(SNP_ID = SNP, MAF = maf) %>%
  mutate(Trial = as.character(Trial)) %>%
  bind_rows(db_end) -> db_end

#===============================================================================
#EMMAX+EMMA200 norfab
#===============================================================================

norfab_emma_res %>%
  rename(Chromosome = chromosomes, Position = positions, P.value = scores, MAF = mafs) %>%
  select(-c("genotype_var_perc", macs)) -> norfab_emma_res_new

norfab_emma_res_new %>%
  mutate(Trait = ifelse(Trait == "Novitron Harm", "Novitron harm", Trait)) %>%
  mutate(Trait = ifelse(Trait == "Internode Length", "Internode_length", Trait)) %>%
  mutate(Trait = ifelse(Trait == "Sterile Tillers", "Sterile tillers", Trait)) -> norfab_emma_res_new

norfab_emma_res_new %>%
  filter(Trait %in% norfab_traits) -> norfab_easy
norfab_easy %>%
  mutate(PDID = get_corresponding_values(norfab_easy, norfab_phen, c("Trial", "Trait"), c("TRID","DescriptionOfTrait"), "PDID")) -> db

norfab_hard <- norfab_emma_res_new %>%
  filter(!(Trait %in% norfab_traits))

norfab_hard %>%
  filter(!((Trait == "White_Gray" & Trial == 26) | (Trait == "Black" & Trial == 30))) %>%
  mutate(PDID = ifelse(Trait %like% "(Per)", 29, NA)) %>%
  mutate(PDID = ifelse(Trait %like% "Rust" & Trial %in% c(25, 26), 28, PDID)) %>%
  mutate(PDID = ifelse(is.na(PDID), 25, PDID)) %>%
  mutate(Trait = ifelse(PDID == 25, str_c(Trait, "seed coat colour", sep = " "), Trait)) %>%
  bind_rows(db) -> db

norfab_phen %>%
  group_by(PDID, TRID) %>%
  summarise(count = length(unique(Name))) -> count_df

myGM <- read_tsv("./Data/myGM.txt")

db %>%
  mutate(Number_of_individuals = get_corresponding_values(db, count_df, c("PDID", "Trial"), c("PDID", "TRID"), "count")) %>%
  mutate(SNP_ID = get_corresponding_values(db, myGM, c("Chromosome", "Position"),c("Chromosome", "Position"), "Name")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Project = "norfab", GWAS_options = "", MAC = "", Phenotype_transformation = "",
         Population = "Core_norfab", Number_of_SNPs = n, GWAS_method = "EMMAX_EMMA200") %>%
  mutate(Trial = as.character(Trial)) %>%
  bind_rows(db_end) -> db_end

#===============================================================================
#EMMAX + EMMA200 ProFABA
#===============================================================================
ProFABA_emma <- read_csv("./EMMA EMMAX/ProFABA_EMMAX_Results.csv")

ProFABA_emma %>%
  rename(Chromosome = chromosomes, Position = positions, P.value = scores, MAF = mafs) %>%
  select(-c("genotype_var_perc", "macs")) -> ProFABA_emma_new

ProFABA_emma_new %>%
  filter(P.value < 2.733585e-06) -> ProFABA_emma_new

ProFABA_emma_new %>%
  filter(Trait %in% c(ProFABA_traits, "Germination")) -> ProFABA_easy
ProFABA_easy %>%
  mutate(PDID = get_corresponding_values(ProFABA_easy, ProFABA_phen, c("Trial", "Trait"), c("TRID","DescriptionOfTrait"), "PDID")) -> db

ProFABA_hard <- ProFABA_emma_new %>%
  filter(!(Trait %in% c(ProFABA_traits, "Germination")))

ProFABA_hard %>%
  mutate(PDID = ifelse(Trait %like% "Seed testa ground", 21, NA)) %>%
  mutate(PDID = ifelse(Trait %like% "Seed shape", 19, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Wing petal colour", 6, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Standard petal colour", 5, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Hilum color", 22, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Leaflet shape", 2, PDID)) %>%
  bind_rows(db) -> db

myGM <- read_csv("./Data/ProFABA_GM.csv")

ProFABA_phen %>%
  group_by(PDID, TRID) %>%
  summarise(count = length(unique(Name))) -> count_df

db %>%
  mutate(Number_of_individuals = get_corresponding_values(db, count_df, c("PDID", "Trial"), c("PDID", "TRID"), "count")) %>%
  mutate(SNP_ID = get_corresponding_values(db, myGM, c("Chromosome", "Position"),c("Chromosome", "Position"), "Name")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Project = "ProFABA", GWAS_options = "", MAC = "", Phenotype_transformation = "",
         Population = "Core_ProFABA", Number_of_SNPs = 21345, GWAS_method = "EMMAX_EMMA200") %>%
  mutate(PDID = as.character(PDID)) %>%
  bind_rows(db_end) -> db_end
#===============================================================================
#Blues
#===============================================================================

ProFABA_blue <- read_csv("./Mixed/ProFABA_Blues_Results.csv")

ProFABA_blue %>%
  dplyr::select(SNP, Position, P.value, maf, nobs, Trait) %>%
  filter(maf < 1 & P.value < 2.733585e-06 ) -> ProFABA_blue
ProFABA_blue %>%
  mutate(Chromosome = get_corresponding_values(ProFABA_blue, myGM, "SNP", "Name", "Chromosome")) %>%
  mutate(PDID = get_corresponding_values(ProFABA_blue, ProFABA_phen, "Trait", "DescriptionOfTrait", "PDID")) -> ProFABA_blue

ProFABA_phen %>%
  group_by(DescriptionOfTrait) %>%
  summarise(Trial = str_c(unique(TRID), collapse = ";")) -> trials

ProFABA_blue %>%
  mutate(Trial = get_corresponding_values(ProFABA_blue, trials, "Trait", "DescriptionOfTrait", "Trial")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Phenotype_transformation = "BLUE", GWAS_options = "", MAC = "",
         Project = "ProFABA", GWAS_method = "FarmCPU", Number_of_SNPs = 21335, Population = "Core_ProFABA") %>%
  rename(MAF = maf, Number_of_individuals = nobs, SNP_ID = SNP) %>%
  bind_rows(db_end) -> db_end

norfab_blue <- read_csv("./Mixed/norfab_Blues_Results.csv")

norfab_blue %>%
  dplyr::select(SNP, Position, P.value, maf, nobs, Trait) %>%
  filter(maf < 1 & P.value < 2.884172e-06 ) -> norfab_blue
norfab_blue %>%
  mutate(Chromosome = get_corresponding_values(norfab_blue, myGM, "SNP", "Name", "Chromosome")) %>%
  filter(Trait %in% norfab_traits & Trait != "Downy Mildew") -> norfab_easy
  
norfab_easy %>%
  mutate(PDID = get_corresponding_values(norfab_easy, norfab_phen, "Trait", "DescriptionOfTrait", "PDID")) -> norfab_easy

norfab_blue %>%
  filter(!(Trait %in% norfab_traits) | Trait == "Downy Mildew") %>%
  mutate(PDID = ifelse(Trait == "Downy Mildew", "36;46", NA)) %>%
  mutate(PDID = ifelse(Trait %like% "(Per)", "10;29", PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "(Cat)", 5, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Rust", "27;28", PDID)) -> norfab_hard
  
norfab_phen %>%
  group_by(PDID) %>%
  summarise(Trial = str_c(unique(TRID), collapse = ";")) -> trials

norfab_easy %>%
  mutate(Trial = get_corresponding_values(norfab_easy, trials, "PDID", "PDID", "Trial")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Phenotype_transformation = "BLUE", GWAS_options = "", MAC = "",
         Project = "norfab", GWAS_method = "FarmCPU", Number_of_SNPs = 21335, Population = "Core_norfab") %>%
  rename(MAF = maf, Number_of_individuals = nobs, SNP_ID = SNP) %>%
  bind_rows(db_end) -> db_end

norfab_hard %>%
  mutate(Trial = ifelse(PDID == 5, trials$Trial[trials$PDID%in%5], NA)) %>%
  mutate(Trial = ifelse(PDID == "36;46", "22;23;25;30", Trial)) %>%
  mutate(Trial = ifelse(PDID == "10;29", "22;25;26", Trial)) %>%
  mutate(Trial = ifelse(PDID == "27;28", "25;26;30", Trial)) %>%
  mutate(Chromosome = get_corresponding_values(norfab_hard, myGM, "SNP", "Name", "Chromosome")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Phenotype_transformation = "BLUE", GWAS_options = "", MAC = "",
         Project = "norfab", GWAS_method = "FarmCPU", Number_of_SNPs = 21335, Population = "Core_norfab") %>%
  rename(MAF = maf, Number_of_individuals = nobs, SNP_ID = SNP) %>%
  bind_rows(db_end) -> db_end

db_end$PDID <- as.character(db_end$PDID)  

#===============================================================================
#AVG
#===============================================================================
ProFABA_AVG <- read_csv("./Data/ProFABA_AVG_cat_Results.csv")
norfab_AVG <- read_csv("./Data/norfab_AVG_cat_Results.csv")

ProFABA_AVG %>%
  select(SNP, Position, P.value, maf, nobs, Trait) %>%
  filter(maf < 1) %>%
  filter(P.value < 2.733585e-06) %>%
  mutate(PDID = ifelse(Trait %like% "Seed testa ground", 21, NA)) %>%
  mutate(PDID = ifelse(Trait %like% "Seed shape", 19, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Wing petal colour", 6, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Standard petal colour", 5, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Hilum color", 22, PDID)) %>%
  mutate(PDID = ifelse(Trait %like% "Leaflet shape", 2, PDID)) -> ProFABA_AVG_new

ProFABA_phen %>%
  group_by(PDID) %>%
  summarise(Trial = str_c(unique(TRID), collapse = ";")) -> trials

ProFABA_AVG_new %>%
  mutate(Chromosome = get_corresponding_values(db, myGM, "SNP", "Name", "Chromosome")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Phenotype_transformation = "AVG", GWAS_options = "", MAC = "",
         Project = "ProFABA", GWAS_method = "FarmCPU", Number_of_SNPs = 21335, Population = "Core_ProFABA") %>%
  rename(MAF = maf, Number_of_individuals = nobs, SNP_ID = SNP) %>%
  mutate(Trial = get_corresponding_values(ProFABA_AVG_new, trials, "PDID", "PDID", "Trial")) %>%
  mutate(PDID = as.character(PDID)) %>%
  bind_rows(db_end) -> db_end

norfab_AVG %>%
  select(SNP, Position, P.value, maf, nobs, Trait) %>%
  filter(maf < 1) %>% 
  filter(P.value < 2.884172e-06) %>%
  mutate(Trait = str_c(Trait, "Seed Coat Colour", sep = " ")) %>%
  mutate(PDID = "25") -> norfab_AVG_new

norfab_phen %>%
  group_by(PDID) %>%
  summarise(Trial = str_c(unique(TRID), collapse = ";")) -> trials

myGM <- read_tsv("./Data/myGM.txt")

norfab_AVG_new %>%
  mutate(Chromosome = get_corresponding_values(norfab_AVG_new, myGM, "SNP", "Name", "Chromosome")) %>%
  mutate(Chromosome = Chr_vector[Chromosome], Phenotype_transformation = "AVG", GWAS_options = "", MAC = "",
         Project = "norfab", GWAS_method = "FarmCPU", Number_of_SNPs = 21335, Population = "Core_norfab") %>%
  rename(MAF = maf, Number_of_individuals = nobs, SNP_ID = SNP) %>%
  mutate(Trial = get_corresponding_values(norfab_AVG_new, trials, "PDID", "PDID", "Trial")) %>%
  bind_rows(db_end) -> db_end

db_end %>%
  mutate(Genotyping_platform = "Axiom 50k", SNP_mapping = "WE_map", row = row_number()) %>%
  write_csv("SNP_database.csv")
