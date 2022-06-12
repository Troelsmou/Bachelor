myGenotypeCount <- function(vector) {
  ret <- c()
  for (i in 1:length(vector)) {
    vector[i] %>%
      strsplit(split = "|", fixed = T) -> intermediate
    intermediate[[1]] %>%
      as.numeric() %>%
      sum() -> gentyp
      ret <- c(ret, gentyp)
  }
  return(ret)
}

genotypes <- read_table(file = "./Data/PROFABA_genotypes.vcf")

genotypes %>%
  mutate(CHROM = ifelse(CHROM == "chr6", 7, CHROM)) %>%
  mutate(CHROM = ifelse(CHROM == "chr5", 6, CHROM)) %>%
  mutate(CHROM = ifelse(CHROM == "chr4", 5, CHROM)) %>%
  mutate(CHROM = ifelse(CHROM == "chr3", 4, CHROM)) %>%
  mutate(CHROM = ifelse(CHROM == "chr2", 3, CHROM)) %>%
  mutate(CHROM = ifelse(CHROM == "chr1L", 2, CHROM)) %>%
  mutate(CHROM = ifelse(CHROM == "chr1S", 1, CHROM)) %>%
  rename(Name = ID) %>%
  rename(Chromosome = CHROM) %>%
  rename(Position = POS) -> genotypes

genotypes %>%
  dplyr::select(Name, Chromosome, Position) %>%
  write_csv(file = "ProFABA_GM.csv")

genotypes %>%
  select(-c(Chromosome, Position, REF, ALT, QUAL, FILTER, INFO, FORMAT)) %>%
  mutate(across(.cols = -Name, .fns = myGenotypeCount)) %>%
  as.data.frame() -> b

rownames(b) <- b$Name

b %>%
  select(-Name) %>%
  t() -> c

write.table(c, file = "ProFABA_GD.txt")

read.table("ProFABA_GD.txt") -> hej

hej %>%
  tibble::rownames_to_column(var = "Taxa") -> ny

write_csv(ny, "ProFABA_GD.csv")

