# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2

library(tidyverse)
library(tidygenomes)

source("src/functions.R")

if (! dir.exists("results")) dir.create("results")
if (! dir.exists("results/parsed")) dir.create("results/parsed")

# adapt species names of genome dataset 2
clusters_named <- read_csv("data/legen_v3_4/clusters_all_named.csv")
clusters_named$species <-
  clusters_named$species %>%
  str_remove(" \\(\\?\\)") %>%
  str_replace("Lactobacillus timonensis", "'Lactobacillus timonensis'") %>%
  str_replace("Lactobacillus terrae", "Lactobacillus metriopterae") %>%
  {
    n <- length(.[str_detect(., "species")])
    .[str_detect(., "species")] <- str_c("Unassigned species", 1:n, sep = " ")
    .
  }
clusters_named %>%
  mutate(species_short = abbreviate_species(species)) %>%
  select(cluster, species, species_short) %>%
  write_csv("results/parsed/clusters_all_named_adapted.csv")

# read genome dataset 1 
genomes_ds1 <- 
  readxl::read_xlsx(
    "data/Table S1 Leuconostococaceae and Lactobacillaceae revised.xlsx"
  ) %>%
  transmute(
    genome = str_c(Genus, Species, sep = " "), 
    phylogroup = `proposed new genus name`,
    genus = Genus
  ) %>%
  mutate(phylogroup = if_else(is.na(phylogroup), genus, phylogroup)) %>%
  select(- genus) %>%
  mutate_at("genome", ~ str_extract(., "^[^ ]+ [^ ]+")) %>%
  distinct() %>%
  mutate_at("genome", abbreviate_species)

# read AAI values of genome dataset 1 
aais_ds1 <- 
  readxl::read_xlsx(
    "data/Table S4 AAI_table_full_May 2019.xlsx",
    range = readxl::cell_limits(c(2, 2), c(NA, NA))
  ) %>%
  rename(genome_1 = `...1`) %>%
  gather(key = "genome_2", value = "aai", - genome_1) %>%
  mutate_at("aai", as.numeric) %>%
  mutate_at("aai", ~ . / 100) %>%
  filter(! is.na(aai)) %>%
  mutate_at(c("genome_1", "genome_2"), ~ str_extract(., "^[^ ]+ [^ .]+")) %>%
  mutate_at(c("genome_1", "genome_2"), abbreviate_species) 

# read cAAI values of genome dataset 1 
caais_ds1 <- 
  readxl::read_xlsx(
    "data/Table S3 cAAI_table_full_Aug 2019.xlsx",
    range = readxl::cell_limits(c(2, 2), c(NA, NA))
  ) %>%
  rename(genome_1 = `...1`) %>%
  gather(key = "genome_2", value = "caai", - genome_1) %>%
  mutate_at("caai", as.numeric) %>%
  mutate_at("caai", ~ . / 100) %>%
  filter(! is.na(caai)) %>%
  mutate_at(c("genome_1", "genome_2"), ~ str_extract(., "^[^ ]+ [^ .]+")) %>%
  mutate_at(c("genome_1", "genome_2"), abbreviate_species) 

# compare presence in tables for all genomes/species of genome dataset 1
genomes_ds1_all <-
  tibble(
    genome = unique(c(aais_ds1$genome_1, aais_ds1$genome_2)), present_aai = T
  ) %>%
  full_join(tibble(
    genome = unique(c(caais_ds1$genome_1, caais_ds1$genome_2)), present_caai = T
  )) %>%
  full_join(tibble(
    genome = unique(genomes_ds1$genome), present_genomes = T
  ))
genomes_ds1_all%>%
  count(present_aai, present_caai, present_genomes)
genomes_ds1_all %>%
  filter(is.na(present_aai) | is.na(present_caai) | is.na(present_genomes)) %>%
  View()

# merge the AAI and cAAI tables
genome_pairs_ds1 <- 
  inner_join(aais_ds1, caais_ds1) %>%
  filter(genome_1 %in% genomes_ds1$genome, genome_2 %in% genomes_ds1$genome)

# write preprocessed data of genome dataset 1
write_csv(genomes_ds1, "results/parsed/genomes_ds1.csv")
write_csv(genome_pairs_ds1, "results/parsed/genome_pairs_ds1.csv")
