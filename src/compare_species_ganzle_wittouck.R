library(tidyverse)
library(readxl)

species_ganzle_1 <-
  read_excel("input/Michael/Strain_lists_for_both_datasets_UCC_VR_190222.xlsx") %>%
  gather(key = "dataset_ganzle", value = "strain") %>%
  filter(strain != "-") %>%
  mutate(present = T) %>%
  spread(key = dataset_ganzle, value = present, fill = FALSE) %>%
  extract(strain, into = c("genus", "species", "subspecies", "strain_1"), 
          regex = "([A-Z][a-z]+)_([a-z]+)_?([a-z]+)?_(.+)") %>%
  filter(species == subspecies | is.na(subspecies)) %>%
  mutate(species = str_c(genus, species, sep = " ")) %>%
  rename(in_genome_dataset = `AAI/POCP/16S dataset (n = 218)`) %>%
  rename(in_16S_dataset = `Expanded 16S dataset (n = 282)`)

species_ganzle_2 <-
  read_excel("input/Michael/Lactobacillus proposed classification.xlsx") %>%
  select(genus = Genus, species = Species, strain_2 = Strain, genus_ganzle = ...11) %>%
  filter(! is.na(species)) %>%
  separate(species, into = c("species", "subspecies"), sep = " ", fill = "right") %>%
  filter(species == subspecies | is.na(subspecies)) %>%
  mutate(species = str_c(genus, species, sep = " "))

species_ganzle <- 
  left_join(species_ganzle_1, species_ganzle_2) %>%
  mutate_at("species", str_replace, "Leuconostoc", "Leuc\\.") %>%
  mutate_at("species", str_replace, "(?<=[A-Z])[a-z]{5,}", "\\.")

species_lebeer <-
  read_csv("input/Lebeer/species_phylogroups.csv") %>%
  mutate(
    status = case_when(
      str_detect(species, " \\(\\?\\)") ~ "best guess",
      str_detect(species, "[0-9]+") ~ "new/unidentified",
      TRUE ~ "normal"
    )
  ) %>%
  mutate_at("species", str_remove, " \\(\\?\\)") %>%
  mutate(in_lebeer = T) %>%
  rename(phylogroup_zheng = phylogroup)

species <- full_join(species_ganzle, species_lebeer)

species %>%
  select(
    species, in_genome_dataset, in_16S_dataset, 
    genus_ganzle, phylogroup_zheng, status
  ) %>%
  View()

species_wittouck <-
  read_csv("input/table_S3_clusters.csv") %>%
  select(species)

species <- full_join(species_ganzle, species_wittouck)

count(species, in_ganzle, in_wittouck)
