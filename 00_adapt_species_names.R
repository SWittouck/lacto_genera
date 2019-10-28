# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2

library(tidyverse)
library(tidygenomes)

update_species_names <- function(tg) {
  
  tg$genomes$species <-
    tg$genomes$species %>%
    str_remove(" \\(\\?\\)") %>%
    str_replace("Lactobacillus timonensis", "'Lactobacillus timonensis'") %>%
    str_replace("Lactobacillus terrae", "Lactobacillus metriopterae") %>%
    {
      n <- length(.[str_detect(., "species")])
      .[str_detect(., "species")] <- str_c("Unassig. species", 1:n, sep = " ")
      .
    }
  
  tg
  
}

if (! dir.exists("parsed")) dir.create("parsed")

load("input/legen_v3_2_repr.rda")
load("input/legen_v3_3_repr_tree_dna.rda")
load("input/legen_v3_3_repr_tree_gc.rda")
load("input/legen_v3_3_repr_tree_protein.rda")
load("input/legen_v3_3_repr_withpairs.rda")
load("input/legen_v3_3_repr.rda")

lgc_repr_v3_2 <-
  lgc_repr_v3_2 %>% 
  update_species_names()
save(lgc_repr_v3_2, file = "parsed/legen_v3_2_repr.rda")
lgc_repr_v3_3_dna_tree <-
  lgc_repr_v3_3_dna_tree %>% 
  update_species_names()
save(lgc_repr_v3_3_dna_tree, file = "parsed/legen_v3_3_repr_tree_dna.rda")
lgc_repr_v3_3_gc_tree <-
  lgc_repr_v3_3_gc_tree %>% 
  update_species_names()
save(lgc_repr_v3_3_gc_tree, file = "parsed/legen_v3_3_repr_tree_gc.rda")
lgc_repr_v3_3_protein_tree <-
  lgc_repr_v3_3_protein_tree %>% 
  update_species_names()
save(
  lgc_repr_v3_3_protein_tree, file = "parsed/legen_v3_3_repr_tree_protein.rda"
)
lgc_repr_v3_3_withpairs <-
  lgc_repr_v3_3_withpairs %>% 
  update_species_names()
save(lgc_repr_v3_3_withpairs, file = "parsed/legen_v3_3_repr_withpairs.rda")
lgc_repr_v3_3 <-
  lgc_repr_v3_3 %>% 
  update_species_names()
save(lgc_repr_v3_3, file = "parsed/legen_v3_3_repr.rda")
