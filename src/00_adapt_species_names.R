# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2

library(tidyverse)
library(tidygenomes)

source("src/functions.R")

if (! dir.exists("results")) dir.create("results")
if (! dir.exists("results/parsed")) dir.create("results/parsed")

load("data/legen_v3_2_repr.rda")
load("data/legen_v3_3_repr_tree_dna.rda")
load("data/legen_v3_3_repr_tree_gc.rda")
load("data/legen_v3_3_repr_tree_protein.rda")
load("data/legen_v3_3_repr_withpairs.rda")
load("data/legen_v3_3_repr.rda")

lgc_repr_v3_2 <-
  lgc_repr_v3_2 %>% 
  update_species_names()
save(lgc_repr_v3_2, file = "results/parsed/legen_v3_2_repr.rda")
lgc_repr_v3_3_dna_tree <-
  lgc_repr_v3_3_dna_tree %>% 
  update_species_names()
save(
  lgc_repr_v3_3_dna_tree, 
  file = "results/parsed/legen_v3_3_repr_tree_dna.rda"
)
lgc_repr_v3_3_gc_tree <-
  lgc_repr_v3_3_gc_tree %>% 
  update_species_names()
save(
  lgc_repr_v3_3_gc_tree, 
  file = "results/parsed/legen_v3_3_repr_tree_gc.rda"
)
lgc_repr_v3_3_protein_tree <-
  lgc_repr_v3_3_protein_tree %>% 
  update_species_names()
save(
  lgc_repr_v3_3_protein_tree, 
  file = "results/parsed/legen_v3_3_repr_tree_protein.rda"
)
lgc_repr_v3_3_withpairs <-
  lgc_repr_v3_3_withpairs %>% 
  update_species_names()
save(
  lgc_repr_v3_3_withpairs, 
  file = "results/parsed/legen_v3_3_repr_withpairs.rda"
)
lgc_repr_v3_3 <-
  lgc_repr_v3_3 %>% 
  update_species_names()
save(lgc_repr_v3_3, file = "results/parsed/legen_v3_3_repr.rda")
