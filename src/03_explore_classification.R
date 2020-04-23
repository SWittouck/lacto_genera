# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 9f2ceb8, ggtree version 1.16.0

library(tidyverse)
library(tidygenomes)
library(ggtree)

source("src/functions.R")

dout <- "results/classification"
if (! dir.exists(dout)) dir.create(dout, recursive = T)
dout_paper <- "results/genus_taxonomy_paper"
if (! dir.exists(dout_paper)) dir.create(dout_paper, recursive = T)

# GENOME DATASET 2 

# load data
load("results/parsed/lgc_withpairs.rda")

# determine exclusivity
lgc_protein <- 
  lgc_protein %>% 
  add_exclusivity(similarity = cni) %>%
  add_phylogroup_measures() %>%  
  modify_at(
    "genomes", mutate_at, c("furthest_within", "closest_between"), 
    ~ translate(., from = genome, to = species)
  )

# write table with exclusivity of genomes
lgc_protein$genomes %>%
  left_join(lgc_protein$nodes) %>%
  write_csv(paste0(dout, "/genomes_exclusivity.csv"))

# write table with exclusivity of phylogroups
lgc_protein$phylogroups %>%
  write_csv(paste0(dout, "/phylogroups_exclusivity.csv"))

# make figure with exclusivities
lgc_protein %>%
  pluck("phylogroups") %>%
  mutate(phylogroup = str_c(phylogroup, " (", pg_genomes, ")")) %>%
  mutate(phylogroup_fct = fct_reorder(
    phylogroup, min_similarity_within, .desc = T
  )) %>%
  ggplot(aes(y = phylogroup_fct, x = min_similarity_within)) +
  geom_point(aes(col = exclusive)) +
  xlim(c(0.5, 1)) +
  scale_color_brewer(palette = "Paired") + 
  xlab("minimum CNI within") +
  ylab("phylogroup and nuber of species") +
  theme_bw() 
ggsave(
  paste0(dout, "/genera_exclusivity.pdf"), 
  units = "cm", width = 30, height = 20
)

# make tree with exclusivities
ggtree_augmented(lgc_protein) +
  geom_tiplab(aes(label = species, col = consensus_phylogroup_member)) +
  theme(legend.position = "bottom")
ggsave(
  paste0(dout, "/tree_full_protein_exclusivity.pdf"), 
  units = "cm", width = 60, height = 100
)

# phylogenetic order of all 31 phylogroups to use in plots
genera_ordered <- 
  read_csv("data/lactobacillaceae_genera_2019.csv") %>%
  arrange(order_phylogeny) %>%
  pull(phylogroup)

# visualize classification based on CNI - all reference genomes
genome_pairs_extended <- genome_pairs_extended(lgc_protein)
classification_plot(genome_pairs_extended, cni) +
  geom_vline(xintercept = 0.70) +
  xlab("core nucleotide identity (cni)")
ggsave(
  paste0(dout, "/classification_cni.pdf"), 
  units = "cm", width = 17, height = 20, scale = 1
)

# visualize classification based on CNI - best reference per genus
cni_cutoff <- 0.70
genome_pairs_extended %>%
  filter_top_hit_per_phylogroup(cni) %>%
  classification_plot(cni, genera_ordered = genera_ordered) +
  geom_vline(xintercept = cni_cutoff) +
  xlab("top CNI to a genus")
ggsave(
  paste0(dout, "/classification_best_cni.pdf"), 
  units = "cm", width = 17, height = 20, scale = 1
)

# number of genomes that can be classified through CNI
eval_cni <- 
  genome_pairs_extended %>% 
  classification_assessment(cni, cni_cutoff) %>%
  rename(cni = n)

# visualize classification based on ANI - best reference per genus
ani_cutoff <- 0.87
genome_pairs_extended %>%
  filter_top_hit_per_phylogroup(ani) %>%
  classification_plot(ani, genera_ordered = genera_ordered) +
  geom_vline(xintercept = ani_cutoff) +
  xlab("top ANI to a genus")
ggsave(
  paste0(dout, "/classification_best_ani.pdf"), 
  units = "cm", width = 17, height = 20, scale = 1
)

# number of genomes that can be classified through ANI
eval_ani <- 
  genome_pairs_extended %>% 
  classification_assessment(ani, ani_cutoff) %>%
  rename(ani = n)

# GENOME DATASET 1

# prepare AAIs from genome dataset 1
genomes_ds1 <- read_csv("results/parsed/genomes_ds1.csv")
genome_pairs_ds1 <- read_csv("results/parsed/genome_pairs_ds1.csv")
genome_pairs_ds1_enriched <-
  genome_pairs_ds1 %>%
  enrich_genome_pairs(genomes_ds1)

# visualize classification based on AAI - best reference per genus
aai_cutoff <- 0.68
genome_pairs_ds1_enriched %>%
  filter_top_hit_per_phylogroup(aai) %>%
  classification_plot(aai, genera_ordered = genera_ordered) +
  geom_vline(xintercept = aai_cutoff) +
  xlab("top AAI to a genus")
ggsave(
  paste0(dout, "/classification_best_aai.pdf"), 
  units = "cm", width = 17, height = 20, scale = 1
)
file.copy(
  paste0(dout, "/classification_best_aai.pdf"),
  paste0(dout_paper, "/figure_S7_aai_classification.pdf")
)

# number of genomes that can be classified through AAI
eval_aai <- 
  genome_pairs_ds1_enriched %>% 
  classification_assessment(aai, aai_cutoff) %>%
  rename(aai = n)

# visualize classification based on cAAI - best reference per genus
caai_cutoff <- 0.74
genome_pairs_ds1_enriched %>%
  filter_top_hit_per_phylogroup(caai) %>%
  classification_plot(caai, genera_ordered = genera_ordered) +
  geom_vline(xintercept = caai_cutoff) +
  xlab("top cAAI to a genus")
ggsave(
  paste0(dout, "/classification_best_caai.pdf"), 
  units = "cm", width = 17, height = 20, scale = 1
)

# number of genomes that can be classified through cAAI
eval_caai <- 
  genome_pairs_ds1_enriched %>% 
  classification_assessment(caai, caai_cutoff) %>%
  rename(caai = n)

# compare classification strategies
eval <- 
  reduce(list(eval_ani, eval_cni, eval_aai, eval_caai), full_join) %>%
  filter(classifiable) %>%
  select(- classifiable) %>%
  gather(key = "method", value = "count", - classification_result) %>%
  spread(key = classification_result, value = count) %>%
  mutate(percentage_classified = map2_dbl(
    `correctly classified`, unclassified, ~ .x / (.x + .y)
  )) %>%
  mutate_at("percentage_classified", round, digits = 3) %>%
  left_join(tribble(
    ~ method, ~ cutoff, ~ dataset,
    "ani", ani_cutoff, "genome dataset 2",
    "cni", cni_cutoff, "genome dataset 2",
    "aai", aai_cutoff, "genome dataset 1",
    "caai", caai_cutoff, "genome dataset 1"
  )) %>%
  select(method, cutoff, dataset, percentage_classified)

eval %>%
  arrange(desc(percentage_classified))
