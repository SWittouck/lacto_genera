# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2, ggtree version 1.16.0

library(tidyverse)
library(tidygenomes)
library(ggtree)

dout <- "results/exclusivity"
if (! dir.exists(dout)) dir.create(dout, recursive = T)

# load data
load("parsed/legen_v3_3_repr_withpairs.rda")

# determine exclusivity
translate <- function(x, from, to) {
  lut <- structure(to, names = from)
  x <- lut[x] %>% unname()
}
lgc_repr_v3_3_withpairs <- 
  lgc_repr_v3_3_withpairs %>% 
  add_exclusivity(similarity = cni) %>%
  add_phylogroup_measures() %>%  
  modify_at(
    "genomes", mutate_at, c("furthest_within", "closest_between"), 
    ~ translate(., from = genome, to = species)
  )

# write table with exclusivity of genomes
lgc_repr_v3_3_withpairs$genomes %>%
  left_join(lgc_repr_v3_3_withpairs$nodes) %>%
  write_csv(paste0(dout, "/genomes_exclusivity.csv"))

# write table with exclusivity of phylogroups
lgc_repr_v3_3_withpairs$phylogroups %>%
  write_csv(paste0(dout, "/phylogroups_exclusivity.csv"))

# make figure with exclusivities
lgc_repr_v3_3_withpairs %>%
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
ggtree_augmented(lgc_repr_v3_3_withpairs) +
  geom_tiplab(aes(label = species, col = consensus_phylogroup_member)) +
  theme(legend.position = "bottom")
ggsave(
  paste0(dout, "/tree_full_protein_exclusivity.pdf"), 
  units = "cm", width = 60, height = 100
)
