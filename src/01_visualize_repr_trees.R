# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2, ggtree version 1.16.0

library(tidyverse)
library(tidygenomes)
library(ggtree)

dout <- "results/trees_and_phylogroups"
if (! dir.exists(dout)) dir.create(dout, recursive = T)
dout_paper <- "results/genus_taxonomy_paper"
if (! dir.exists(dout_paper)) dir.create(dout_paper, recursive = T)

load("results/parsed/legen_v3_3_repr_tree_protein.rda")
load("results/parsed/legen_v3_3_repr_tree_dna.rda")
load("results/parsed/legen_v3_3_repr_tree_gc.rda")

# write table with phylogroup membership of species
species <-
  full_join(
    lgc_repr_v3_3_protein_tree$genomes %>% 
      left_join(lgc_repr_v3_3_protein_tree$nodes),
    lgc_repr_v3_3_gc_tree$genomes %>% 
      left_join(lgc_repr_v3_3_gc_tree$nodes),
    by = "species",
    suffix = c("_protein", "_gc")
  ) %>% 
    select(species, phylogroup_protein, phylogroup_gc) %>%
    mutate(same_phylogroup = phylogroup_protein == phylogroup_gc) 
write_csv(species, path = paste0(dout, "/species_phylogroups.csv"))

# visualize full trees
ggtree_augmented(
  lgc_repr_v3_3_protein_tree, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species, col = is_phylogroup_type)) +
  geom_nodelab(aes(label = node_label)) +
  xlim(c(0, 2.4)) +
  scale_color_brewer(palette = "Paired")
ggsave(
  paste0(dout, "/tree_full_protein.pdf"), 
  units = "cm", width = 40, height = 100
)
ggtree_augmented(
  lgc_repr_v3_3_dna_tree, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species, col = is_phylogroup_type)) +
  geom_nodelab(aes(label = node_label)) +
  xlim(c(0, 2.4)) +
  scale_color_brewer(palette = "Paired")
ggsave(
  paste0(dout, "/tree_full_dna.pdf"), 
  units = "cm", width = 40, height = 100
)
ggtree_augmented(
  lgc_repr_v3_3_gc_tree, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species, col = is_phylogroup_type)) +
  geom_nodelab(aes(label = node_label)) +
  xlim(c(0, 0.04)) +
  scale_color_brewer(palette = "Paired")
ggsave(
  paste0(dout, "/tree_full_gc.pdf"), 
  units = "cm", width = 40, height = 100
)

# visualize genus taxonomy figure 6
phylogroups_leuconostocaceae <- c(
  "Leuconostoc", "Weissella", "Fructobacillus", "Oenococcus", "Convivina"
)
lgc_repr_v3_3_protein_tree$phylogroups <-
  lgc_repr_v3_3_protein_tree$phylogroups %>%
  mutate(family = case_when(
    phylogroup %in% phylogroups_leuconostocaceae ~ "Leuconostocaceae",
    TRUE ~ "Lactobacillaceae"
  ))
lgc_repr_v3_3_protein_tree$genomes <-
  lgc_repr_v3_3_protein_tree$genomes %>%
  mutate(fontface = if_else(is_phylogroup_type, "bold.italic", "italic")) 
lgc_repr_v3_3_protein_tree$nodes <-
  lgc_repr_v3_3_protein_tree$nodes %>%
  mutate(support_bs = str_extract(node_label, "[0-9]+$")) %>%
  mutate_at("support_bs", as.integer) 
lgc_repr_v3_3_protein_tree$tree <-
  lgc_repr_v3_3_protein_tree$tree %>% add_rootbranch()
lgc_repr_v3_3_protein_tree %>%
  ggtree_augmented(
    layout = "circular", col = "grey", size = 0.1
  ) +
  geom_tiplab2(
    aes(angle = angle, label = species, col = family, fontface = fontface),
    size = 1.5, align = T, linesize = 0.2, offset = 0.05
  ) +
  geom_nodelab(
    aes(label = if_else(support_bs < 90, support_bs, as.integer(NA))),
    size = 1
  ) +
  xlim(c(0, 3.4)) +
  scale_color_brewer(palette = "Paired", na.value = "#b2df8a")
ggsave(
  paste0(dout_paper, "/figure_S4_tree_full_protein.pdf"), 
  units = "cm", width = 16, height = 16
)

# write table with genomes
lgc_repr_v3_3_protein_tree$genomes %>%
  left_join(lgc_repr_v3_3_protein_tree$nodes) %>%
  left_join(lgc_repr_v3_3_protein_tree$phylogroups) %>%
  select(assembly_accession = genome, species, phylogroup) %>%
  write_csv(paste0(dout_paper, "/table_S2_genomes.csv"))

# create overview trees
lgc_repr_v3_3_protein_tree_overview <-
  lgc_repr_v3_3_protein_tree %>%
  filter_genomes(is_phylogroup_type)
lgc_repr_v3_3_dna_tree_overview <-
  lgc_repr_v3_3_dna_tree %>%
  filter_genomes(is_phylogroup_type)
lgc_repr_v3_3_gc_tree_overview <-
  lgc_repr_v3_3_gc_tree %>%
  filter_genomes(is_phylogroup_type)

# visualize overview trees
ggtree_augmented(
  lgc_repr_v3_3_protein_tree_overview, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species)) +
  xlim(c(0, 1.8))
ggsave(
  paste0(dout, "/tree_overview_protein.pdf"), 
  units = "cm", width = 20, height = 15
)
ggtree_augmented(
  lgc_repr_v3_3_dna_tree_overview, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species)) +
  xlim(c(0, 1.8))
ggsave(
  paste0(dout, "/tree_overview_dna.pdf"), 
  units = "cm", width = 20, height = 15
)
ggtree_augmented(
  lgc_repr_v3_3_gc_tree_overview, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species)) +
  xlim(c(0, 0.03))
ggsave(
  paste0(dout, "/tree_overview_gc.pdf"), 
  units = "cm", width = 20, height = 15
)

# visualize genus taxonomy figure 7b
ggtree_augmented(
  lgc_repr_v3_3_protein_tree_overview, layout = "rectangular", col = "grey"
) +
  geom_tiplab(aes(label = species, col = family), align = T, offset = 0.05) +
  xlim(c(0, 2.2)) +
  scale_color_brewer(palette = "Paired", na.value = "#b2df8a")
ggsave(
  paste0(dout_paper, "/figure_7b_tree_overview_protein.pdf"), 
  units = "cm", width = 20, height = 15
)

# vizualize protein overview tree with metabolism type
heterofermentative <- 
  c(
    "L. buchneri", "L. fructivorans", "L. kunkeei", "L. brevis", 
    "L. malefermentans", "L. fermentum", "L. vaccinostercus", "L. rossiae",
    "W. viridescens", "O. oeni", "Leuc. mesenteroides", "C. intestini", 
    "F. fructosus"
  )
lgc_repr_v3_3_protein_tree_overview %>%
  modify_at("nodes", mutate, support_bs = str_extract(node_label, "^[0-9]+") %>% as.integer()) %>%
  ggtree_augmented(col = "grey50") +
  geom_tiplab(aes(label = species_short, col = species_short %in% heterofermentative)) +
  xlim(c(0, 2.5)) +
  geom_point(
    aes(shape = support_bs %>% `>=`(90) %>% as.character()), 
    size = 3, color = "grey50"
  ) +
  scale_shape_manual(
    name = "support",
    values = c("TRUE" = 16, "FALSE" = 1),
    label = c("TRUE" = "bs >= 90", "FALSE" = "bs < 90"),
    na.translate = F
  ) +
  scale_color_manual(
    name = "metabolism",
    values = c("TRUE" = "#336B87", "FALSE" = "#763626"),
    label = c("TRUE" = "heterofermentative", "FALSE" = "homofermentative")
  ) +
  theme(legend.position = "bottom")
ggsave(
  paste0(dout, "/tree_overview_protein_metabolism.pdf"), 
  bg = "white", units = "cm", width = 20, height = 20
)
