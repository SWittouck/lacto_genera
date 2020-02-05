# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2, ggtree version 1.16.0

library(tidyverse)
library(tidygenomes)
library(ggtree)

source("src/functions.R")

dout <- "results/trees_and_phylogroups"
if (! dir.exists(dout)) dir.create(dout, recursive = T)
dout_paper <- "results/genus_taxonomy_paper"
if (! dir.exists(dout_paper)) dir.create(dout_paper, recursive = T)

# import data 
genomes_clusters <- read_csv("data/legen_v3_4/genomes_clusters.csv")
clusters_species <- read_csv("results/parsed/clusters_all_named_adapted.csv")
outgroups <- 
  "data/legen_v3_4/outgroup_genomes.tsv" %>%
  read_tsv(col_names = c("genome", "species")) %>%
  mutate(species_short = species)
tree_dna <- ape::read.tree("data/legen_v3_4/lacto_dna.treefile")
tree_protein <- ape::read.tree("data/legen_v3_4/lacto_protein.treefile")
tree_gc <- ape::read.tree("data/legen_v3_4/lacto_gc.treefile")
# remark: we remove the phylogroup Acetilactobacillus because it's not present
# in genome dataset 2
phylogroups <- 
  read_csv("data/lactobacillaceae_genera_2019.csv") %>%
  filter(phylogroup != "Acetilactobacillus")

# preprocess genomes and add outgroups
genomes <- 
  genomes_clusters %>%
  left_join(clusters_species) %>%
  bind_rows(outgroups)

# define the root location of the tree
root <- c(
  "L. casei", "Listeria monocytogenes", 
  "Brochothrix thermosphacta"
)

# construct tidygenomes object per type of tree (protein, dna, gene content)
phylogroups <- rename(phylogroups, genome_type = species_type)
lgc_protein <- 
  prepare_tidygenomes(
    genomes = genomes, tree = tree_protein, phylogroups = phylogroups, 
    root = root, genome_identifier = species_short
  )
lgc_dna <- 
  prepare_tidygenomes(
    genomes = genomes, tree = tree_dna, phylogroups = phylogroups, 
    root = root, genome_identifier = species_short
  )
lgc_gc <- 
  prepare_tidygenomes(
    genomes = genomes, tree = tree_gc, phylogroups = phylogroups, 
    root = root, genome_identifier = species_short
  )
save(lgc_protein, file = "results/parsed/lgc_protein.rda")

# write table with phylogroup membership of species
species <-
  full_join(
    lgc_protein$genomes %>% left_join(lgc_protein$nodes),
    lgc_gc$genomes %>% left_join(lgc_gc$nodes),
    by = "species",
    suffix = c("_protein", "_gc")
  ) %>% 
  select(species, phylogroup_protein, phylogroup_gc) %>%
  mutate(same_phylogroup = phylogroup_protein == phylogroup_gc) 
write_csv(species, path = paste0(dout, "/species_phylogroups.csv"))

# visualize full trees
ggtree_augmented(
  lgc_protein, layout = "rectangular", col = "grey"
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
  lgc_dna, layout = "rectangular", col = "grey"
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
  lgc_gc, layout = "rectangular", col = "grey"
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
lgc_protein$phylogroups <-
  lgc_protein$phylogroups %>%
  mutate(family = case_when(
    phylogroup %in% phylogroups_leuconostocaceae ~ "Leuconostocaceae",
    TRUE ~ "Lactobacillaceae"
  ))
lgc_protein$genomes <-
  lgc_protein$genomes %>%
  mutate(fontface = if_else(is_phylogroup_type, "bold.italic", "italic")) 
lgc_protein$nodes <-
  lgc_protein$nodes %>%
  mutate(support_bs = str_extract(node_label, "[0-9]+$")) %>%
  mutate_at("support_bs", as.integer) 
lgc_protein$tree <-
  lgc_protein$tree %>% add_rootbranch()
lgc_protein %>%
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
lgc_protein$genomes %>%
  left_join(lgc_protein$nodes) %>%
  left_join(lgc_protein$phylogroups) %>%
  select(assembly_accession = genome, species, phylogroup) %>%
  write_csv(paste0(dout_paper, "/table_S2_genomes.csv"))

# create overview trees
lgc_protein_overview <- lgc_protein %>% filter_genomes(is_phylogroup_type)
lgc_dna_overview <- lgc_dna %>% filter_genomes(is_phylogroup_type)
lgc_gc_overview <- lgc_gc %>% filter_genomes(is_phylogroup_type)

# visualize overview trees
ggtree_augmented(
  lgc_protein_overview, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species)) +
  xlim(c(0, 1.8))
ggsave(
  paste0(dout, "/tree_overview_protein.pdf"), 
  units = "cm", width = 20, height = 15
)
ggtree_augmented(
  lgc_dna_overview, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species)) +
  xlim(c(0, 1.8))
ggsave(
  paste0(dout, "/tree_overview_dna.pdf"), 
  units = "cm", width = 20, height = 15
)
ggtree_augmented(
  lgc_gc_overview, layout = "rectangular", col = "grey"
  ) +
  geom_tiplab(aes(label = species)) +
  xlim(c(0, 0.03))
ggsave(
  paste0(dout, "/tree_overview_gc.pdf"), 
  units = "cm", width = 20, height = 15
)

# visualize genus taxonomy figure 7b
ggtree_augmented(
  lgc_protein_overview, layout = "rectangular", col = "grey"
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
    "W. viridescens", "O. oeni", "Lc. mesenteroides", "C. intestini", 
    "Fb. fructosus"
  )
lgc_protein_overview %>%
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
