# dependencies: R version 3.6.1, tidyverse version 1.2.1, tidygenomes version
# 0.1.2, ggtree version 1.16.0, castor version 1.4.1 (for RED)

library(tidyverse)
library(tidygenomes)
library(ggtree)

dout <- "results/signature_genes"
if (! dir.exists(dout)) dir.create(dout, recursive = T)
dout_paper <- "results/genus_taxonomy_paper"
if (! dir.exists(dout_paper)) dir.create(dout_paper, recursive = T)

# load data
load("results/parsed/legen_v3_3_repr_withpairs.rda")

# add and map pres/abs patterns
lgc_repr_v3_3_withpairs <- 
  lgc_repr_v3_3_withpairs %>%
  add_orthogroup_measures() %>%
  add_phylogroup_measures() %>%
  filter_orthogroups(og_genomes > 1, og_genomes < 239) %>%
  add_patterns() %>%
  map_patterns()

# type species of phylogroups 
phylogroups_types <-
  lgc_repr_v3_3_withpairs$genomes %>%
  filter(is_phylogroup_type) %>%
  left_join(lgc_repr_v3_3_withpairs$nodes) %>%
  select(phylogroup, type_species = species)

# number of genomes and signature genes per phylogroup
lgc_repr_v3_3_withpairs$nodes %>%
  filter(is_phylogroup_ancestor) %>%
  left_join(lgc_repr_v3_3_withpairs$phylogroups) %>%
  left_join(lgc_repr_v3_3_withpairs$patterns) %>%
  left_join(phylogroups_types) %>%
  select(
    phylogroup, type_species, n_species = pg_genomes, 
    n_signature_genes = frequency
  ) %>%
  mutate(n_signature_genes = if_else(
    is.na(n_signature_genes) & n_species > 1, 0L, n_signature_genes
  )) %>%
  write_csv(paste0(dout_paper, "/table_S5_genera_signature_genes.csv"))

# tree with numbers of signature genes (no labels)
lgc_repr_v3_3_withpairs %>%
  ggtree_augmented() +
  geom_point(aes(size = frequency)) +
  scale_size(range = c(1, 15))
ggsave(
  paste0(dout, "/tree_full_protein_signature_genes.pdf"), 
  units = "cm", width = 40, height = 40
)

# tree with numbers of signature genes
lgc_repr_v3_3_withpairs %>% 
  add_phylogroup_color(n = 12) %>%
  modify_at("tree", ~ castor::date_tree_red(.) %>% .$tree) %>%
  modify_at(
    "genomes", mutate,
    fontface = if_else(is_phylogroup_type, "bold", "plain"),
    alpha = if_else(is_phylogroup_type, 2, 1)
  ) %>%
  ggtree_augmented(color = "grey50", alpha = 0.5) +
  geom_point(aes(size = frequency), alpha = 0.5) +
  geom_tiplab(aes(
    label = species_short, 
    col = if_else(is_phylogroup_type, as.character(NA), phylogroup_color)
  )) +
  xlim(c(0, 1.1)) +
  scale_color_brewer(palette = "Paired", guide = F, na.value = "black") + 
  scale_alpha(range = c(0.5, 1), guide = FALSE) +
  scale_size(range = c(1, 30)) + 
  theme(rect = element_blank()) + 
  theme(legend.position = c(0.1, 0.9)) 
ggsave(
  paste0(dout, "/tree_full_protein_signature_genes2.pdf"), 
  units = "cm", width = 50, height = 100, dpi = 300, bg = "transparent"
)

# circular version of signature genes tree
lgc_repr_v3_3_withpairs %>%
  modify_at("tree", ~ date_tree_red(.) %>% .$tree) %>%
  modify_at("nodes", mutate, branch_type = case_when(
    phylogroup == "no phylogroup" | is_phylogroup_ancestor ~ "ancestral",
    str_detect(phylogroup, "group") ~ "Lactobacillus\nphylogroup",
    TRUE ~ "other genus"
  )) %>%
  ggtree_augmented(
    aes(color = branch_type), 
    layout = "circular"
  ) +
  geom_point(aes(size = frequency)) +
  xlim(c(0, 1)) +
  scale_color_brewer(palette = "Paired", name = "branch type") + 
  scale_size(range = c(1, 10), guide = F) + 
  theme(
    rect = element_blank(),
    legend.position = "right"
  )
ggsave(
  paste0(dout, "/tree_full_protein_signature_genes_circular.pdf"), 
  units = "cm", width = 20, height = 16, bg = "transparent"
)

# frequent patterns of gene presence/absence
pdf(
  paste0(dout_paper, "/figure_6_upset_plot.pdf"), 
  width = 60 * 0.39, height = 80 * 0.39, onefile = F
)
lgc_repr_v3_3_withpairs %>%
  upset_plot(
    genome_name = species, 
    genome_bold = is_phylogroup_type, 
    genome_col = is_phylogroup_ancestor & ! is.na(is_phylogroup_ancestor),
    n = 300,
    freq_min = 4,
    phylogroups = T
  )
dev.off()

# exploration of signature genes with some exceptions
load("parsed/legen_v3_3_repr_withpairs.rda")
phylogroups <- lgc_repr_v3_3_withpairs$phylogroups$phylogroup
dir.create(paste0(dout, "/phylogroup_cores"))
for (phylogroup in phylogroups) {
  
  pg <-
    lgc_repr_v3_3_withpairs %>%
    filter_genomes(phylogroup == !! phylogroup) %>%
    add_orthogroup_measures() %>%
    filter_orthogroups(og_genomes == max(og_genomes))
  
  ogs <- pg$orthogroups$orthogroup
  
  pg_genomes <- nrow(pg$genomes)
  
  filename <- paste0(dout, "/phylogroup_cores/", phylogroup, ".pdf")
  
  pdf(filename, width = 30 * 0.39, height = 80 * 0.39, onefile = F)
  lgc_repr_v3_3_withpairs %>%
    filter_orthogroups(orthogroup %in% !! ogs) %>%
    add_orthogroup_measures() %>%
    filter_orthogroups(
      (og_genomes <= !! pg_genomes + 3) | og_genomes == max(og_genomes)
    ) %>%
    add_patterns() %>%
    add_phylogroup_color(n = 2) %>%
    upset_plot(
      genome_name = species, 
      genome_col = phylogroup_color, 
      genome_bold = is_phylogroup_type, 
      color_scale = scale_color_manual(
        values = c("1" = "#1f78b4", "2" = "#33a02c"), guide = "none"
      ),
      n = 300,
      frequency_labels = T
    )
  dev.off()
  
}
