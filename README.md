# Lactobacillus genus taxonomy

This repository contains code to analyze the genus-level (and family-level) taxonomy of the Lactobacillus Genus Complex (families Lactobacillaceae and Leuconostocaceae). 

## Data

From the [legen](https://github.com/SWittouck/legen_pipeline) pipeline, version 3.4:

* genome_pairs.csv.zip: data_v3/similarities/genome_pairs.csv.zip
* genomes_clusters.csv: data_v3/genome_clusters/genomes_clusters.csv
* clusters_all_named.csv: data_v3/taxonomy/clusters_all_named.csv
* outgroup_genomes.tsv: data_v3/outgroups/outgroup_genomes.tsv
* pangenome_legen_v3_3: data_v3/representatives_v3_3/pangenome/OrthoFinder/Orthogroups/
* lacto_dna.treefile: data_v3/representatives_v3_3/tree_dna/lacto_dna.treefile
* lacto_protein.treefile: data_v3/representatives_v3_3/tree_protein/lacto_protein.treefile
* lacto_gc.treefile: data_v3/representatives_v3_3/tree_gc/lacto_gc.treefile

Other data:

* lactobacillaceae_genera_2019.csv: constructed manually using the manuscript; the column order_phylogeny is based on figure 5
* Tables S1, S3 and S4 are supplementary tables of the manuscript (with some small types corrected, which will normally also happen in the "official" versions)

## Scripts (src folder)

* 00_adapt_species names.R: make some slight changes to the species names (compared to the legen pipeline)
* 01_visualize_repr_trees.R: visualize various phylogenetic trees of representative genomes of LGC species
* 02_explore_signature_genes.R: explore signature genes (core and exclusive genes for a clade) in the LGC
* 03_explore_classification.R: explore how new genomes can be classified on the genus level in the new taxonomy
