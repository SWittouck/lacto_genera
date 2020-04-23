# Lactobacillus genus taxonomy

This repository contains code to analyze the genus-level (and family-level) taxonomy of the Lactobacillus Genus Complex (families Lactobacillaceae and Leuconostocaceae).

The results of these analyses have been published in IJSEM:

[Zheng, J., Wittouck, S., Salvetti, E., Franz, C. M. A. P., Harris, H. M. B., Mattarelli, P., O’Toole, P. W., Pot, B., Vandamme, P., Walter, J., Watanabe, K., Wuyts, S., Felis, G. E., Gänzle, M. G., & Lebeer, S. (2020). A taxonomic note on the genus Lactobacillus: Description of 23 novel genera, emended description of the genus Lactobacillus Beijerinck 1901, and union of Lactobacillaceae and Leuconostocaceae. International Journal of Systematic and Evolutionary Microbiology. https://doi.org/https://doi.org/10.1099/ijsem.0.004107](https://doi.org/10.1099/ijsem.0.004107)

## Data

From the [legen pipeline](https://github.com/SWittouck/legen_pipeline), version 3.4:

* genome_pairs.csv.zip: data_v3/similarities/genome_pairs.csv.zip
* genomes_clusters.csv: data_v3/genome_clusters/genomes_clusters.csv
* clusters_all_named.csv: data_v3/taxonomy/clusters_all_named.csv
* outgroup_genomes.tsv: data_v3/outgroups/outgroup_genomes.tsv
* pangenome_legen_v3_3: data_v3/representatives_v3_3/pangenome/OrthoFinder/Orthogroups/
* lacto_dna.treefile: data_v3/representatives_v3_3/tree_dna/lacto_dna.treefile
* lacto_protein.treefile: data_v3/representatives_v3_3/tree_protein/lacto_protein.treefile
* lacto_gc.treefile: data_v3/representatives_v3_3/tree_gc/lacto_gc.treefile

Other data:

* lactobacillaceae_genera_2019.csv: constructed manually using the IJSEM manuscript; the column order_phylogeny is based on figure 5
* table_S1_corrected:
    * downloaded from <https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.004107#supplementary_data>
    * deleted row 158 (exact duplicate of row 149; Lactobacillus mulieris)
    * cell C53: changed to "Lactobacillus algidus" (added genus name)
    * cell O53: changed to "Dellaglioa algida" (added genus name)
    * cell C245: changed to "Lactobacillus aquaticus" (removed the single quote)
    * cell O324: changed to "Pediococcus acidilactici" (added genus name)
    * cells C300 - C304: fixed double spaces
* table_S3:
    * downloaded from <https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.004107#supplementary_data>
* table_S4:
    * downloaded from <https://www.microbiologyresearch.org/content/journal/ijsem/10.1099/ijsem.0.004107#supplementary_data>

## Scripts

All scripts can be found in the src folder: 

* 00_adapt_species names.R: make some slight changes to the species names (compared to the legen pipeline)
* 01_visualize_repr_trees.R: visualize various phylogenetic trees of representative genomes of LGC species
* 02_explore_signature_genes.R: explore signature genes (core and exclusive genes for a clade) in the LGC
* 03_explore_classification.R: explore how new genomes can be classified on the genus level in the new taxonomy
