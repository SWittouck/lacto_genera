#' Abbreviate LGS genus names
#'
#' This function abbreviates the genus names in species names of the
#' Lactobacillus Genus Complex.
#' 
#' @param species A character vector of species names 
#' 
#' @return A character vector of abbreviated species names
abbreviate_species <- function(species) {
  species %>%
    str_replace_all("Lactobacillus", "L.") %>%
    str_replace_all("Pediococcus", "P.") %>%
    str_replace_all("Leuconostoc", "Lc.") %>%
    str_replace_all("Oenococcus", "O.") %>%
    str_replace_all("Fructobacillus", "Fb.") %>%
    str_replace_all("Weissella", "W.") %>%
    str_replace_all("Convivina", "C.") %>%
    str_replace_all("Unassigned", "Unassig.") %>%
    str_replace_all("pseudomesenteroides", "pseudomesent.")
}

#' Root tree given three tips
#'
#' This function roots a phylogenetic tree given three tip labels.
#'
#' Tips a, b and c define exactly one internal node in the unrooted tree. The
#' tree will be rooted on the branch leading from this node to tip a.
#' 
#' @param tree An object of class phylo
#' @param tips Three tip labels 
#' 
#' @return An object of class phylo
root_tree.phylo <- function(tree, tips) {
  
  if (! "phylo" %in% class(tree)) {
    stop("tree should be of class phylo")
  }
  
  # root on node defined by tips
  root_new <- 
    ape::mrca(tree) %>%
    {.[tips, tips]} %>%
    {.[. > length(tree$tip.label)]} %>%
    {.[. == max(.)]} %>%
    {.[1]}
  tree <- tree %>% ape::root.phylo(node = root_new)
  
  # resolve root node such that first tip is (part of) outgroup
  outgroup <-
    ape::mrca(tree) %>%
    {.[tips[1], ]} %>%
    {.[. != length(tree$tip.label) + 1]} %>%
    names()
  tree <- tree %>% ape::root.phylo(outgroup = outgroup, resolve.root = T)
  
  # divide root branch length
  if ("edge.length" %in% names(tree)) {
    n_tips <- length(tree$tip.label) 
    l <- sum(tree$edge.length[tree$edge[, 1] == n_tips + 1])
    tree$edge.length[tree$edge[, 1] == n_tips + 1] <- l / 2
  }
  
  # return tree
  tree
  
}

#' Root tree given three genomes
#'
#' This applies [root_tree.phylo] to the tree component of a tidygenomes object.
#' 
#' @param tg An tidygenomes object
#' @param root Three tips that identify the root (see [root_tree.phylo])
#' @param genome_identifier Variable of the genome table that corresponds to the
#'   given genomes
#' 
#' @return A tidygenomes object
root_tree <- function(tg, root, genome_identifier = genome) {
  
  if (! "tree" %in% names(tg)) {
    stop("Tg should contain a tree")
  }
  
  genome_identifier <- rlang::enexpr(genome_identifier)
  
  tips <-
    tg$genomes %>%
    mutate(genome_identifier = !! genome_identifier) %>%
    filter(genome_identifier %in% !! root) %>%
    left_join(tg$nodes, by = "node") %>%
    pull(node) 
  
  if (! length(tips) == 3) {
    stop("Not all nodes were found")
  }
  
  tg$tree <- tg$tree %>% root_tree.phylo(tips) 
  
  tg
  
}

#' Prepare tidygenomes object
#'
#' This function prepares a tidygenomes object from a genome table, a
#' phylogenetic tree and a root location in the tree.
#' 
#' @param genomes A genome table with a column `genome`
#' @param tree An object of class phylo with tips corresponding to genomes
#' @param root Three tips that identify the root (see [root_tree.phylo])
#' 
#' @return A tidygenomes object
prepare_tidygenomes <- function(genomes, tree, phylogroups, root) {
  
  phylogroups <- rename(phylogroups, genome_type = species_type)
  
  as_tidygenomes(genomes) %>%
    add_tidygenomes(tree) %>%
    root_tree(root, genome_identifier = species) %>%
    add_phylogroups(phylogroups, genome_identifier = species_short)
  
}

#' Translate elements of a character vector
#'
#' This function relaces the elements of a character vector using a translation
#' table given by a "from" and a "to" vector.
#' 
#' @param x Character vector to translate
#' @param from Character vector with unique elements to be replaced
#' @param to Character vector with unique elements to replace the elements of
#'   "from"
#' 
#' @return A translated character vector
translate <- function(x, from, to) {
  lut <- structure(to, names = from)
  x <- lut[x] %>% unname()
}

#' Complement the pairs of a pair table
#'
#' For each unique pair of objects (a, b), the function adds a row for the pair
#' (b, a) with the same information.
#' 
#' @param pairs Data frame where each row represents a pair of objects
#' @param object_1 Variable defining the first object (unquoted)
#' @param object_2 Variable defining the second object (unquoted)
#' 
#' @return A complemented table
complete_pairs <- function(pairs, object_1, object_2) {
  
  object_1 <- rlang::enexpr(object_1)
  object_2 <- rlang::enexpr(object_2)
  
  pairs_2 <-
    pairs %>%
    rename(object_1_new = !! object_2, object_2_new = !! object_1) %>%
    rename(!! object_1 := object_1_new, !! object_2 := object_2_new)
  
  bind_rows(pairs, pairs_2)
  
}

#' Enrich genome pair table with genome metadata
#'
#' This function completes a genome pair table (see [complete_pairs]) and adds
#' genome metadata supplied by a genome table.
#' 
#' @param genome_pairs Data frame containing the columns `genome_1` and
#'   `genome_2`
#' @param genomes Data frame containing the column `genome`
#' 
#' @return An enriched genome pair table
enrich_genome_pairs <- function(genome_pairs, genomes) {
  
  genome_pairs %>%
    complete_pairs(genome_1, genome_2) %>%
    rename(genome_query = genome_1, genome_reference = genome_2) %>%
    left_join(
      genomes %>% 
        rename(genome_query = genome, phylogroup_query = phylogroup),
      by = "genome_query"
    ) %>%
    left_join(
      genomes %>% 
        rename(genome_reference = genome, phylogroup_reference = phylogroup),
      by = "genome_reference"
    ) %>% 
    mutate(same_phylogroup = phylogroup_query == phylogroup_reference) 
  
}

#' Return an extended genome pair table
#'
#' This function returns an extended genome pair table from a tidygenomes
#' object, using the function [enrich_genome_pairs].
#' 
#' @param tg Tidygenomes object
#' 
#' @return An enriched genome pair table
genome_pairs_extended <- function(tg) {
  
  genomes <-
    tg$genomes %>%
    left_join(tg$nodes) %>%
    select(genome, phylogroup)
  
  enrich_genome_pairs(tg$pairs, genomes)
  
}

#' Plot the classification of genomes
#'
#' This function plots similarity values of query genomes to reference genomes,
#' grouped per phylogroup/genus of the query genomes. Similarity values to
#' reference genomes of the correct phylogroup are plotted in dark blue.
#'
#' @param genome_pairs A data frame containing the columns `genome_query`,
#'   `phylogroup_query` and `same_phylogroup`
#' @param similarity A genome-genome similarity expression
#' @param genera_ordered An optional vector defining the order in which all
#'   unique phylogroups/genera should appear in the plot
#'
#' @return A ggplot object
classification_plot <- 
  function(genome_pairs, similarity, genera_ordered = NULL) {
  
  similarity <- rlang::enexpr(similarity)
  
  genome_pairs %>%
    arrange(same_phylogroup) %>%
    {if (! is.null(genera_ordered)) {
      mutate_at(., "phylogroup_query", factor, levels = genera_ordered)
    } else {
      .
    }} %>%
    ggplot(aes(x = !! similarity, y = genome_query, col = same_phylogroup)) +
    geom_point(size = 0.1) +
    facet_grid(
      rows = vars(phylogroup_query), scales = "free_y"
    ) +
    scale_color_brewer(palette = "Paired", name = "own genus") +
    ylab("genome") +
    theme_bw() +
    theme(
      text = element_text(size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(angle = 0, face = "italic"),
      legend.position = "bottom"
    )
  
}

#' Filter top hits to reference phylogroups
#'
#' This function filters the similarity values of query genomes to reference
#' genomes, to retain only the best hit of each query genome to each reference
#' phylogroup.
#'
#' @param genome_pairs A data frame containing the columns `genome_query`,
#'   `phylogroup_query` and `phylogroup_reference`
#' @param similarity A genome-genome similarity expression
#'
#' @return A filtered genome pair table 
filter_top_hit_per_phylogroup <- function(genome_pairs, similarity) {
  
  similarity <- rlang::enexpr(similarity)
  
  genome_pairs %>%
    group_by(genome_query, phylogroup_query, phylogroup_reference) %>%
    arrange(desc(!! similarity)) %>%
    slice(1) %>%
    ungroup() 
  
}

#' Assess classification performance
#'
#' This function assesses the classification performance of a genome-genome
#' similarity measure and a cutoff. It counts how many genomes had a reference
#' genome of the correct genus/phylogroup to enable classification, and how many
#' genomes were correctly classified, incorrectly classified or unclassified.
#'
#' When the best hit of a genome scores above the cutoff to a reference genus,
#' it is classified to that genus. If it scores below the cutoff to all
#' reference genera, it is considered "unclassified".
#'
#' @param genome_pairs A data frame containing the columns `genome_query`,
#'   `phylogroup_query` and `phylogroup_reference`
#' @param similarity A genome-genome similarity expression
#' @param min_value A minimum cutoff similarity value for classification
#'
#' @return A genome count table
classification_assessment <- function(genome_pairs, similarity, min_value) {
  
  similarity <- rlang::enexpr(similarity)
  
  genome_pairs %>%
    group_by(genome_query) %>%
    mutate(classifiable = phylogroup_query[1] %in% phylogroup_reference) %>%
    arrange(desc(!! similarity)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(classification_result = case_when(
      !! similarity < min_value ~ "unclassified",
      phylogroup_query == phylogroup_reference ~ "correctly classified",
      TRUE ~ "incorrectly classified"
    )) %>%
    count(classifiable, classification_result)
  
}
