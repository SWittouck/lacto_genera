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
#' @param genomes Three genomes
#' @param genome_identifier Variable of the genome table that corresponds to the
#'   given genomes
#' 
#' @return A tidygenomes object
#' 
#' @export
root_tree <- function(tg, genomes, genome_identifier = genome) {
  
  if (! "tree" %in% names(tg)) {
    stop("Tg should contain a tree")
  }
  
  genome_identifier <- rlang::enexpr(genome_identifier)
  
  tips <-
    tg$genomes %>%
    mutate(genome_identifier = !! genome_identifier) %>%
    filter(genome_identifier %in% !! genomes) %>%
    left_join(tg$nodes, by = "node") %>%
    pull(node) 
  
  if (! length(tips) == 3) {
    stop("Not all nodes were found")
  }
  
  tg$tree <- tg$tree %>% root_tree.phylo(tips) 
  
  tg
  
}

translate <- function(x, from, to) {
  lut <- structure(to, names = from)
  x <- lut[x] %>% unname()
}

complete_pairs <- function(pairs, object_1, object_2) {
  
  object_1 <- rlang::enexpr(object_1)
  object_2 <- rlang::enexpr(object_2)
  
  pairs_2 <-
    pairs %>%
    rename(object_1_new = !! object_2, object_2_new = !! object_1) %>%
    rename(!! object_1 := object_1_new, !! object_2 := object_2_new)
  
  bind_rows(pairs, pairs_2)
  
}

enrich_genome_pairs <- function(genome_pairs, genomes) {
  
  singleton_phylogroups <-
    genomes %>%
    count(phylogroup, name = "n_species") %>%
    filter(n_species == 1) %>%
    pull(phylogroup)
  
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
    mutate(same_phylogroup = phylogroup_query == phylogroup_reference) %>%
    mutate(phylogroup_query = if_else(
      phylogroup_query %in% singleton_phylogroups, "new genus", 
      phylogroup_query
    ))
  
}

genome_pairs_extended <- function(tg) {
  
  genomes <-
    tg$genomes %>%
    left_join(tg$nodes) %>%
    select(genome, phylogroup)
  
  enrich_genome_pairs(tg$pairs, genomes)
  
}

classification_plot <- function(genome_pairs, similarity) {
  
  similarity <- rlang::enexpr(similarity)
  
  genome_pairs %>%
    arrange(same_phylogroup) %>%
    ggplot(aes(x = !! similarity, y = genome_query, col = same_phylogroup)) +
    geom_point(size = 0.1)+
    facet_grid(rows = vars(phylogroup_query), scales = "free_y", space = "free_y") +
    scale_color_brewer(palette = "Paired", name = "correct reference genus") +
    ylab("genome") +
    theme_bw() +
    theme(
      text = element_text(size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom"
    )
  
}

prepare_tidygenomes <- function(genomes, tree, root_location) {
  
  as_tidygenomes(genomes) %>%
    add_tidygenomes(tree_protein) %>%
    root_tree(genomes = root_location, genome_identifier = species) %>%
    add_phylogroups(
      phylogroups %>% rename(genome_type = species_type), 
      genome_identifier = species_short
    )
  
}

classification_assessment <- function(genome_pairs_extended, similarity, min_value) {
  
  similarity <- rlang::enexpr(similarity)
  
  genome_pairs_extended %>%
    group_by(genome_query) %>%
    arrange(desc(!! similarity)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(genus_in_db = if_else(
      phylogroup_query == "new genus", "no", "yes")
    ) %>%
    mutate(classification_result = case_when(
      !! similarity < min_value ~ "unclassified",
      phylogroup_query == phylogroup_reference ~ "correctly classified",
      TRUE ~ "incorrectly classified"
    )) %>%
    count(genus_in_db, classification_result)
  
}

