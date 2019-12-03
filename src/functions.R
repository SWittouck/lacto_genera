update_species_names <- function(tg) {
  
  tg$genomes$species <-
    tg$genomes$species %>%
    str_remove(" \\(\\?\\)") %>%
    str_replace("Lactobacillus timonensis", "'Lactobacillus timonensis'") %>%
    str_replace("Lactobacillus terrae", "Lactobacillus metriopterae") %>%
    {
      n <- length(.[str_detect(., "species")])
      .[str_detect(., "species")] <- str_c("Unassig. species", 1:n, sep = " ")
      .
    }
  
  tg
  
}