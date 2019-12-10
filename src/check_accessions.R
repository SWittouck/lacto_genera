library(tidyverse)
library(RCurl)

genomes <- read_csv("input/table_S2_genomes.csv")

make_urls <- function(accession) {
  ftp_prefix <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
  paste0(
    ftp_prefix, str_sub(accession, 5, 7), "/", str_sub(accession, 8, 10), "/", 
    str_sub(accession, 11, 13), "/"
  )
}

results <-
  genomes %>%
  mutate_at("assembly_accession", ~ str_replace(., "GCF", "GCA")) %>%
  mutate(ftp_url = make_urls(assembly_accession)) %>%
  mutate(result = map(ftp_url, safely(getURL))) %>%
  mutate(ftp_contents = map_chr(result, ~ .[[1]]))

results %>%
  select(assembly_accession, ftp_contents) %>%
  View()
