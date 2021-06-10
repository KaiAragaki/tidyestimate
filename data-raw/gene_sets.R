# Prepares "gene_sets" dataset ----
library(dplyr)
library(utils)

utils::download.file("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/estimate/data/SI_geneset.RData?root=estimate",
                     destfile = "./data-raw/SI_geneset_old.RData", 
                     method = "curl")
load("./data-raw/SI_geneset_old.RData")

gene_sets <- SI_geneset |>
  t()  |> 
  as.data.frame() |>
  dplyr::slice(-1) |>
  dplyr::rename(stromal_signature = StromalSignature,
                immune_signature = ImmuneSignature) |> 
  dplyr::mutate(stromal_signature = case_when(stromal_signature == "WISP1" ~ "CCN4",
                                              stromal_signature == "GPR124" ~ "ADGRA2",
                                              stromal_signature == "LPPR4" ~ "PLPPR4", 
                                              stromal_signature == "TXNDC3" ~ "NME8",
                                              stromal_signature == "ODZ4" ~ "TENM4",
                                              TRUE ~ stromal_signature),
                immune_signature = case_when(immune_signature == "FYB" ~ "FYB1",
                                             TRUE ~ immune_signature))
rownames(gene_sets) <- NULL

usethis::use_data(gene_sets, overwrite = TRUE)
