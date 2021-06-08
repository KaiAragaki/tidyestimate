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
                      immune_signature = ImmuneSignature)
rownames(gene_sets) <- NULL

usethis::use_data(gene_sets, overwrite = TRUE)
