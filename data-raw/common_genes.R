# Prepares "common_genes" dataset ------
library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(utils)

utils::download.file("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/estimate/data/common_genes.RData?revision=2&root=estimate",
              destfile = "./data-raw/common_genes_old.RData", 
              method = "curl")
load("./data-raw/common_genes_old.RData")

common_genes_old <- common_genes

get_identifiers <- function(x) {
        bm <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                               dataset = "hsapiens_gene_ensembl")
        biomaRt::getBM(
                attributes = c("entrezgene_id", "hgnc_symbol", "external_synonym"),
                filters = "entrezgene_id",
                values = x$EntrezID, 
                mart = bm)
}

# Takes the common genes list used in the estimate package
common_genes <- get_identifiers(common_genes_old)

dupe_check <- unique(common_genes[,1:2])

dupes_ent <- dupe_check$entrezgene_id[duplicated(dupe_check$entrezgene_id)] 
ent_dupes <- dupe_check[dupe_check$entrezgene_id %in% dupes_ent,]

dupes_hgnc <- dupe_check$hgnc_symbol[duplicated(dupe_check$hgnc_symbol)] 
# Only dupes are those that have a blank name.
# Getting rid of them

common_genes <- filter(common_genes, hgnc_symbol != "")

dupe_check <- unique(common_genes[,1:2])
dupes_hgnc <- dupe_check$hgnc_symbol[duplicated(dupe_check$hgnc_symbol)] 
# No more HGNC duplicates

dupes_ent <- dupe_check$entrezgene_id[duplicated(dupe_check$entrezgene_id)] 
ent_dupes <- dupe_check[dupe_check$entrezgene_id %in% dupes_ent,]

# Three duplicates remain. In all cases, one of the two is a two-symbol gene.
# I'm going to filter FOR the dashed genes, then use a filtering join to get 
# rid of them in the main dataset.

remove <- filter(ent_dupes, str_detect(hgnc_symbol, "-"))

common_genes <- anti_join(common_genes, remove)

dupe_check <- unique(common_genes[,1:2])
dupes_ent <- dupe_check$entrezgene_id[duplicated(dupe_check$entrezgene_id)] 
# No more dupes!

# Now, who's missing?

aj <- anti_join(common_genes_old, common_genes, by = c("EntrezID" = "entrezgene_id"))
# 4104   - pseudogene
# 5003   - ??
# 8123   - ??
# 10896  - discontinued
# 23285  - discontinued
# 23766  - pseudogene
# 27004  - ??
# 28363  - pseudogene
# 50514  - ??
# 54581  - pseudogene
# 54661  - pseudogene
# 55000  - ??
# 65012  - ??
# 79857  - uncharacterized LOC
# 80094  - ??
# 80761  - discontinued
# 83955  - pseudogene
# 91355  - ??
# 117153 - replaced by 4253
# 255926 - pseudogene
# 643314 - ??

# So, 12 out of 21 genes are justifiably removed
# This leaves 9 genes that are, for reasons unknown, not looked up by biomaRt.
# 9/10412 = ~0.00086

common_genes <- mutate(common_genes, across(everything(), as.character))

# It is very important that gene set genes find a match. Since the find alias
# function is rather conservative, we need to 'force' it to make matches by
# removing extraneous aliases that are probably not going to be found in the
# given dataframe

# At the time of writing, the stromal_signature has 5 HGNC symbols that are out
# of date, and only 2 have multiple aliases. Additionally, the immune_signature
# has 1 symbols that is out of date, with multiple aliases.

# These will be updated in the genesets, but their aliases also need tuning to
# ensure a proper match can be made.

# We are essentially forcing a "one to one" relationship with the most common
# alias and the proper HGNC symbol

common_genes <- filter(common_genes, (hgnc_symbol != "CCN4") | (hgnc_symbol == "CCN4" & external_synonym == "WISP1"),
                       (hgnc_symbol != "PLPPR4") | (hgnc_symbol == "PLPPR4" & external_synonym == "LPPR4"),
                       (hgnc_symbol != "FYB1") | (hgnc_symbol == "FYB1" & external_synonym == "FYB"))

common_genes <- mutate(common_genes,
                       external_synonym = stringi::stri_trans_general(external_synonym, "latin-ascii"))

usethis::use_data(common_genes, overwrite = TRUE)
