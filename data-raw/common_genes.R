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
                attributes = c('entrezgene_id', 'hgnc_symbol'),
                filters = 'entrezgene_id',
                values = x$EntrezID, 
                mart = bm)
}

# Takes the common genes list used in the estimate package
common_genes <- get_identifiers(common_genes_old)

dupes_ent <- common_genes$entrezgene_id[duplicated(common_genes$entrezgene_id)] 
ent_dupes <- common_genes[common_genes$entrezgene_id %in% dupes_ent,]

dupes_hgnc <- common_genes$hgnc_symbol[duplicated(common_genes$hgnc_symbol)] 

blank <- filter(common_genes, hgnc_symbol == "")

# Blanks that do not have another entry with a name.
anti_join(blank, ent_dupes)

# There are just two: 23766, and 79857. The first is a pseudogene, the second is
# an ncRNA gene. Both feel safe to remove.

common_genes <- filter(common_genes, hgnc_symbol != "")

dupes_hgnc <- common_genes$hgnc_symbol[duplicated(common_genes$hgnc_symbol)] 
# No more HGNC duplicates

dupes_ent <- common_genes$entrezgene_id[duplicated(common_genes$entrezgene_id)] 
ent_dupes <- common_genes[common_genes$entrezgene_id %in% dupes_ent,]

# Three duplicates remain. In all cases, one of the two is a two-symbol gene.
# I'm going to filter FOR the dashed genes, then use a filtering join to get 
# rid of them in the main dataset.

remove <- filter(ent_dupes, str_detect(hgnc_symbol, "-"))

common_genes <- anti_join(common_genes, remove)

dupes_ent <- common_genes$entrezgene_id[duplicated(common_genes$entrezgene_id)] 
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

usethis::use_data(common_genes, overwrite = TRUE)
