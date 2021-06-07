library(dplyr)
library(tidyr)
library(stringr)

a <- get_identifiers(common_genes)

dupes_ent <- a$entrezgene_id[duplicated(a$entrezgene_id)] 

ent_dupes <- a[a$entrezgene_id %in% dupes_ent,]

dupes_hgnc <- a$hgnc_symbol[duplicated(a$hgnc_symbol)] 

blank <- filter(a, hgnc_symbol == "")

# Blanks that do not have another entry with a name.
anti_join(blank, ent_dupes)

# There are just two: 23766, and 79857. The first is a pseudogene, the second is
# an ncRNA gene. Both feel safe to remove.

a <- filter(a, hgnc_symbol != "")

dupes_hgnc <- a$hgnc_symbol[duplicated(a$hgnc_symbol)] 
# No more HGNC duplicates

dupes_ent <- a$entrezgene_id[duplicated(a$entrezgene_id)] 

ent_dupes <- a[a$entrezgene_id %in% dupes_ent,]

# Three duplicates remain. In all cases, one of the two is a two-symbol gene.
# I'm going to filter FOR the dashed genes, then use a filtering join to get 
# rid of them in the main dataset.

remove <- filter(ent_dupes, str_detect(hgnc_symbol, "-"))

a <- anti_join(a, remove)

dupes_ent <- a$entrezgene_id[duplicated(a$entrezgene_id)] 
# No more dupes!

# Now, who's missing?

aj <- anti_join(common_genes, a_fin, by = c("EntrezID" = "entrezgene_id"))
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

# this leaves 9 genes that are, for reasons unknown, not looked up by biomaRt.

# 9/10412 = ~0.00086 = ~0.08%. I'd say more that good enough.