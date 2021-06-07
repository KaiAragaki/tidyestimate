library(dplyr)

al2 <- aliases |> 
        group_by(entrezgene_id, hgnc_symbol) |> 
        summarize(syn = list(external_synonym)) |> 
        rowwise() |> 
        mutate(all = list(c(syn, hgnc_symbol))) |> 
        select(entrezgene_id, all) |> 
        unnest(cols = c(all)) |> 
        distinct()