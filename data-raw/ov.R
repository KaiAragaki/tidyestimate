# Prepares "ov" dataset ----
library(dplyr)
library(utils)

utils::download.file("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/estimate/inst/extdata/sample_input.txt?root=estimate",
              destfile = "./data-raw/input_data_old.txt", 
              method = "curl")
ov <- read.delim("./data-raw/input_data_old.txt") |>
        apply(2, round, digit = 4) # To reduce file size.
        # I will show that estimate and tidyestimate still produce the same
        # results later.

usethis::use_data(ov, overwrite = TRUE)
