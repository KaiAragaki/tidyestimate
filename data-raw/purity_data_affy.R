# Prepares "purity_data_affy" dataset ----
library(dplyr)
library(utils)

utils::download.file("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/estimate/data/PurityDataAffy.RData?root=estimate",
              destfile = "./data-raw/purity_data_affy_old.RData", 
              method = "curl")
load("./data-raw/purity_data_affy_old.RData")

# Rename cols for better homogeneity and descriptive value
purity_data_affy <- PurityDataAffy |>
        dplyr::rename(purity_observed = tumor.purity,
                      stromal = StromalScore,
                      immune = ImmuneScore,
                      estimate = ESTIMATEScore,
                      purity_predicted = fit,
                      ci_95_low = lwr.p,
                      ci_95_high  = upr.p
        )

usethis::use_data(purity_data_affy, overwrite = TRUE)
