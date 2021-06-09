#' Plot Affymetrix purity scores against ESTIMATE study purity scores
#'
#' @param scores a data frame, usually one output from `estimate_score`
#' @param is_affymetrix logical. Are these data from an Affymetrix experiment?
#'   Must be TRUE - this is essentially a verification from the user
#'
#' @return a ggplot
#' @export
#'
#' @importFrom ggplot2 aes
#'
#' @examples
#' filter_common_genes(ov, id = "hgnc_symbol", tidy = FALSE) |> 
#'   estimate_score(is_affymetrix = TRUE) |>
#'   plot_purity(is_affymetrix = TRUE)
#' 
plot_purity <- function(scores, is_affymetrix) {
        if(!is_affymetrix){
                stop("Non-Affymetrix data is not supported")
        }
        
        scores_nona <- dplyr::filter(scores, !is.na(.data$purity))
        
        global_min <- min(scores_nona$estimate, tidyestimate::purity_data_affy$estimate, na.rm = T)
        global_max <- max(scores_nona$estimate, tidyestimate::purity_data_affy$estimate, na.rm = T)
        purity_line <- dplyr::tibble(estimate = seq(from = global_min, 
                                                    to = global_max, 
                                                    by = 10),
                                     purity = cos(0.6049872018 + 0.0001467884 * .data$estimate))
        
        ggplot2::ggplot(tidyestimate::purity_data_affy, aes(x = .data$estimate)) + 
                ggplot2::geom_point(aes(y = .data$purity_observed), alpha = 0.2, shape = 1) +
                ggplot2::geom_line(aes(y = .data$ci_95_low), linetype = "dashed", alpha = 0.3) +
                ggplot2::geom_line(aes(y = .data$ci_95_high), linetype = "dashed", alpha = 0.3) +
                ggplot2::geom_line(data = purity_line, aes(y = .data$purity), alpha = 0.4) +
                ggplot2::geom_point(data = scores, aes(y = .data$purity), color = "red", size = 2, shape = 1) +
                ggrepel::geom_text_repel(data = scores, aes(y = .data$purity, label = .data$sample)) + 
                ggplot2::coord_cartesian(ylim = c(0, 1)) +
                ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA, 
                                                                        color = "black"), 
                               panel.grid = ggplot2::element_blank()) +
                ggplot2::labs(x = "ESTIMATE score", y = "Tumor purity")
}
