#' @title Plot leapR pathway bars (single panel)
#' @description
#' This plotting helper expects leapR generated results to plot.
#' It will use BH_pvalue if present, otherwise pvalue.
#'
# REQUIRED ARGUMENTS
#' @param res_df A leapR df containing BH_pvalue (or pvalue) and a pathway/term label column.
#' @param title Plot title.
# OPTIONAL ARGUMENTS 
#' @param top_n Number of top pathways/genes to display.
#' @param star_thresholds list of numeric significance thresholds for star annotations.
#' @param wrap Wrap width for pathway labels (helps formatting).
#' @param max_stars Maximum number of stars to draw per bar (default 5).
#' @param fill_sig Fill color for significant bars
#' @param fill_ns Fill color for non-significant bars.
#' @param outline Bar border color.
#' @param axis_text_y_size Font size for y-axis (category) labels.
#' @param axis_text_x_size  Font size for x-axis (numeric) labels.
#'
#' @return A \pkg{ggplot2} object (or \code{NULL} if nothing to plot).
#' @export
plot_leapr_bar <- function(res_df,
                           title = NULL,
                           top_n = 15,
                           star_thresholds = c(0.05, 0.01, 0.001),
                           wrap = 42,
                           max_stars = 5L,
                           # colors
                           fill_sig = "#2C7BB6",   # darker for significant
                           fill_ns  = "#BFD7FF",   # lighter for non-significant
                           outline  = NA,
                           # text sizes
                           axis_text_y_size = 8,
                           axis_text_x_size = 9) {
  df <- .prep_bar_df(
    res_df           = res_df,
    top_n            = top_n,
    star_thresholds  = star_thresholds,
    wrap             = wrap,
    max_stars        = max_stars
  )
  if (is.null(df)) return(NULL)

  df_sig <- df[df$signif %in% TRUE, , drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(feature, score), y = score)) +
    ggplot2::geom_col(fill = fill_ns, width = 0.85, color = outline) +
    ggplot2::geom_col(data = df_sig, fill = fill_sig, width = 0.85, color = outline) +
    ggplot2::geom_text(ggplot2::aes(y = score + 0.12, label = stars),
                       size = 3, hjust = 0) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = axis_text_y_size),
      axis.text.x        = ggplot2::element_text(size = axis_text_x_size),
      plot.margin        = ggplot2::margin(5.5, 30, 5.5, 5.5)
    ) +
    ggplot2::labs(title = title, x = NULL,
                  y = expression(-log[10]("BH-adjusted p")~"(FDR)"))
}

# =========================
# Internal helpers (not exported)
# =========================

#' @keywords internal
#' @noRd
.pick_pcol <- function(df) {
  n <- colnames(df)
  if ("BH_pvalue"       %in% n) return(list(col = "BH_pvalue"))
  if ("pvalue" %in% n) return(list(col = "pvalue"))
  NULL
}

#' @keywords internal
#' @noRd
.stars_from_p <- function(p, cutoffs, max_stars = 5L) {
  if (is.na(p)) return("")
  n <- sum(p < cutoffs)
  n <- max(0L, min(as.integer(n), as.integer(max_stars)))
  if (n == 0L) "" else paste(rep("*", n), collapse = "")
}

#' @keywords internal
#' @noRd
.prep_bar_df <- function(res_df,
                         top_n = 15,
                         star_thresholds = c(0.05, 0.01, 0.001),
                         wrap = 42,
                         max_stars = 5L) {
  if (is.null(res_df)) return(NULL)
  df <- as.data.frame(res_df)
  if (!nrow(df)) return(NULL)

  term_col <- intersect(
    c("feature","term","Term","pathway","Pathway","set","Set","geneset","gene_set","Category"),
    colnames(df)
  )
  if (length(term_col)) {
    df$feature <- as.character(df[[term_col[1]]])
  } else if (!("feature" %in% names(df))) {
    df <- tibble::rownames_to_column(df, "feature")
  }

  picked <- .pick_pcol(df)
  if (is.null(picked)) {
    message("No BH-adjusted p-value column found (expected BH_pvalue or pvalue from leapR).")
    return(NULL)
  }

  adjp <- suppressWarnings(as.numeric(df[[picked$col]]))
  star_thresholds <- sort(unique(as.numeric(star_thresholds)))
  star_thresholds <- star_thresholds[is.finite(star_thresholds) & star_thresholds > 0]
  if (!length(star_thresholds)) star_thresholds <- c(0.05, 0.01, 0.001)

  out <- df |>
    dplyr::mutate(
      adj_p   = adjp,
      score   = -log10(pmax(adj_p, 1e-300)),
      stars   = vapply(adjp, .stars_from_p, character(1),
                       cutoffs = star_thresholds, max_stars = max_stars),
      signif  = !is.na(adj_p) & is.finite(adj_p) & (adj_p < max(star_thresholds)),
      feature = stringr::str_wrap(as.character(feature), width = wrap)
    ) |>
    dplyr::filter(is.finite(adj_p)) |>
    dplyr::arrange(adj_p, dplyr::desc(score)) |>
    dplyr::slice_head(n = top_n)

  if (!nrow(out)) return(NULL)
  out
}
