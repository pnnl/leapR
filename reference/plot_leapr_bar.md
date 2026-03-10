# Plot leapR pathway bars (single panel)

This plotting helper expects leapR generated results to plot. It will
use BH_pvalue if present, otherwise pvalue.

## Usage

``` r
plot_leapr_bar(
  res_df,
  title = NULL,
  top_n = 15,
  star_thresholds = c(0.05, 0.01, 0.001),
  wrap = 42,
  max_stars = 5L,
  fill_sig = "#2C7BB6",
  fill_ns = "#BFD7FF",
  outline = NA,
  axis_text_y_size = 8,
  axis_text_x_size = 9
)
```

## Arguments

- res_df:

  A leapR df containing BH_pvalue (or pvalue) and a pathway/term label
  column.

- title:

  Plot title.

- top_n:

  Number of top pathways/genes to display.

- star_thresholds:

  list of numeric significance thresholds for star annotations.

- wrap:

  Wrap width for pathway labels (helps formatting).

- max_stars:

  Maximum number of stars to draw per bar (default 5).

- fill_sig:

  Fill color for significant bars

- fill_ns:

  Fill color for non-significant bars.

- outline:

  Bar border color.

- axis_text_y_size:

  Font size for y-axis (category) labels.

- axis_text_x_size:

  Font size for x-axis (numeric) labels.

## Value

A ggplot2 object (or `NULL` if nothing to plot).
