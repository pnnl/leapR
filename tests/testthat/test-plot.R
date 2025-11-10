# Example data
set.seed(1)
n <- 20
df <- data.frame(
feature    = paste0("Path_", seq_len(n)),
BH_pvalue  = runif(n, min = 1e-4, max = 0.2),
stringsAsFactors = FALSE
)

# Plot basic
p <- plot_leapr_bar(df, title = "Test", top_n = 10)
expect_s3_class(p, "ggplot")

# Check that its a ggplot object
# Check that the base layer has exactly 10( top_n) rows
gb <- ggplot2::ggplot_build(p)
expect_gte(length(gb$data), 2)
base_n <- nrow(gb$data[[1]])
expect_equal(base_n, 10)
