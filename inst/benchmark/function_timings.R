# Benchmark script used to populate vignettes/articles/timing.Rmd.
# Run with: Rscript inst/benchmark/function_timings.R
# Writes a CSV to inst/benchmark/function_timings.csv.

suppressMessages({
  devtools::load_all(quiet = TRUE)
  library(phyloseq)
})

set.seed(1)

data(data_fungi)
data(data_fungi_mini)

# Drop NA-height samples for functions that need a complete grouping factor.
data_fungi_woNA <- subset_samples(
  data_fungi,
  !is.na(data_fungi@sam_data$Height)
)

# Number of replicate timings per function (low to keep the benchmark feasible)
n_reps <- 1

# Helper: time `expr` and return elapsed seconds (3rd value from system.time).
time_one <- function(label, dataset, expr_fun) {
  res <- tryCatch(
    {
      gc(verbose = FALSE)
      elapsed_t <- system.time(
        replicate(n_reps, expr_fun(), simplify = FALSE)
      )
      unname(elapsed_t["elapsed"]) / n_reps
    },
    error = function(e) {
      message(label, " (", dataset, ") failed: ", conditionMessage(e))
      NA_real_
    }
  )
  message(sprintf("%-32s %-20s %7.2f s", label, dataset, res))
  data.frame(
    fonction = label,
    dataset = dataset,
    elapsed_s = round(res, 2),
    stringsAsFactors = FALSE
  )
}

rows <- list()

# ---- Alpha diversity ----
rows <- c(
  rows,
  list(time_one("hill_pq", "data_fungi_mini", function() {
    hill_pq(data_fungi_mini, fact = "Height")
  })),
  list(time_one("hill_pq", "data_fungi", function() {
    hill_pq(data_fungi_woNA, fact = "Height")
  })),
  list(time_one("hill_tuckey_pq", "data_fungi_mini", function() {
    hill_tuckey_pq(data_fungi_mini, "Height")
  })),
  list(time_one("hill_tuckey_pq", "data_fungi", function() {
    hill_tuckey_pq(data_fungi_woNA, "Height")
  })),
  list(time_one("profile_hill_pq", "data_fungi_mini", function() {
    profile_hill_pq(data_fungi_mini)
  })),
  list(time_one("profile_hill_pq", "data_fungi", function() {
    profile_hill_pq(data_fungi_woNA)
  })),
  list(time_one("hill_acc_pq[sample,n=10]", "data_fungi_mini", function() {
    hill_acc_pq(data_fungi_mini, type = "sample", n_permutations = 10)
  })),
  list(time_one("hill_acc_pq[sample,n=10]", "data_fungi", function() {
    hill_acc_pq(data_fungi_woNA, type = "sample", n_permutations = 10)
  }))
)

# ---- Beta diversity / multivariate tests ----
rows <- c(
  rows,
  list(time_one("adonis_pq", "data_fungi_mini", function() {
    adonis_pq(data_fungi_mini, "Height")
  })),
  list(time_one("adonis_pq", "data_fungi", function() {
    adonis_pq(data_fungi_woNA, "Height")
  })),
  list(time_one("graph_test_pq", "data_fungi_mini", function() {
    graph_test_pq(data_fungi_mini, "Height")
  })),
  list(time_one("graph_test_pq", "data_fungi", function() {
    graph_test_pq(data_fungi_woNA, "Height")
  })),
  list(time_one("plot_tsne_pq", "data_fungi_mini", function() {
    plot_tsne_pq(data_fungi_mini, fact = "Height")
  })),
  list(time_one("plot_tsne_pq", "data_fungi", function() {
    plot_tsne_pq(data_fungi_woNA, fact = "Height")
  }))
)

# ---- Cleanup / verification ----
rows <- c(
  rows,
  list(time_one("verify_pq", "data_fungi_mini", function() {
    verify_pq(data_fungi_mini, check_taxonomy = FALSE)
  })),
  list(time_one("verify_pq", "data_fungi", function() {
    verify_pq(data_fungi, check_taxonomy = FALSE)
  })),
  list(time_one("verify_tax_table", "data_fungi_mini", function() {
    suppressWarnings(verify_tax_table(data_fungi_mini))
  })),
  list(time_one("verify_tax_table", "data_fungi", function() {
    suppressWarnings(verify_tax_table(data_fungi))
  }))
)

# ---- Plotting ----
rows <- c(
  rows,
  list(time_one("summary_plot_pq", "data_fungi_mini", function() {
    summary_plot_pq(data_fungi_mini)
  })),
  list(time_one("summary_plot_pq", "data_fungi", function() {
    summary_plot_pq(data_fungi)
  })),
  list(time_one("ggvenn_pq", "data_fungi_mini", function() {
    ggvenn_pq(data_fungi_mini, fact = "Height")
  })),
  list(time_one("ggvenn_pq", "data_fungi", function() {
    ggvenn_pq(data_fungi_woNA, fact = "Height")
  }))
)

out <- do.call(rbind, rows)
out_path <- file.path("inst", "benchmark", "function_timings.csv")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_path, row.names = FALSE)
message("Wrote ", out_path)
