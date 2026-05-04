#!/usr/bin/env Rscript
# Compute RF distances on a .trees file and output the upper-triangle values.
#
# Usage:
#   Rscript phangorn.R <trees_file> <metric> <output_file>
#
# Reads ALL trees from the file, computes the specified metric, and writes:
#   - Lines 2+: upper-triangle distances, one per line, row-major order
#     (i.e. d[1,2], d[1,3], ..., d[1,n], d[2,3], ..., d[n-1,n])
#
# Progress is printed to stderr.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript phangorn.R <trees_file> <metric> <output_file>\n",
      file = stderr())
  quit(status = 1)
}

trees_file <- args[1]
metric <- args[2]
output_file <- args[3]

suppressPackageStartupMessages({
  library(phangorn)
  library(ape)
})

# --- Load ---
cat(sprintf("Reading %s ...\n", trees_file), file = stderr())
t0 <- proc.time()
trees <- read.tree(trees_file)
load_time <- (proc.time() - t0)[3]
n <- length(trees)
n_tips <- length(trees[[1]]$tip.label)
cat(sprintf("Read %d trees (%d tips) in %.2fs\n", n, n_tips, load_time),
    file = stderr())

# --- Compute Distance ---
cat(sprintf("Computing %s on %d trees ...\n", metric, n), file = stderr())
t0 <- proc.time()
if (metric == "RF") {
  d <- RF.dist(trees)
} else if (metric == "wRF") {
    d <- wRF.dist(trees)
} else if (metric == "KF") {
    d <- KF.dist(trees)
} else {
  cat(sprintf("Unsupported metric: %s\n", metric), file = stderr())
  quit(status = 1)
}

# --- Write output ---
# Upper triangle in row-major order: (1,2), (1,3), ..., (1,n), (2,3), ...
m <- as.matrix(d)
idx <- which(upper.tri(m), arr.ind = TRUE)
idx <- idx[order(idx[, 1], idx[, 2]), ]
vals <- m[idx]

con <- file(output_file, "w")
writeLines(formatC(vals, digits = 17, format = "fg"), con)
close(con)

cat(sprintf("Wrote %d distances to %s\n", length(vals), output_file),
    file = stderr())