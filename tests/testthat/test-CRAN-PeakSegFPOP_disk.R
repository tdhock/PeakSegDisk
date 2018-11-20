library(testthat)
context("PeakSegFPOP_disk")
library(PeakSegDisk)

r <- function(chrom, chromStart, chromEnd, coverage){
  data.frame(chrom, chromStart, chromEnd, coverage)
}
four <- rbind(
  r("chr1", 0, 10,  2),
  r("chr1", 10, 20, 10),
  r("chr1", 20, 30, 14),
  r("chr1", 30, 40, 13))
write.table(
  four, tmp <- tempfile(),
  sep="\t", row.names=FALSE, col.names=FALSE)

test_that("character penalty works", {
  names.list <- PeakSegFPOP_disk(tmp, "10.5")
  unlink(names.list$db)
  seg.df <- read.table(names.list$segments)
  names(seg.df) <- c("chrom", "chromStart", "chromEnd", "status", "mean")
  expect_equal(paste(seg.df$chrom), c("chr1", "chr1", "chr1"))
  expect_equal(seg.df$chromStart, c(30, 10, 0))
  expect_equal(seg.df$chromEnd, c(40, 30, 10))
  expect_equal(paste(seg.df$status), c("background", "peak", "background"))
  m <- mean(four$coverage[-1])
  expected.mean <- c(m, m, four$coverage[1])
  expect_equal(seg.df$mean, expected.mean, tol=1e-3)
  loss.df <- read.table(names.list$loss)
  names(loss.df) <- c(
    "penalty", "segments", "peaks", "bases",
    "mean.pen.cost", "total.loss", "equality.constraints",
    "mean.intervals", "max.intervals")
})

test_that("error for numeric penalty", {
  expect_error({
    PeakSegFPOP_disk(tmp, 10.5)
  }, "pen.str must be a character string that can be converted to a non-negative numeric scalar")
})

data(Mono27ac, envir=environment())
data.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11:60000-580000")
dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
write.table(
  Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
pen.str <- "1952.6"
lock.src <- system.file("extdata", "__db.lock.file", package="PeakSegDisk", mustWork=TRUE)
lock.dest <- file.path(data.dir, paste0("__db.coverage.bedGraph_penalty=", pen.str, ".db"))
file.copy(lock.src, lock.dest)
coverage.bedGraph <- file.path(data.dir, "coverage.bedGraph")

test_that("error for lock file", {
  expect_error({
    PeakSegFPOP_disk(coverage.bedGraph, pen.str)
  }, lock.dest, fixed=TRUE)
})
