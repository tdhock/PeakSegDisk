library(testthat)
context("PeakSegFPOP_disk")
library(PeakSegDisk)

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
test_that("PeakSegFPOP_disk errors for lock file", {
  expect_error({
    PeakSegFPOP_disk(coverage.bedGraph, pen.str)
  }, "database locked; delete lock file __db.*.db to unlock", fixed=TRUE)
})

test_that("problem.PeakSegFPOP deletes lock file", {
  fit <- problem.PeakSegFPOP(data.dir, pen.str)
  expect_is(fit, "list")
  expect_is(fit$loss, "data.table")
  expect_is(fit$segments, "data.table")
})
