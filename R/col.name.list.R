### Named list of character vectors (column names of bed/bedGraph/tsv
### files).
col.name.list <- list(
  loss=c(
    "penalty", "segments", "peaks", "bases", "bedGraph.lines",
    "mean.pen.cost", "total.loss", "equality.constraints",
    "mean.intervals", "max.intervals"),
  segments=c("chrom","chromStart", "chromEnd", "status", "mean"),
  coverage=c("chrom", "chromStart", "chromEnd", "count"))

