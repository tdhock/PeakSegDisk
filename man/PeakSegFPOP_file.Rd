\name{PeakSegFPOP_file}
\alias{PeakSegFPOP_file}
\title{PeakSegFPOP using disk storage}
\description{Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
using a file on disk to store the O(N) function piece lists,
each of size O(log N).
This is a low-level function that just runs the algo
and produces the result files (without reading them into R),
so normal users are recommended to instead use \code{\link{PeakSegFPOP_dir}},
which calls this function then reads the result files into R.}
\usage{PeakSegFPOP_file(bedGraph.file, pen.str)}
\arguments{
  \item{bedGraph.file}{character scalar: tab-delimited tabular text file with four
columns: chrom, chromStart, chromEnd, coverage. The algorithm
creates a large temporary file in the same directory, so make sure
that there is disk space available on that device.}
  \item{pen.str}{character scalar that can be converted to a numeric scalar via
as.numeric: non-negative penalty. More penalty means fewer
peaks. "0" and "Inf" are OK. Character is required rather than
numeric, so that the user can reliably find the results in the
output files, which are in the same directory as \code{bedGraph.file},
and named using the penalty value,
e.g. coverage.bedGraph_penalty=136500650856.439_loss.tsv }
}

\value{A list of input parameters (\code{bedGraph.file}, penalty) and result
files (segments, db, loss).}

\author{Toby Dylan Hocking}




\examples{

r <- function(chrom, chromStart, chromEnd, coverage){
  data.frame(chrom, chromStart, chromEnd, coverage)
}
four <- rbind(
  r("chr1", 0, 10,  2),
  r("chr1", 10, 20, 10),
  r("chr1", 20, 30, 14),
  r("chr1", 30, 40, 13))
dir.create(prob.dir <- tempfile())
coverage.bedGraph <- file.path(prob.dir, "coverage.bedGraph")
write.table(
  four, coverage.bedGraph,
  sep="\t", row.names=FALSE, col.names=FALSE)
pstr <- "10.5"
result.list <- PeakSegDisk::PeakSegFPOP_file(coverage.bedGraph, pstr)
dir(prob.dir)

## segments file can be read to see optimal segment means.
outf <- function(suffix){
  paste0(coverage.bedGraph, "_penalty=", pstr, suffix)
}
segments.bed <- outf("_segments.bed")
seg.df <- read.table(segments.bed)
names(seg.df) <- col.name.list$segments
seg.df

## loss file can be read to see optimal Poisson loss, etc.
loss.tsv <- outf("_loss.tsv")
loss.df <- read.table(loss.tsv)
names(loss.df) <- col.name.list$loss
loss.df

## The temporary db file is used in the underlying C++ code to store
## the optimal cost function for each candidate changepoint. It can
## be safely be deleted.
tmp.db <- outf(".db")
unlink(tmp.db)
dir(prob.dir)

}
