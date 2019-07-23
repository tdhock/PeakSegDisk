\name{PeakSegFPOP_file}
\alias{PeakSegFPOP_file}
\title{PeakSegFPOP using disk storage}
\description{Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
using a file on disk to store the O(N) function piece lists,
each of size O(log N).
This is a low-level function that just runs the algo
and produces the result files (without reading them into R),
so normal users are recommended to instead use PeakSegFPOP_dir,
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
output files, which are in the same directory as bedGraph.file,
and named using the penalty value,
e.g. coverage.bedGraph_penalty=136500650856.439_loss.tsv }
}

\value{A list of input parameters (bedGraph.file, penalty) and result
files (segments, db, loss).}

\author{Toby Dylan Hocking}




\examples{

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
pstr <- "10.5"
PeakSegFPOP_file(tmp, pstr)

outf <- function(suffix){
  paste0(tmp, "_penalty=", pstr, "_", suffix)
}
segments.bed <- outf("segments.bed")
seg.df <- read.table(segments.bed)
names(seg.df) <- col.name.list$segments
seg.df

loss.tsv <- outf("loss.tsv")
loss.df <- read.table(loss.tsv)
names(loss.df) <- col.name.list$loss
loss.df

}
