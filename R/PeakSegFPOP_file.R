PeakSegFPOP_file <- structure(function # PeakSegFPOP using disk storage
### Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
### using a file on disk to store the O(N) function piece lists,
### each of size O(log N).
### This is a low-level function that just runs the algo
### and produces the result files (without reading them into R),
### so normal users are recommended to instead use PeakSegFPOP_dir,
### which calls this function then reads the result files into R.
(bedGraph.file,
### character scalar: tab-delimited tabular text file with four
### columns: chrom, chromStart, chromEnd, coverage. The algorithm
### creates a large temporary file in the same directory, so make sure
### that there is disk space available on that device.
  pen.str,
### character scalar that can be converted to a numeric scalar via
### as.numeric: non-negative penalty. More penalty means fewer
### peaks. "0" and "Inf" are OK. Character is required rather than
### numeric, so that the user can reliably find the results in the
### output files, which are in the same directory as bedGraph.file,
### and named using the penalty value,
### e.g. coverage.bedGraph_penalty=136500650856.439_loss.tsv
  db.file=NULL
### character scalar: file for writing temporary cost function
### database. Default NULL means to write the same disk where the
### bedGraph file is stored, but this may be useful use tempfile() if
### you don't want to perform the intensive disk i/o elsewhere.
){
  if(!(
    is.character(bedGraph.file) &&
    length(bedGraph.file)==1 &&
    file.exists(bedGraph.file)
  )){
    stop(
      "bedGraph.file=", bedGraph.file,
      " must be the name of a data file to segment")
  }
  if(!is.character(pen.str)){
    stop(paste(
      "pen.str must be a character string",
      "that can be converted to a non-negative numeric scalar"
    ))
  }
  penalty <- as.numeric(pen.str)
  if(!(
    is.numeric(penalty) &&
    length(penalty)==1 &&
    0 <= penalty && penalty <= Inf
  )){
    stop("as.numeric(pen.str)=", penalty, " but it must be a non-negative numeric scalar")
  }
  norm.file <- normalizePath(bedGraph.file, mustWork=TRUE)
  db.file <- if(is.null(db.file)){
    sprintf("%s_penalty=%s.db", norm.file, pen.str)
  }
  if(!(
    is.character(db.file) &&
    length(db.file)==1
  )){
    stop(
      "db.file=", db.file,
      " must be a temporary file name where cost function db can be written")
  }
  result <- .C(
    "PeakSegFPOP_interface",
    bedGraph.file=as.character(norm.file),
    penalty=pen.str,
    db.file=db.file,
    PACKAGE="PeakSegDisk")
  prefix <- paste0(bedGraph.file, "_penalty=", pen.str)
  loss.tsv <- paste0(prefix, "_loss.tsv")
  if(file.size(loss.tsv)==0){
    stop(
      "unable to write to loss output file ",
      loss.tsv,
      " (disk is probably full)"
    )
  }
  result
### A list of input parameters (bedGraph.file, penalty) and result
### files (segments, db, loss).
}, ex=function(){
  
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
  
})

