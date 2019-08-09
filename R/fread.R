wc2int <- function
### Convert wc output to integer number of lines.
(wc.output
  ){
  no.initial.spaces <- sub("^ *", "", wc.output)
  lines.chr <- sub(" .*", "", no.initial.spaces)
  result <- as.integer(lines.chr)
  if(!(
    is.integer(result) &&
    length(result)==1 &&
    is.finite(result)
  )){
    print(wc.output)
    stop("could not extract line count")
  }
  result
### integer
}

fread.first <- structure(function # Quickly read first line
### Read the first line of a text file. Useful for quickly checking if
### the coverage.bedGraph_penalty=VALUE_segments.bed file is
### consistent with the coverage.bedGraph_penalty=VALUE_loss.tsv
### file. (used by the PeakSegFPOP_dir caching mechanism)
(file.name,
### Name of file to read.
  col.name.vec
### Character vector of column names.
){
  dt <- fread(file.name, nrows=1L, col.names=col.name.vec)
  dt
### Data table with one row.
}, ex=function(){

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
  first.seg.line <- fread.first(segments.bed, col.name.list$segments)
  last.seg.line <- fread.last(segments.bed, col.name.list$segments)

  loss.tsv <- outf("loss.tsv")
  loss.row <- fread.first(loss.tsv, col.name.list$loss)

  seg.bases <- first.seg.line$chromEnd - last.seg.line$chromStart
  loss.row$bases == seg.bases

})

fread.last <- structure(function # Quickly read last line
### Read the last line of a text file. Useful for quickly checking if
### the coverage.bedGraph_penalty=VALUE_segments.bed file is
### consistent with the coverage.bedGraph_penalty=VALUE_loss.tsv
### file. (used by the PeakSegFPOP_dir caching mechanism)
(file.name,
### Name of file to read.
  col.name.vec
### Character vector of column names.
){
  wc.cmd <- paste("wc -l", file.name)
  wc.output <- system(wc.cmd, intern=TRUE)
  lines.int <- wc2int(wc.output)
  dt <- fread(file.name, skip=lines.int-1L, col.names=col.name.vec)
  dt
### Data table with one row.
}, ex=function(){

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
  first.seg.line <- fread.first(segments.bed, col.name.list$segments)
  last.seg.line <- fread.last(segments.bed, col.name.list$segments)

  loss.tsv <- outf("loss.tsv")
  loss.row <- fread.first(loss.tsv, col.name.list$loss)

  seg.bases <- first.seg.line$chromEnd - last.seg.line$chromStart
  loss.row$bases == seg.bases

})

