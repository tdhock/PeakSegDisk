PeakSegFPOP_df <- structure(function
### Write data frame to disk then run PeakSegFPOP solver.
(count.df, 
### data.frame with columns count, chromStart, chromEnd.
  pen.num,
### Non-negative numeric scalar.
  base.dir=tempdir()
){
  if(!(
    is.numeric(pen.num) &&
    length(pen.num)==1 &&
    0 <= pen.num)){
    stop("pen.num must be non-negative numeric scalar")
  }
  if(!is.data.frame(count.df)){
    stop("count.df must be data.frame")
  }
  exp.names <- c("chrom", "chromStart", "chromEnd", "count")
  if(!identical(names(count.df), exp.names)){
    stop("count.df must have names ", paste(exp.names, collapse=", "))
  }
  n.data <- nrow(count.df)
  if(n.data < 3){
    stop("must have at least three rows in count.df")
  }
  if(!is.integer(count.df$chromStart)){
    stop("count.df$chromStart must be integer")
  }
  if(!is.integer(count.df$chromEnd)){
    stop("count.df$chromEnd must be integer")
  }
  if(!is.integer(count.df$count)){
    stop("count.df$count must be integer")
  }
  if(!all(count.df$chromStart < count.df$chromEnd)){
    stop("chromStart must be less than chromEnd for all rows of count.df")
  }
  if(any(count.df$chromStart < 0)){
    stop("count.df$chromStart must always be non-negative")
  }
  if(any(count.df$count < 0)){
    stop("count.df$count must always be non-negative")
  }
  data.dir <- file.path(
    base.dir,
    with(count.df, sprintf(
      "%s-%d-%d", chrom[1], min(chromStart), max(chromEnd))))
  dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
  coverage.bedGraph <- file.path(data.dir, "coverage.bedGraph")
  write.table(
    count.df, coverage.bedGraph,
    sep="\t", row.names=FALSE, col.names=FALSE)
  problem.PeakSegFPOP(data.dir, paste(pen.num))
### List of solver results, same as problem.PeakSegFPOP.
}, ex=function(){

  library(PeakSegDisk)
  sim.seg <- function(seg.mean, size.mean=15){
    seg.size <- rpois(1, size.mean)
    rpois(seg.size, seg.mean)
  }
  set.seed(1)
  seg.mean.vec <- c(1.5, 3.5, 0.5, 4.5, 2.5)
  z.list <- lapply(seg.mean.vec, sim.seg)
  z.rep.vec <- unlist(z.list)
  library(ggplot2)
  count.df <- data.frame(
    position=seq_along(z.rep.vec),
    count=z.rep.vec)
  gg.count <- ggplot()+
    geom_point(aes(
      position, count),
      shape=1,
      data=count.df)
  gg.count
  
  n.segs <- length(seg.mean.vec)
  seg.size.vec <- sapply(z.list, length)
  seg.end.vec <- cumsum(seg.size.vec)
  change.vec <- seg.end.vec[-n.segs]+0.5
  change.df <- data.frame(
    changepoint=change.vec)
  gg.change <- gg.count+
    geom_vline(aes(
      xintercept=changepoint),
      data=change.df)
  gg.change
  
  z.rle.vec <- rle(z.rep.vec)
  chromEnd <- cumsum(z.rle.vec$lengths)
  coverage.df <- data.frame(
    chrom="chrUnknown",
    chromStart=c(0L, chromEnd[-length(chromEnd)]),
    chromEnd,
    count=z.rle.vec$values)
  gg.rle <- gg.change+
    geom_segment(aes(
      chromStart+0.5, count, xend=chromEnd+0.5, yend=count),
      data=coverage.df)
  gg.rle
  
  fit <- PeakSegFPOP_df(coverage.df, 10.5)
  gg.rle+
    geom_segment(aes(
      chromStart+0.5, mean, xend=chromEnd+0.5, yend=mean),
      color="green",
      data=fit$segments)
  
})
  
PeakSegFPOP_vec <- structure(function
### Convert rle to integer vector, write to disk, then run PeakSegFPOP solver.
(count.vec, 
### data.frame with columns count, chromStart, chromEnd.
  pen.num,
### Non-negative numeric scalar.
  base.dir=tempdir()
){
  if(!(
    is.numeric(pen.num) &&
    length(pen.num)==1 &&
    0 <= pen.num)){
    stop("pen.num must be non-negative numeric scalar")
  }
  if(!is.integer(count.vec)){
    stop("count.vec must be integer")
  }
  z.rle.vec <- rle(count.vec)
  chromEnd <- cumsum(z.rle.vec$lengths)
  coverage.df <- data.frame(
    chrom="chrUnknown",
    chromStart=c(0L, chromEnd[-length(chromEnd)]),
    chromEnd,
    count=z.rle.vec$values)
  PeakSegFPOP_df(coverage.df, 10.5)
### List of solver results, same as problem.PeakSegFPOP.
}, ex=function(){

  library(PeakSegDisk)
  sim.seg <- function(seg.mean, size.mean=15){
    seg.size <- rpois(1, size.mean)
    rpois(seg.size, seg.mean)
  }
  set.seed(1)
  seg.mean.vec <- c(1.5, 3.5, 0.5, 4.5, 2.5)
  z.list <- lapply(seg.mean.vec, sim.seg)
  z.rep.vec <- unlist(z.list)
  library(ggplot2)
  count.df <- data.frame(
    position=seq_along(z.rep.vec),
    count=z.rep.vec)
  gg.count <- ggplot()+
    geom_point(aes(
      position, count),
      shape=1,
      data=count.df)
  gg.count
  
  n.segs <- length(seg.mean.vec)
  seg.size.vec <- sapply(z.list, length)
  seg.end.vec <- cumsum(seg.size.vec)
  change.vec <- seg.end.vec[-n.segs]+0.5
  change.df <- data.frame(
    changepoint=change.vec)
  gg.change <- gg.count+
    geom_vline(aes(
      xintercept=changepoint),
      data=change.df)
  gg.change
  
  fit <- PeakSegFPOP_vec(z.rep.vec, 10.5)
  gg.change+
    geom_segment(aes(
      chromStart+0.5, mean, xend=chromEnd+0.5, yend=mean),
      color="green",
      data=fit$segments)
  
})
  
