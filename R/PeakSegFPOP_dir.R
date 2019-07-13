PeakSegFPOP_dir <- structure(function
### Run PeakSegFPOP_file on one genomic segmentation problem
### directory, and read the result files into R. Actually, this
### function will first check if the result files are already present
### (and consistent), and if so, it will simply read them into R
### (without running PeakSegFPOP_file) -- this is a caching mechanism
### that can save a lot of time. 
(problem.dir,
### Path to a directory like sampleID/problems/problemID which
### contains a coverage.bedGraph file with the aligned read counts for
### one genomic segmentation problem.
 penalty.str
### character which can be interpreted as a non-negative numeric
### penalty parameter (larger values for fewer peaks). "0" means max
### peaks, "Inf" means no peaks. Needs to be a character because that
### is used to create files which cache/store the optimal solution.
 ){
   if(!(
     is.character(problem.dir) && 
     length(problem.dir)==1 &&
     dir.exists(problem.dir))){
     stop(
       "problem.dir=", problem.dir,
       " must be the name of a directory containing a file named coverage.bedGraph")
   }
  stopifnot(is.character(penalty.str))
  stopifnot(length(penalty.str)==1)
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  penalty_loss.tsv <- paste0(pre, "_loss.tsv")
  penalty_timing.tsv <- paste0(pre, "_timing.tsv")
  already.computed <- tryCatch({
    timing <- fread(
      penalty_timing.tsv,
      col.names=c("penalty", "megabytes", "seconds"))
    first.seg.line <- fread.first(penalty_segments.bed, col.name.list$segments)
    last.seg.line <- fread.last(penalty_segments.bed, col.name.list$segments)
    first.cov.line <- fread.first(prob.cov.bedGraph, col.name.list$coverage)
    last.cov.line <- fread.last(prob.cov.bedGraph, col.name.list$coverage)
    penalty.loss <- fread(penalty_loss.tsv, col.names=col.name.list$loss)
    nrow.ok <- nrow(timing)==1 && nrow(penalty.loss)==1 &&
      nrow(first.seg.line)==1 && nrow(last.seg.line)==1 &&
      nrow(first.cov.line)==1 && nrow(last.cov.line)==1
    loss.segments.consistent <-
      first.seg.line$chromEnd-last.seg.line$chromStart == penalty.loss$bases
    ## segments files are written by decoding/backtracking after
    ## dynamic progamming, so it is normal that the first line of the
    ## segment file is actually the last segment in terms of position
    ## on the chromosome.
    start.ok <- first.cov.line$chromStart == last.seg.line$chromStart
    end.ok <- last.cov.line$chromEnd == first.seg.line$chromEnd
    nrow.ok && loss.segments.consistent && start.ok && end.ok
  }, error=function(e){
    FALSE
  })
  if(!already.computed){
    penalty.db <- paste0(pre, ".db")
    unlink(penalty.db)#in case interrupted previously.
    base.db <- basename(penalty.db)
    base.lock <- paste0("__db.", base.db)
    path.lock <- file.path(problem.dir, base.lock)
    unlink(path.lock)
    seconds <- system.time({
      PeakSegFPOP_file(prob.cov.bedGraph, penalty.str)
    })[["elapsed"]]
    megabytes <- if(file.exists(penalty.db)){
      file.size(penalty.db)/1024/1024
    }else{
      0
    }
    timing <- data.table(
      penalty=as.numeric(penalty.str),
      megabytes,
      seconds)
    write.table(
      timing,
      penalty_timing.tsv,
      row.names=FALSE, col.names=FALSE,
      quote=FALSE, sep="\t")
    unlink(penalty.db)
    penalty.loss <- fread(penalty_loss.tsv, col.names=col.name.list$loss)
  }
  penalty.segs <- fread(penalty_segments.bed, col.names=col.name.list$segments)
  list(
    segments=penalty.segs,
    loss=data.table(
      penalty.loss,
      timing[, list(megabytes, seconds)]))
### Named list of two data.tables: segments has one row for every
### segment in the optimal model, and loss has one row and contains
### the following columns. penalty=same as input, segments=number of
### segments in optimal model, peaks=number of peaks in optimal model,
### bases=number of positions described in bedGraph file,
### bedGraph.lines=number of lines in bedGraph file, total.loss=total
### Poisson loss=sum_i
### m_i-z_i*log(m_i)=mean.pen.cost*bases-penalty*peaks,
### mean.pen.cost=mean penalized
### cost=(total.loss+penalty*peaks)/bases, equality.constraints=number
### of adjacent segment means that have equal values in the optimal
### solution, mean.intervals=mean number of intervals/candidate
### changepoints stored in optimal cost functions -- useful for
### characterizing the computational complexity of the algorithm,
### max.intervals=maximum number of intervals, megabytes=disk usage of
### *.db file, seconds=timing of PeakSegFPOP_file.
}, ex=function(){

  library(PeakSegDisk)
  data(Mono27ac, envir=environment())
  data.dir <- file.path(
    tempfile(),
    "H3K27ac-H3K4me3_TDHAM_BP",
    "samples",
    "Mono1_H3K27ac",
    "S001YW_NCMLS",
    "problems",
    "chr11-60000-580000")
  dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
  write.table(
    Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

  ## Compute one model with penalty=1952.6
  fit <- PeakSegFPOP_dir(data.dir, "1952.6")

  ## Visualize that model.
  ann.colors <- c(
    noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
  library(ggplot2)
  lab.min <- Mono27ac$labels[1, chromStart]
  lab.max <- Mono27ac$labels[.N, chromEnd]
  changes <- fit$segments[, list(
    constraint=ifelse(diff(mean)==0, "equality", "inequality"),
    chromStart=chromEnd[-1],
    chromEnd=chromEnd[-1])]
  gg <- ggplot()+theme_bw()
  if(require(penaltyLearning)){
    gg <- gg+
      penaltyLearning::geom_tallrect(aes(
        xmin=chromStart, xmax=chromEnd,
        fill=annotation),
        color="grey",
        data=Mono27ac$labels)+
      scale_fill_manual("label", values=ann.colors)
  }
  gg <- gg+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=Mono27ac$coverage)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=fit$segments)+
    geom_segment(aes(
      chromStart, mean,
      xend=chromEnd, yend=mean),
      color="green",
      size=1,
      data=fit$segments)+
    geom_vline(aes(
      xintercept=chromEnd, linetype=constraint),
      color="green",
      data=changes)+
    scale_linetype_manual(values=c(inequality="dotted", equality="solid"))
  gg

  gg+
    coord_cartesian(xlim=c(lab.min, lab.max))

})

