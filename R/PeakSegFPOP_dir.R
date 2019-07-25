PeakSegFPOP_dir <- structure(function # PeakSeg penalized solver with caching
### Main function/interface for the PeakSegDisk package.
### Run the low-level solver, PeakSegFPOP_file,
### on one genomic segmentation problem
### directory, and read the result files into R. Actually, this
### function will first check if the result files are already present
### (and consistent), and if so, it will simply read them into R
### (without running PeakSegFPOP_file) -- this is a caching mechanism
### that can save a lot of time.
### To run the algo on an integer vector, use PeakSegFPOP_vec;
### for a data.frame, use PeakSegFPOP_df.
### To compute the optimal model for a given number of peaks,
### use sequentialSearch_dir.
(problem.dir,
### Path to a directory like sampleID/problems/chrXX-start-end which
### contains a coverage.bedGraph file with the aligned read counts for
### one genomic segmentation problem. Note that the standard
### coverage.bedGraph file name is required; for full flexibility the
### user can run the algo on an arbitrarily named file via
### PeakSegFPOP_file.
 penalty.param
### non-negative numeric penalty parameter (larger values for fewer
### peaks), or character scalar which can be interpreted as such. 0
### means max peaks, Inf means no peaks. 
){
  ##details<<
  ## Finds the optimal change-points using the Poisson loss and the
  ## PeakSeg constraint. For N data points, the functional pruning
  ## algorithm is O(N log N) time and disk space, and O(log N) memory.
  ## It computes the exact
  ## solution to the following optimization problem. Let Z be an
  ## N-vector of count data, typically the coverage, number of aligned
  ## DNA sequence reads in part of the genome
  ## (the fourth column of bedGraph.file, non-negative integers). Let W
  ## be an N-vector of positive weights
  ## (number of bases with the given amount of count/coverage,
  ## chromEnd - chromStart,
  ## third column of bedGraph.file - second column). Let penalty
  ## be a non-negative real number
  ## (larger for fewer peaks, smaller for more peaks).
  ## Find the N-vector M of real numbers
  ## (segment means) and (N-1)-vector C of change-point indicators in
  ## {-1,0,1} which minimize the penalized Poisson Loss,
  ## penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
  ## w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
  ## change is up and the next change is down, etc (sum_{i=1}^t c_i in
  ## {0,1} for all t<N-1), and (2) the last change is down
  ## 0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
  ## variable has an equal segment mean after: c_i=0 implies
  ## m_i=m_{i+1}, (4) every positive-valued change-point variable may
  ## have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
  ## negative-valued change-point variable may have a down change
  ## after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
  ## constraints are active for non-zero change-point variables, the
  ## recovered model is not feasible for the strict inequality
  ## constraints of the PeakSeg problem, and the optimum of the PeakSeg
  ## problem is undefined.
  if(!(
    is.character(problem.dir) && 
    length(problem.dir)==1 &&
    dir.exists(problem.dir))){
    stop(
      "problem.dir=", problem.dir,
      " must be the name of a directory",
      " containing a file named coverage.bedGraph")
  }
  ##alias<< PeakSegDisk
  if(!(
    (is.numeric(penalty.param) || is.character(penalty.param)) &&
    length(penalty.param)==1 &&
    (!is.na(penalty.param))
  )){
    stop("penalty.param must be numeric or character, length 1, not missing") 
  }
  penalty.str <- paste(penalty.param)
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
  L <- list(
    segments=penalty.segs,
    loss=data.table(
      penalty.loss,
      timing[, list(megabytes, seconds)]))
  class(L) <- c("PeakSegFPOP_dir", "list")
  L
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

  data(Mono27ac, package="PeakSegDisk", envir=environment())
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
  (fit <- PeakSegDisk::PeakSegFPOP_dir(data.dir, 1952.6))

  ## Visualize that model.
  ann.colors <- c(
    noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
  library(ggplot2)
  lab.min <- Mono27ac$labels[1, chromStart]
  lab.max <- Mono27ac$labels[.N, chromEnd]

  plist <- coef(fit)
  gg <- ggplot()+
    theme_bw()+
    geom_rect(aes(
      xmin=chromStart/1e3, xmax=chromEnd/1e3,
      ymin=-Inf, ymax=Inf,
      fill=annotation),
      color="grey",
      alpha=0.5,
      data=Mono27ac$labels)+
    scale_fill_manual("label", values=ann.colors)+
    geom_step(aes(
      chromStart/1e3, count),
      color="grey50",
      data=Mono27ac$coverage)+
    geom_segment(aes(
      chromStart/1e3, mean,
      xend=chromEnd/1e3, yend=mean),
      color="green",
      size=1,
      data=plist$segments)+
    geom_vline(aes(
      xintercept=chromEnd/1e3, linetype=constraint),
      color="green",
      data=plist$changes)+
    scale_linetype_manual(
      values=c(
        inequality="dotted",
        equality="solid"))
  gg

  gg+
    coord_cartesian(xlim=c(lab.min, lab.max)/1e3, ylim=c(0, 10))

  ## Default plotting method only shows model.
  (gg <- plot(fit))

  ## Data can be added on top of model.
  gg+
    geom_step(aes(
      chromStart, count),
      color="grey50",
      data=Mono27ac$coverage)
    
})

### Compute changes and peaks to display/plot.
coef.PeakSegFPOP_dir <- function(object, ...){
  chromEnd <- status <- NULL
  ## above to avoid CRAN check NOTE.
  object$changes <- object$segments[, list(
    type="segmentation",
    constraint=ifelse(
      diff(mean)==0, "equality", "inequality"),
    chromEnd=chromEnd[-1])]
  object$peaks <- data.table(
    type="peaks",
    object$segments[status=="peak"])
  object$segments <- data.table(type="segmentation", object$segments)
  object
### model list with additional named elements peaks and changes.
}
  
### Plot a PeakSeg model with attached data.
plot.PeakSegFPOP_dir <- function(x, ...){
  chromStart <- type <- chromEnd <- constraint <- NULL
  ## above to avoid CRAN check NOTE.
  if(!requireNamespace("ggplot2")){
    stop("install ggplot2 for plotting functionality")
  }
  L <- coef(x)
  ggplot2::ggplot()+
    ggplot2::theme_bw()+
    ggplot2::scale_color_manual(values=c(
      data="grey50",
      peaks="deepskyblue",
      segmentation="green"
    ))+
    ggplot2::geom_segment(ggplot2::aes(
      chromStart+0.5, mean,
      color=type,
      xend=chromEnd+0.5, yend=mean),
      size=1,
      data=L$segments)+
    ggplot2::geom_segment(ggplot2::aes(
      chromStart+0.5, Inf,
      color=type,
      xend=chromEnd+0.5, yend=Inf),
      size=2,
      data=L$peaks)+
    ggplot2::geom_point(ggplot2::aes(
      chromStart+0.5, Inf,
      color=type),
      shape=1,
      size=2,
      data=L$peaks)+
    ggplot2::geom_vline(ggplot2::aes(
      xintercept=chromEnd+0.5,
      color=type,
      linetype=constraint),
      data=L$changes)+
    ggplot2::scale_linetype_manual(
      values=c(
        inequality="dotted",
        equality="solid"))+
    ggplot2::xlab("position")
### a ggplot.
}

