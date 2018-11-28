### Named list of character vectors (column names of bed/bedGraph/tsv
### files).
col.name.list <- list(
  loss=c(
    "penalty", "segments", "peaks", "bases", "bedGraph.lines",
    "mean.pen.cost", "total.loss", "equality.constraints",
    "mean.intervals", "max.intervals"),
  segments=c("chrom","chromStart", "chromEnd", "status", "mean"),
  coverage=c("chrom", "chromStart", "chromEnd", "count"))

### mclapply with error checking.
mclapplyError <- function(...){
  result.list <- parallel::mclapply(...)
  is.error <- sapply(result.list, inherits, "try-error")
  if(any(is.error)){
    print(result.list[is.error])
    stop("errors in mclapply")
  }
  result.list
}

PeakSegFPOP_disk <- structure(function # PeakSegFPOP on disk
### Run the PeakSeg Functional Pruning Optimal Partitioning algorithm,
### using a file on disk (rather than in memory as in
### PeakSegOptimal::PeakSegFPOP) to store the O(N) function piece lists,
### each of size O(log N).
### This is a low-level function that just runs the algo
### and produces the result files (without reading them into R),
### so normal users are recommended to instead use problem.PeakSegFPOP,
### which calls this function then reads the result files into R.
### Finds the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points, the functional pruning
### algorithm is O(N log N) time and disk space, and O(log N) memory.
### It computes the exact
### solution to the following optimization problem. Let Z be an
### N-vector of count data, typically the coverage, number of aligned
### DNA sequence reads in part of the genome
### (the fourth column of bedGraph.file, non-negative integers). Let W
### be an N-vector of positive weights
### (number of bases with the given amount of count/coverage,
### chromEnd - chromStart,
### third column of bedGraph.file - second column). Let penalty
### be a non-negative real number
### (larger for fewer peaks, smaller for more peaks).
### Find the N-vector M of real numbers
### (segment means) and (N-1)-vector C of change-point indicators in
### {-1,0,1} which minimize the penalized Poisson Loss,
### penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
### w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
### change is up and the next change is down, etc (sum_{i=1}^t c_i in
### {0,1} for all t<N-1), and (2) the last change is down
### 0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
### variable has an equal segment mean after: c_i=0 implies
### m_i=m_{i+1}, (4) every positive-valued change-point variable may
### have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
### negative-valued change-point variable may have a down change
### after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
### constraints are active for non-zero change-point variables, the
### recovered model is not feasible for the strict inequality
### constraints of the PeakSeg problem, and the optimum of the PeakSeg
### problem is undefined.
(bedGraph.file,
### character scalar: tab-delimited tabular text file with four
### columns: chrom, chromStart, chromEnd, coverage. The algorithm
### creates a large temporary file in the same directory, so make sure
### that there is disk space available on that device.
  pen.str
### character scalar that can be converted to a numeric scalar via
### as.numeric: non-negative penalty. More penalty means fewer
### peaks. "0" and "Inf" are OK. Character is required rather than
### numeric, so that the user can reliably find the results in the
### output files, which are in the same directory as bedGraph.file,
### and named using the penalty value,
### e.g. coverage.bedGraph_penalty=136500650856.439_loss.tsv
){
  if(!(
    is.character(bedGraph.file) &&
    length(bedGraph.file)==1 &&
    file.exists(bedGraph.file)
  )){
    stop("bedGraph.file must be the name of a data file to segment")
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
    stop("penalty=", penalty, " but it must be a non-negative numeric scalar")
  }
  norm.file <- normalizePath(bedGraph.file, mustWork=TRUE)
  result <- .C(
    "PeakSegFPOP_interface",
    bedGraph.file=as.character(norm.file),
    penalty=pen.str,
    PACKAGE="PeakSegDisk")
  prefix <- paste0(bedGraph.file, "_penalty=", pen.str)
  result$segments <- paste0(prefix, "_segments.bed")
  result$db <- paste0(prefix, ".db")
  result$loss <- paste0(prefix, "_loss.tsv")
  if(file.size(result$loss)==0){
    stop(
      "unable to write to loss output file ",
      result$loss,
      " (disk is probably full)"
    )
  }
  result
### A list of input parameters (bedGraph.file, penalty) and result
### files (segments, db, loss).
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
  names.list <- PeakSegFPOP_disk(tmp, "10.5")
  unlink(names.list$db)
  seg.df <- read.table(names.list$segments)
  names(seg.df) <- col.name.list$segments
  seg.df
  loss.df <- read.table(names.list$loss)
  names(loss.df) <- col.name.list$loss
  loss.df

})

fread.first <- function
### Read the first line of a text file.
(file.name,
### Name of file to read.
  col.name.vec
### Character vector of column names.
){
  dt <- fread(file.name, nrows=1L, col.names=col.name.vec)
  dt
### Data table with one row.
}

fread.last <- function
### Read the last line of a text file.
(file.name,
### Name of file to read.
  col.name.vec
### Character vector of column names.
){
  wc.cmd <- paste("wc -l", file.name)
  wc.output <- system(wc.cmd, intern=TRUE)
  lines.chr <- sub(" .*", "", wc.output)
  lines.int <- as.integer(lines.chr)
  dt <- fread(file.name, skip=lines.int-1L, col.names=col.name.vec)
  dt
### Data table with one row.
}

problem.PeakSegFPOP <- structure(function
### Run PeakSegFPOP_disk on one genomic segmentation problem
### directory, and read the result files into R. Actually, this
### function will first check if the result files are already present
### (and consistent), and if so, it will simply read them into R
### (without running PeakSegFPOP_disk) -- this is a caching mechanism
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
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
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
      PeakSegFPOP_disk(prob.cov.bedGraph, penalty.str)
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
### *.db file, seconds=timing of PeakSegFPOP_disk.
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
  fit <- problem.PeakSegFPOP(data.dir, "1952.6")

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

problem.sequentialSearch <- structure(function
### Compute the most likely peak model with at most the number of
### peaks given by peaks.int. This function repeated calls
### problem.PeakSegFPOP with different penalty values, until either
### (1) it finds the peaks.int model, or (2) it concludes that there
### is no peaks.int model, in which case it returns the next simplest
### model (with fewer peaks than peaks.int).
(problem.dir,
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
  peaks.int,
### int: target number of peaks.
  verbose=0
### Print messages?
){
  stopifnot(
    is.integer(peaks.int) &&
    length(peaks.int)==1 &&
    0 <= peaks.int)
  ## above to avoid "no visible binding for global variable" NOTEs in
  ## CRAN check.
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  model.list <- list()
  next.pen <- c(0, Inf)
  iteration <- 0
  under <- over <- data.table(peaks=NA)
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "mc.cores=", getOption("mc.cores"),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    model.list[next.str] <- mclapplyError(
      next.str, function(penalty.str){
        L <- problem.PeakSegFPOP(problem.dir, penalty.str)
        L$loss$iteration <- iteration
        L$loss$under <- under$peaks
        L$loss$over <- over$peaks
        L
      }
    )
    if(iteration==1){
      under <- model.list[["Inf"]]$loss
      over <- model.list[["0"]]$loss
      max.peaks <- floor((over$bases-1)/2)
      if(max.peaks < peaks.int){
        stop(
          "peaks.int=",
          peaks.int,
          " but max=",
          max.peaks,
          " peaks for N=",
          over$bases,
          " data")
      }
    }else{
      Mnew <- model.list[[next.str]]$loss
      if(Mnew$peaks %in% c(under$peaks, over$peaks)){## not a new model.
        candidate <- under ##pick the simpler one.
        next.pen <- NULL
      }else{#new model.
        if(Mnew$peaks < peaks.int){
          under <- Mnew
        }else{
          over <- Mnew
        }
      }
    }
    if(peaks.int==under$peaks){
      candidate <- under
      next.pen <- NULL
    }
    if(peaks.int==over$peaks){
      candidate <- over
      next.pen <- NULL
    }
    if(!is.null(next.pen)){
      next.pen <- (over$total.loss-under$total.loss)/(under$peaks-over$peaks)
      if(next.pen<0){
        ## sometimes happens for a large number of peaks -- cost is
        ## numerically unstable so we don't get a good penalty to try --
        ## anyways these models are way too big, so just return under.
        candidate <- under
        next.pen <- NULL
      }
    }
  }#while(!is.null(pen))
  out <- model.list[[paste(candidate$penalty)]]
  loss.list <- lapply(model.list, "[[", "loss")
  out$others <- do.call(rbind, loss.list)[order(iteration)]
  out
### Same result list from problem.PeakSegFPOP, with an additional
### component "others" describing the other models that were computed
### before finding the optimal model with peaks.int (or fewer)
### peaks. Additional loss columns are as follows: under=number of
### peaks in smaller model during binary search; over=number of peaks
### in larger model during binary search; iteration=number of times
### PeakSegFPOP has been run.
}, ex=function(){

  ## Create simple 6 point data set discussed in supplementary
  ## materials. GFPOP/GPDPA computes up-down model with 2 peaks, but
  ## neither CDPA (PeakSegDP::cDPA) nor PDPA (jointseg)
  r <- function(chrom, chromStart, chromEnd, coverage){
    data.frame(chrom, chromStart, chromEnd, coverage)
  }
  supp <- rbind(
    r("chr1", 0, 1,  3),
    r("chr1", 1, 2, 9),
    r("chr1", 2, 3, 18),
    r("chr1", 3, 4, 15),
    r("chr1", 4, 5, 20),
    r("chr1", 5, 6, 2)
  )
  data.dir <- file.path(tempfile(), "chr1-0-6")
  dir.create(data.dir, recursive=TRUE)
  write.table(
    supp, file.path(data.dir, "coverage.bedGraph"),
    sep="\t", row.names=FALSE, col.names=FALSE)

  ## Compute optimal up-down model with 2 peaks via sequential search.
  fit <- PeakSegDisk::problem.sequentialSearch(data.dir, 2L)

  library(ggplot2)
  ggplot()+
    theme_bw()+
    geom_point(aes(
      chromEnd, coverage),
      data=supp)+
    geom_segment(aes(
      chromStart+0.5, mean,
      xend=chromEnd+0.5, yend=mean),
      data=fit$segments,
      color="green")

})
