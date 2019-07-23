sequentialSearch_dir <- structure(function # Compute PeakSeg model with given number of peaks
### Compute the most likely peak model with at most the number of
### peaks given by peaks.int. This function repeated calls
### PeakSegFPOP_dir with different penalty values, until either
### (1) it finds the peaks.int model, or (2) it concludes that there
### is no peaks.int model, in which case it returns the next simplest
### model (with fewer peaks than peaks.int).
### The first pair of penalty values (0, Inf) is run in parallel
### via the user-specified future plan,
### if the future.apply package is available.
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
  LAPPLY <- if(requireNamespace("future.apply")){
    future.apply::future_lapply
  }else{
    lapply
  }
  while(length(next.pen)){
    if(verbose)cat(
      "Next =", paste(next.pen, collapse=", "),
      "\n")
    next.str <- paste(next.pen)
    iteration <- iteration+1
    model.list[next.str] <- LAPPLY(
      next.str, function(penalty.str){
        L <- PeakSegFPOP_dir(problem.dir, penalty.str)
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
### Same result list from PeakSegFPOP_dir, with an additional
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
  library(PeakSegDisk)
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
  fit <- PeakSegDisk::sequentialSearch_dir(data.dir, 2L)

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

  data(ChIPreads, envir=environment())
  end.counts <- ChIPreads[, list(
    count=.N #ignores dup reads, sum(count) would not.
  ), by=list(experiment, chrom, chromEnd)]

  aligned.dt <- rbind(
    ChIPreads[, .(
      data.type="each", experiment, chrom, chromStart=chromStart, chromEnd,
      count=1)], #ignore duplicate reads.
    end.counts[, .(
      data.type="end", experiment, chrom,
      chromStart=chromEnd-1L, chromEnd, count)])
  seq.dt <- aligned.dt[, {
    event.dt <- rbind(
      data.table(count, pos=chromStart+1L),
      data.table(count=-count, pos=chromEnd+1L))
    total.dt <- event.dt[, .(
      count=sum(count)
      ), by=list(pos)][order(pos)]
    total.dt[, cum := cumsum(count)]
    ## it is somewhat confusing because total.dt pos is the first base
    ## with cum, and cum goes all the way up to but not including the
    ## pos of the next row.
    total.dt[, data.table(
      chromStart=pos[-.N]-1L,
      chromEnd=pos[-1]-1L,
      count=cum[-.N])]
  }, by=list(data.type, experiment, chrom)]

  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(data.type ~ experiment, scales="free")+
    geom_step(aes(
      chromEnd, count),
      data=seq.dt)

})
