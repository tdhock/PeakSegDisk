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

  ## Longer example on a bigger data set.
  if(interactive()){

    ## Load data set with one row for every genomic region with a
    ## unique aligned read, and compute mean read size in bases.
    data(ChIPreads, package="PeakSegDisk", envir=environment())
    experiments <- ChIPreads[, .(
      mean.bases=mean(chromEnd-chromStart),
      median.bases=median(chromEnd-chromStart),
      chromStart=min(chromStart)
    ), by=list(experiment)]
    
    ## Compute data set with two representations of these aligned
    ## reads: count each read at each aligned base in the read, or
    ## just the end/last base of the read.
    end.counts <- ChIPreads[, list(
      count=.N #ignores dup reads, sum(count) would not.
    ), by=list(experiment, chrom, chromEnd)]
    aligned.dt <- rbind(
      ChIPreads[, .(
        bases.counted="each", experiment, chrom,
        chromStart, chromEnd,
        count=1)], #ignore duplicate reads.
      end.counts[, .(
        bases.counted="end", experiment, chrom,
        chromStart=chromEnd-1L, chromEnd,
        count)])

    ## Compute count profile for each base in these genomic regions.
    library(data.table)
    seq.dt <- aligned.dt[, {
      event.dt <- rbind(
        data.table(count, pos=chromStart+1L),
        data.table(count=-count, pos=chromEnd+1L))
      edge.vec <- event.dt[, {
        as.integer(seq(min(pos), max(pos), l=100))
      }]
      event.bins <- rbind(
        event.dt,
        data.table(count=0L, pos=edge.vec))
      total.dt <- event.bins[, .(
        count=sum(count)
      ), by=list(pos)][order(pos)]
      total.dt[, cum := cumsum(count)]
      total.dt[, bin.i := cumsum(pos %in% edge.vec)]
      ## it is somewhat confusing because total.dt pos is the first base
      ## with cum, and cum goes all the way up to but not including the
      ## pos of the next row.
      total.dt[, data.table(
        chromStart=pos[-.N]-1L,
        chromEnd=pos[-1]-1L,
        count=cum[-.N],
        bin.i=bin.i[-.N])]
    }, by=list(bases.counted, experiment, chrom)]
    gg.data <- ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_grid(
        bases.counted ~ experiment,
        scales="free",
        labeller=label_both)+
      geom_step(aes(
        chromStart/1e3, count, color=data.type),
        data=data.table(seq.dt, data.type="exact"))+
      scale_color_manual(values=c(
        exact="black",
        bins="red",
        model="deepskyblue"
      ))+
      scale_x_continuous("Position on hg19 chrom (kb = kilo bases)")
    print(gg.data)

    ## Compute mean profile in bins.
    bin.dt <- seq.dt[, {
      bases <- chromEnd - chromStart
      data.table(
        binStart=min(chromStart),
        binEnd=max(chromEnd),
        mean.count=sum(count*bases)/sum(bases),
        bases=sum(bases)
      )}, by=list(bases.counted, experiment, bin.i)]
    gg.bins <- gg.data+
      geom_step(aes(
        binStart/1e3, mean.count, color=data.type),
        alpha=0.75,
        size=1,
        data=data.table(bin.dt, data.type="bins"))+
      scale_y_log10("Aligned DNA sequence reads (log scale)")
    print(gg.bins)

    ## Compute optimal segmentation model with 2 peaks.
    segs.dt <- seq.dt[, {
      data.dir <- file.path(tempdir(), bases.counted, experiment)
      dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
      coverage.bedGraph <- file.path(data.dir, "coverage.bedGraph")
      fwrite(
        .SD[, .(chrom, chromStart, chromEnd, count)],
        coverage.bedGraph,
        sep="\t",
        quote=FALSE,
        col.names=FALSE)
      fit <- PeakSegDisk::sequentialSearch_dir(data.dir, 2L, verbose=1)
      data.table(fit$segments, data.type="model")
    }, by=list(bases.counted, experiment)]
    changes.dt <- segs.dt[, {
      .SD[-1]
    }, by=list(bases.counted, experiment, data.type)]
    gg.model <- gg.bins+
      geom_segment(aes(
        chromStart/1e3, mean,
        xend=chromEnd/1e3, yend=mean,
        color=data.type),
        data=segs.dt)+
      geom_vline(aes(
        xintercept=chromEnd/1e3,
        color=data.type),
        data=changes.dt)
    print(gg.model)

    ## Compute difference between peak positions of two models.
    peaks.dt <- segs.dt[status=="peak"]
    peaks.dt[, peak.i := rep(1:2, l=.N)]
    peak.pos.tall <- melt(
      peaks.dt,
      measure.vars=c("chromStart", "chromEnd"))
    peak.pos.wide <- dcast(
      peak.pos.tall,
      experiment + variable + peak.i ~ bases.counted)
    peak.pos.wide[, diff.bases := abs(each-end)]

    read.size.panel <- "each"
    bases.max.dt <- seq.dt[, .(max.count=max(count)), by=list(bases.counted)]
    read.size.y <- bases.max.dt[
      read.size.panel, max.count, on=list(bases.counted)]
    diff.panel <- "end"
    diff.y <- bases.max.dt[
      diff.panel, max.count, on=list(bases.counted)]
    diff.y <- Inf
    diff.vjust <- 1.1
    gg.model+
      geom_text(aes(
        chromStart/1e3, read.size.y, label=sprintf(
          "Median read size:\n%.0f bases",
          median.bases)),
        hjust=0,
        vjust=1,
        data=data.table(experiments, bases.counted=read.size.panel))+
      geom_text(aes(
        end/1e3, diff.y,
        label=diff.bases,
        color=data.type),
        data=data.table(
          bases.counted=diff.panel,
          data.type="model",
          peak.pos.wide),
        vjust=diff.vjust,
        hjust=0)+
      geom_text(aes(
        chromStart/1e3, diff.y,
        label="Peak position\ndifference in bases:",
        color=data.type,
      ),
      hjust=0,
      vjust=diff.vjust,
      data=data.table(
        data.type="model",
        bases.counted=diff.panel,
        experiments["H3K36me3", on=list(experiment)]))
        

  }

})
