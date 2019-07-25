\name{PeakSegFPOP_df}
\alias{PeakSegFPOP_df}
\title{PeakSeg penalized solver for data.frame}
\description{Write data frame to disk then run PeakSegFPOP_dir solver.}
\usage{PeakSegFPOP_df(count.df, pen.num, base.dir = tempdir())}
\arguments{
  \item{count.df}{data.frame with columns count, chromStart, chromEnd.}
  \item{pen.num}{Non-negative numeric scalar.}
  \item{base.dir}{base.dir/chrXX-start-end/coverage.bedGraph will be written, where
chrXX is the chrom column, start is the first chromStart position,
and end is the last chromEnd position.}
}

\value{List of solver results, same as PeakSegFPOP_dir.}

\author{Toby Dylan Hocking}




\examples{

## Simulate a sequence of Poisson count data.
sim.seg <- function(seg.mean, size.mean=15){
  seg.size <- rpois(1, size.mean)
  rpois(seg.size, seg.mean)
}
set.seed(1)
seg.mean.vec <- c(1.5, 3.5, 0.5, 4.5, 2.5)
z.list <- lapply(seg.mean.vec, sim.seg)
z.rep.vec <- unlist(z.list)

## Plot the simulated data sequence.
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

## Plot the true changes.
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

## Plot the run-length encoding of the same data.
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

## Fit a peak model and plot the segment means.
fit <- PeakSegDisk::PeakSegFPOP_df(coverage.df, 10.5)
gg.rle+
  geom_segment(aes(
    chromStart+0.5, mean, xend=chromEnd+0.5, yend=mean),
    color="green",
    data=fit$segments)

## Default plot method shows data as geom_step.
(gg <- plot(fit))

## Plot data as points to verify the step representation.
gg+
  geom_point(aes(
    position, count),
    color="grey",
    shape=1,
    data=count.df)

## The time complexity of the algorithm depends on the number of
## intervals (candidate changepoints stored). Here we compute the
## mean number of intervals for real Mono27ac data, and synthetic
## count data which are always increasing.
data(Mono27ac, package="PeakSegDisk", envir=environment())
library(data.table)
loss.list <- list()
type.colors <- c(
  "Theoretical\nworst case"="grey",
  "Increasing"="red",
  "Mono27ac"="black")
dfac <- function(x){
  factor(x, names(type.colors))
}
N.data.vec <- 10^seq(1, 3)
for(penalty in c(0, 1e2, 1e4, 1e6)){
  for(N.data in N.data.vec){
    some.cov <- Mono27ac$coverage[1:N.data]
    some.inc <- data.table(some.cov)
    some.inc[, count := 1:.N]
    data.list <- list(
      Mono27ac=some.cov,
      Increasing=some.inc)
    for(data.type in names(data.list)){
      df <- data.list[[data.type]]
      fit <- PeakSegDisk::PeakSegFPOP_df(df, penalty)
      loss.list[[paste(penalty, N.data, data.type)]] <- data.table(
        N.data,
        data.type=dfac(data.type),
        fit$loss)
    }
  }
}
(loss <- do.call(rbind, loss.list))

## Theoretically the most intervals that could be stored is O(i) for
## data point i\in{1, ..., N}. Therefore the largest total number of
## intervals is sum(1:N), which can also be computed by
## N*(N+1)/2. The largest mean is mean(1:N), which can be computed
## via sum(1:N)/N = (N+1)/2 = N/2 + 0.5.
worst.dt <- data.table(
  N.data=N.data.vec,
  mean.intervals=(N.data.vec+1)/2,
  data.type=dfac("Theoretical\nworst case"))

## The plot below shows that the algorithm achieves the theoretical
## worst case time complexity for the synthetic increasing data,
## when the penalty is large. But the number of intervals is always
## much smaller for real Mono27ac ChIP-seq data.
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ penalty, labeller=label_both)+
  geom_line(aes(
    N.data, mean.intervals, color=data.type),
    data=worst.dt)+
  scale_color_manual(
    values=type.colors,
    breaks=names(type.colors))+
  geom_line(aes(
    bedGraph.lines, mean.intervals, color=data.type),
    data=loss)+
  scale_x_log10()+
  scale_y_log10()

}
