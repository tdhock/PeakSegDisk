\name{PeakSegFPOP_df}
\alias{PeakSegFPOP_df}
\title{PeakSegFPOP df}
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

## plot method.
plot(fit)

}
