\name{PeakSegFPOP_vec}
\alias{PeakSegFPOP_vec}
\title{PeakSegFPOP vec}
\description{Convert rle to integer vector, write to disk, then run PeakSegFPOP solver.}
\usage{PeakSegFPOP_vec(count.vec, pen.num, base.dir = tempdir())}
\arguments{
  \item{count.vec}{data.frame with columns count, chromStart, chromEnd.}
  \item{pen.num}{Non-negative numeric scalar.}
  \item{base.dir}{
}
}

\value{List of solver results, same as problem.PeakSegFPOP.}

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

fit <- PeakSegFPOP_vec(z.rep.vec, 10.5)
gg.change+
  geom_segment(aes(
    chromStart+0.5, mean, xend=chromEnd+0.5, yend=mean),
    color="green",
    data=fit$segments)

}
