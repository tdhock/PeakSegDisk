\name{PeakSegFPOP_dir}
\alias{PeakSegFPOP_dir}
\title{PeakSegFPOP dir}
\description{Run PeakSegFPOP_file on one genomic segmentation problem
directory, and read the result files into R. Actually, this
function will first check if the result files are already present
(and consistent), and if so, it will simply read them into R
(without running PeakSegFPOP_file) -- this is a caching mechanism
that can save a lot of time. }
\usage{PeakSegFPOP_dir(problem.dir, penalty.param)}
\arguments{
  \item{problem.dir}{Path to a directory like sampleID/problems/problemID which
contains a coverage.bedGraph file with the aligned read counts for
one genomic segmentation problem.}
  \item{penalty.param}{non-negative numeric penalty parameter (larger values for fewer
peaks), or character scalar which can be interpreted as such. 0
means max peaks, Inf means no peaks. }
}
\details{Finds the optimal change-points using the Poisson loss and the
PeakSeg constraint. For N data points, the functional pruning
algorithm is O(N log N) time and disk space, and O(log N) memory.
It computes the exact
solution to the following optimization problem. Let Z be an
N-vector of count data, typically the coverage, number of aligned
DNA sequence reads in part of the genome
(the fourth column of bedGraph.file, non-negative integers). Let W
be an N-vector of positive weights
(number of bases with the given amount of count/coverage,
chromEnd - chromStart,
third column of bedGraph.file - second column). Let penalty
be a non-negative real number
(larger for fewer peaks, smaller for more peaks).
Find the N-vector M of real numbers
(segment means) and (N-1)-vector C of change-point indicators in
{-1,0,1} which minimize the penalized Poisson Loss,
penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
change is up and the next change is down, etc (sum_{i=1}^t c_i in
{0,1} for all t<N-1), and (2) the last change is down
0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
variable has an equal segment mean after: c_i=0 implies
m_i=m_{i+1}, (4) every positive-valued change-point variable may
have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
negative-valued change-point variable may have a down change
after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
constraints are active for non-zero change-point variables, the
recovered model is not feasible for the strict inequality
constraints of the PeakSeg problem, and the optimum of the PeakSeg
problem is undefined.}
\value{Named list of two data.tables: segments has one row for every
segment in the optimal model, and loss has one row and contains
the following columns. penalty=same as input, segments=number of
segments in optimal model, peaks=number of peaks in optimal model,
bases=number of positions described in bedGraph file,
bedGraph.lines=number of lines in bedGraph file, total.loss=total
Poisson loss=sum_i
m_i-z_i*log(m_i)=mean.pen.cost*bases-penalty*peaks,
mean.pen.cost=mean penalized
cost=(total.loss+penalty*peaks)/bases, equality.constraints=number
of adjacent segment means that have equal values in the optimal
solution, mean.intervals=mean number of intervals/candidate
changepoints stored in optimal cost functions -- useful for
characterizing the computational complexity of the algorithm,
max.intervals=maximum number of intervals, megabytes=disk usage of
*.db file, seconds=timing of PeakSegFPOP_file.}

\author{Toby Dylan Hocking}




\examples{

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
fit <- PeakSegFPOP_dir(data.dir, 1952.6)

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

}
