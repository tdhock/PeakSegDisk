library(PeakSegDisk)
X <- rbind(
  c(0, 0),
  c(0, 0),
  c(0, 0))
res <- .C(
  "Eigen_interface",
  X=as.double(X),
  nrow(X),
  ncol(X),
  PACKAGE="PeakSegDisk",
  DUP=FALSE)
matrix(res$X, nrow(X), ncol(X))
