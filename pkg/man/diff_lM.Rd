% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errorchecks.R
\name{diff_lM}
\alias{diff_lM}
\title{Calculate Difference Matrices Between lefkoMat Objects of Equal Dimensions}
\usage{
diff_lM(mpm1, mpm2)
}
\arguments{
\item{mpm1}{The first \code{lefkoMat} object.}

\item{mpm2}{The second \code{lefkoMat} object.}
}
\value{
An object of class \code{lefkoDiff}, which is a set of \code{A},
\code{U}, and \code{F} matrices corresponding to the differences between each
set of matrices, followed by the \code{hstages}, \code{ahstages}, and
\code{labels} elements from each input \code{lefkoMat} object. Elements
labelled with a \code{1} at the end refer to \code{mpm1}, while those
labelled \code{2} at the end refer to \code{mpm2}.
}
\description{
Function \code{diff_lM()} takes two \code{lefkoMat} objects with completely
equal dimensions, including both the size and number of matrices, and
gives the matrix differences between each corresponding set.
}
\section{Notes}{

The exact difference is calculated as the respective matrix in \code{mpm1}
minus the corresponding matrix in \code{mpm2}.

This function first checks to see if the number of matrices is the same, and
then whether the matrix dimensions are the same. If the two sets differ in at
least one of these characteristics, then the function will yield a fatal
error.

If the lengths and dimensions of the input \code{lefkoMat} objects are the
same, then this will check if the \code{labels} element is essentially the
same. If not, then the function will yield a warning, but will still operate.
}

\examples{
sizevector <- c(0, 0, 0, 0, 0, 0, 1, 3, 6, 11, 19.5)
stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
  "XLg")
repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1.5, 1.5, 3.5, 5)
comments <- c("Dormant seed", "1st yr protocorm", "2nd yr protocorm",
  "3rd yr protocorm", "Seedling", "Dormant adult",
  "Extra small adult (1 shoot)", "Small adult (2-4 shoots)",
  "Medium adult (5-7 shoots)", "Large adult (8-14 shoots)",
  "Extra large adult (>14 shoots)")
cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
  repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
  propstatus = propvector, immstatus = immvector, indataset = indataset, 
  binhalfwidth = binvec, comments = comments)

cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
  patchidcol = "patch", individcol = "plantid", blocksize = 4,
  sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
  repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
  stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
  NRasRep = TRUE)

seeds_per_pod <- 5000

cypsupp2_raw <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
    "XSm", "SD", "P1"),
  stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep"),
  eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", NA, NA),
  eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
  givenrate = c(0.03, 0.15, 0.1, 0.1, 0.1, 0.05, NA, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, (0.5 * seeds_per_pod),
    (0.5 * seeds_per_pod)),
  type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
  stageframe = cypframe_raw, historical = FALSE)
cypsupp3_raw <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3",
    "SL", "SL", "SL", "D", "D", "SD", "P1"),
  stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
    "rep", "rep"),
  stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "SL", "P3",
    "SL", "mat", "mat"),
  eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "D", NA, NA),
  eststage2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
  eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", NA, NA),
  givenrate = c(0.01, 0.05, 0.10, 0.20, 0.1, 0.1, 0.05, 0.05, 0.05, NA, NA,
    NA, NA),
  multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    (0.5 * seeds_per_pod), (0.5 * seeds_per_pod)),
  type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
  type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  stageframe = cypframe_raw, historical = TRUE)

cypmatrix2rp <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw,
  year = "all", patch = "all", stages = c("stage3", "stage2"),
  size = c("size3added", "size2added"), supplement = cypsupp2_raw, 
  yearcol = "year2", patchcol = "patchid", indivcol = "individ")

cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw,
  year = "all", stages = c("stage3", "stage2"),
  size = c("size3added", "size2added"), supplement = cypsupp2_raw, 
  yearcol = "year2", patchcol = "patchid", indivcol = "individ")

cypmatrix3rp <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
  year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"), 
  size = c("size3added", "size2added", "size1added"), supplement = cypsupp3_raw, 
  yearcol = "year2", patchcol = "patchid", indivcol = "individ")

cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
  year = "all", stages = c("stage3", "stage2", "stage1"), 
  size = c("size3added", "size2added", "size1added"), supplement = cypsupp3_raw, 
  yearcol = "year2", patchcol = "patchid", indivcol = "individ")

cypmatrix2r_3 <- hist_null(cypmatrix2r)
cypmatrix2r_3 <- delete_lM(cypmatrix2r_3, year = 2004)
diff_r <- diff_lM(cypmatrix3r, cypmatrix2r_3)

cypmatrix2rp_3 <- hist_null(cypmatrix2rp)
cypmatrix2rp_3 <- delete_lM(cypmatrix2rp_3, year = 2004)
diff_rp <- diff_lM(cypmatrix3rp, cypmatrix2rp_3)

}
