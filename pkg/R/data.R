#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Horizontal
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009.
#' 
#' @docType data
#' 
#' @usage data(cypdata)
#' 
#' @format A data frame with 77 individuals and 29 variables. Each row 
#' corresponds to an unique individual, and each variable from \code{size.04} 
#' on refers to the state of the individual in a particular year.
#' 
#' \describe{
#'   \item{plantid}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{X}{An X coordinate for the plant within the population.}
#'   \item{Y}{A Y coordinate for the plant within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{Inf2.04}{Number of double inflorescences in 2004.}
#'   \item{Inf.04}{Number of inflorescences in 2004.}
#'   \item{Veg.04}{Number of stems without inflorescences in 2004.}
#'   \item{Pod.04}{Number of fruits in 2004.}
#'   \item{Inf2.05}{Number of double inflorescences in 2005.}
#'   \item{Inf.05}{Number of inflorescences in 2005.}
#'   \item{Veg.05}{Number of stems without inflorescences in 2005.}
#'   \item{Pod.05}{Number of fruits in 2005.}
#'   \item{Inf2.06}{Number of double inflorescences in 2006.}
#'   \item{Inf.06}{Number of inflorescences in 2006.}
#'   \item{Veg.06}{Number of stems without inflorescences in 2006.}
#'   \item{Pod.06}{Number of fruits in 2006.}
#'   \item{Inf2.07}{Number of double inflorescences in 2007.}
#'   \item{Inf.07}{Number of inflorescences in 2007.}
#'   \item{Veg.07}{Number of stems without inflorescences in 2007.}
#'   \item{Pod.07}{Number of fruits in 2007.}
#'   \item{Inf2.08}{Number of double inflorescences in 2008.}
#'   \item{Inf.08}{Number of inflorescences in 2008.}
#'   \item{Veg.08}{Number of stems without inflorescences in 2008.}
#'   \item{Pod.08}{Number of fruits in 2008.}
#'   \item{Inf2.09}{Number of double inflorescences in 2009.}
#'   \item{Inf.09}{Number of inflorescences in 2009.}
#'   \item{Veg.09}{Number of stems without inflorescences in 2009.}
#'   \item{Pod.09}{Number of fruits in 2009.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' # Cypripedium example using blocksize
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#'                        
#' cypmatrix2r$A[[intersect(which(cypmatrix2r$labels$patch == "A"), 
#'   which(cypmatrix2r$labels$year2 == 2004))]]
#' 
#' lambda3(cypmatrix2r)
#' 
#' # Cypripedium example using partial repeat patterns with blocksize and part
#' # explicit variable name cast
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = c("Inf.04", "Inf.05", "Inf.06", "Inf.07", "Inf.08", "Inf.09"),
#'   repstrbcol = c("Inf2.04", "Inf2.05", "Inf2.06", "Inf2.07", "Inf2.08", "Inf2.09"), 
#'   fecacol = "Pod.04", stageassign = cypframe_raw, stagesize = "sizeadded",
#'   NAas0 = TRUE, NRasRep = TRUE)
#' 
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#'                        
#' cypmatrix2r$A[[intersect(which(cypmatrix2r$labels$patch == "A"), 
#'   which(cypmatrix2r$labels$year2 == 2004))]]
#' 
#' lambda3(cypmatrix2r)
"cypdata"

#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Vertical
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009. Same dataset as \code{cypdata}, but arranged in an ahistorical
#' vertical format.
#' 
#' @docType data
#' 
#' @usage data(cypvert)
#' 
#' @format A data frame with 77 individuals, 322 rows, and 14 variables. Each
#' row corresponds to a specific two-year transition for a specific individual.
#' Variable codes are similar to those for \code{cypdata}, but use \code{.2} to
#' identify occasion \emph{t} and \code{.3} to identify occasion \emph{t}+1.
#' 
#' \describe{
#'   \item{plantid}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{X}{An X coordinate for the plant within the population.}
#'   \item{Y}{A Y coordinate for the plant within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{year2}{Year in occasion \emph{t}.}
#'   \item{Inf2.2}{Number of double inflorescences in occasion \emph{t}.}
#'   \item{Inf.2}{Number of inflorescences in occasion \emph{t}.}
#'   \item{Veg.2}{Number of stems without inflorescences in occasion \emph{t}.}
#'   \item{Pod.2}{Number of fruits in occasion \emph{t}.}
#'   \item{Inf2.3}{Number of double inflorescences in occasion \emph{t}+1.}
#'   \item{Inf.3}{Number of inflorescences in occasion \emph{t}+1.}
#'   \item{Veg.3}{Number of stems without inflorescences in occasion \emph{t}+1.}
#'   \item{Pod.3}{Number of fruits in occasion \emph{t}+1.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' data(cypvert)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypframe_raw
#' 
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", 
#'   individcol = "plantid", year2col = "year2", sizea2col = "Inf2.2", 
#'   sizea3col = "Inf2.3", sizeb2col = "Inf.2", sizeb3col = "Inf.3", 
#'   sizec2col = "Veg.2", sizec3col = "Veg.3", repstra2col = "Inf2.2", 
#'   repstra3col = "Inf2.3", repstrb2col = "Inf.2", repstrb3col = "Inf.3", 
#'   feca2col = "Pod.2", feca3col = "Pod.3", repstrrel = 2, 
#'   stageassign = cypframe_raw, stagesize = "sizeadded", censorcol = "censor",
#'   censor = FALSE, NAas0 = TRUE, NRasRep = TRUE, reduce = TRUE)
#'   
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v2, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#'   
#' cypmatrix2r$A[[intersect(which(cypmatrix2r$labels$patch == "A"), 
#'   which(cypmatrix2r$labels$year2 == 2004))]]
#' 
#' lambda3(cypmatrix2r)
"cypvert"

#' Demographic Dataset of \emph{Lathyrus vernus} Population
#' 
#' A dataset containing the states and fates of \emph{Lathyrus vernus} (spring
#' vetch), family Fabaceae, from a population in Sweden monitored annually
#' from 1988 to 1991 in six study plots.
#' 
#' @docType data
#' 
#' @usage data(lathyrus)
#' 
#' @format A data frame with 1119 individuals and 34 variables. Each row
#' corresponds to a unique individual, and each variable from \code{Volume88}
#' on refers to the state of the individual in a given year.
#' 
#' \describe{
#'   \item{SUBPLOT}{A variable refering to patch within the population.}
#'   \item{GENET}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{Volume88}{Aboveground volume in cubic mm in 1988.}
#'   \item{lnVol88}{Natural logarithm of \code{Volume88}.}
#'   \item{FCODE88}{Equals 1 if flowering and 0 if not flowering in 1988.}
#'   \item{Flow88}{Number of flowers in 1988.}
#'   \item{Intactseed88}{Number of intact mature seeds produced in 1988.
#'   Not always an integer, as in some cases seed number was estimated via 
#'   linear modeling.}
#'   \item{Dead1988}{Marked as 1 if known to be dead in 1988.}
#'   \item{Dormant1988}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1988.}
#'   \item{Missing1988}{Marked as 1 if not found in 1988.}
#'   \item{Seedling1988}{Marked as 1, 2, or 3 if observed as a seedling in year
#'   \emph{t}. Numbers refer to certainty of assignment: 1 = certain that plant
#'   is a seedling in 1988, 2 = likely that plant is a seedling in 1988,
#'   3 = probable that plant is a seedling in 1988.}
#'   \item{Volume89}{Aboveground volume in cubic mm in 1989.}
#'   \item{lnVol89}{Natural logarithm of \code{Volume89}.}
#'   \item{FCODE89}{Equals 1 if flowering and 0 if not flowering in 1989.}
#'   \item{Flow89}{Number of flowers in 1989.}
#'   \item{Intactseed89}{NZumber of intact mature seeds produced in 1989.
#'   Not always an integer, as in some cases seed number was estimated via
#'   linear modeling.}
#'   \item{Dead1989}{Marked as 1 if known to be dead in 1989.}
#'   \item{Dormant1989}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1989.}
#'   \item{Missing1989}{Marked as 1 if not found in 1989.}
#'   \item{Seedling1989}{Marked as 1, 2, or 3 if observed as a seedling in
#'   year \emph{t}. Numbers refer to certainty of assignment: 1 = certain 
#'   that plant is a seedling in 1989, 2 = likely that plant is a seedling 
#'   in 1989, 3 = probable that plant is a seedling in 1989.}
#'   \item{Volume90}{Aboveground volume in mm<sup>3</sup> in 1990.}
#'   \item{lnVol90}{Natural logarithm of \code{Volume90}.}
#'   \item{FCODE90}{Equals 1 if flowering and 0 if not flowering in 1990.}
#'   \item{Flow90}{Number of flowers in 1990.}
#'   \item{Intactseed90}{NZumber of intact mature seeds produced in 1990.
#'   Not always an integer, as in some cases seed number was estimated via 
#'   linear modeling.}
#'   \item{Dead1990}{Marked as 1 if known to be dead in 1990.}
#'   \item{Dormant1990}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1990.}
#'   \item{Missing1990}{Marked as 1 if not found in 1990.}
#'   \item{Seedling1990}{Marked as 1, 2, or 3 if observed as a seedling in
#'   year \emph{t}. Numbers refer to certainty of assignment: 1 = certain 
#'   that plant is a seedling in 1990, 2 = likely that plant is a seedling
#'   in 1990, 3 = probable that plant is a seedling in 1990.}
#'   \item{Volume91}{Aboveground volume in mm<sup>3</sup> in 1991.}
#'   \item{lnVol91}{Natural logarithm of \code{Volume91}.}
#'   \item{FCODE91}{Equals 1 if flowering and 0 if not flowering in 1991.}
#'   \item{Flow91}{Number of flowers in 1991.}
#'   \item{Intactseed91}{NZumber of intact mature seeds produced in 1991.
#'   Not always an integer, as in some cases seed number was estimated via
#'   linear modeling.}
#'   \item{Dead1991}{Marked as 1 if known to be dead in 1991.}
#'   \item{Dormant1991}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1991.}
#'   \item{Missing1991}{Marked as 1 if not found in 1991.}
#'   \item{Seedling1991}{Marked as 1, 2, or 3 if observed as a seedling 
#'   in year \emph{t}. Numbers refer to certainty of assignment: 
#'   1 = certain that plant is a seedling in 1991, 2 = likely that plant 
#'   is a seedling in 1991, 3 = probable that plant is a seedling in 
#'   1991.}
#' }
#' 
#' @source Ehrlen, J. 2000. The dynamics of plant populations: does the 
#' history of individuals matter? \emph{Ecology} 81(6):1675-1684.
#' 
#' @examples
#' # Lathyrus example using blocksize - when repeated patterns exist in variable
#' # order
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' ehrlen3mean$A[[1]]
#' 
#' lambda3(ehrlen3mean)
#' 
#' # Lathyrus example without blocksize - when no repeated patterns exist in
#' # variable order and all variables names are specified
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET",
#'   juvcol = c("Seedling1988", "Seedling1989", "Seedling1990", "Seedling1991"),
#'   sizeacol = c("Volume88", "Volume89", "Volume90", "Volume91"),
#'   repstracol = c("FCODE88", "FCODE89", "FCODE90", "FCODE91"),
#'   fecacol = c("Intactseed88", "Intactseed89", "Intactseed90", "Intactseed91"),
#'   deadacol = c("Dead1988", "Dead1989", "Dead1990", "Dead1991"),
#'   nonobsacol = c("Dormant1988", "Dormant1989", "Dormant1990", "Dormant1991"),
#'   censorcol = c("Missing1988", "Missing1989", "Missing1990", "Missing1991"), 
#'   stageassign = lathframe, stagesize = "sizea",
#'   censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' ehrlen3mean$A[[1]]
#' 
#' lambda3(ehrlen3mean)
"lathyrus"

#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Vertical
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009. Same dataset as \code{cypdata}, but arranged in an ahistorical
#' vertical format.
#' 
#' @docType data
#' 
#' @usage data(cypvert)
#' 
#' @format A data frame with 77 individuals, 322 rows, and 14 variables. Each
#' row corresponds to a specific two-year transition for a specific individual.
#' Variable codes are similar to those for \code{cypdata}, but use \code{.2} to
#' identify occasion \emph{t} and \code{.3} to identify occasion \emph{t}+1.
#' 
#' \describe{
#'   \item{plantid}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{X}{An X coordinate for the plant within the population.}
#'   \item{Y}{A Y coordinate for the plant within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{year2}{Year in occasion \emph{t}.}
#'   \item{Inf2.2}{Number of double inflorescences in occasion \emph{t}.}
#'   \item{Inf.2}{Number of inflorescences in occasion \emph{t}.}
#'   \item{Veg.2}{Number of stems without inflorescences in occasion \emph{t}.}
#'   \item{Pod.2}{Number of fruits in occasion \emph{t}.}
#'   \item{Inf2.3}{Number of double inflorescences in occasion \emph{t}+1.}
#'   \item{Inf.3}{Number of inflorescences in occasion \emph{t}+1.}
#'   \item{Veg.3}{Number of stems without inflorescences in occasion \emph{t}+1.}
#'   \item{Pod.3}{Number of fruits in occasion \emph{t}+1.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' data(cypvert)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypframe_raw
#' 
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", 
#'   individcol = "plantid", year2col = "year2", sizea2col = "Inf2.2", 
#'   sizea3col = "Inf2.3", sizeb2col = "Inf.2", sizeb3col = "Inf.3", 
#'   sizec2col = "Veg.2", sizec3col = "Veg.3", repstra2col = "Inf2.2", 
#'   repstra3col = "Inf2.3", repstrb2col = "Inf.2", repstrb3col = "Inf.3", 
#'   feca2col = "Pod.2", feca3col = "Pod.3", repstrrel = 2, 
#'   stageassign = cypframe_raw, stagesize = "sizeadded", censorcol = "censor",
#'   censor = FALSE, NAas0 = TRUE, NRasRep = TRUE, reduce = TRUE)
#'   
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v2, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#'   
#' cypmatrix2r$A[[intersect(which(cypmatrix2r$labels$patch == "A"), 
#'   which(cypmatrix2r$labels$year2 == 2004))]]
#' 
#' lambda3(cypmatrix2r)
"cypvert"

#' Matrix Set of \emph{Anthyllis vulneraria} Populations in Belgium
#' 
#' A \code{lefkoMat} object containing projection matrices developed from
#' demographic data gathered on nine \emph{Anthyllis vulneraria} populations
#' from 2003 to 2006 in southwestern Belgium. These matrices were originally
#' published in Davison et al. (2010; Journal of Eccology 98(2): 255-267).
#' 
#' @docType data
#' 
#' @usage data(anthyllis)
#' 
#' @format A \code{lefkoMat} object holding 27 matrices. The structure of the
#' object is as below:
#' 
#' \describe{
#'   \item{A}{The 27 A matrices.}
#'   \item{U}{The 27 survival-transition matrices used to develop the A
#'   matrices.}
#'   \item{F}{The 27 fecundity matrices used to develop the A matrices.}
#'   \item{hstages}{Not used, so set to \code{NA}.}
#'   \item{agestages}{Not used, so set to \code{NA}.}
#'   \item{ahstages}{The edited stageframe describing the life history of the
#'   study organism as interpreted in the original demographic study.}
#'   \item{labels}{The order of the matrices, where each population is treated
#'   as a separate patch and each matrix corresponds to a different combination
#'   of population and year in time \emph{t}.}
#'   \item{matrixqc}{A vector of integers used in the quality control section of
#'   \code{lefkoMat} summary statements.}
#' }
#' 
#' @source Davison, R. et al. 2010. Demographic effects of extreme weather
#' events on a short-lived calcareous grassland species: stochastic life table
#' response experiments. \emph{Journal of Ecology} 98(2):255-267.
#' 
#' @examples
#' data(anthyllis)
#' 
#' summary(anthyllis)
#' 
#' lambda3(anthyllis)
"anthyllis"

