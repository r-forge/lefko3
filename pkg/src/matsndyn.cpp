// [[Rcpp::depends(RcppArmadillo)]]
#define BOOST_DISABLE_ASSERTS

#include <RcppArmadilloExtensions/sample.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using namespace arma;

// Pop Char

//' Append NumericVector to the End of Another NumericVector
//' 
//' This function appends one NumericVector fully to another.
//' 
//' @param A Any NumericVector.
//' @param B Any other NumericVector.
//' 
//' @return Returns a new NumericVector with elements of vector A followed by
//' elements of vector B.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.concat_dbl)]]
NumericVector concat_dbl(NumericVector x, NumericVector y) {
  
  std::vector<double> xconv = as<std::vector<double> >(x);
  std::vector<double> yconv = as<std::vector<double> >(y);
  std::vector<double> z(xconv.size() + yconv.size());
  
  std::copy(xconv.begin(), xconv.end(), z.begin());
  std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
  
  NumericVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Append IntegerVector to the End of Another IntegerVector
//' 
//' Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @param A Any IntegerVector.
//' @param B Any other IntegerVector.
//' 
//' @return Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.concat_int)]]
IntegerVector concat_int(IntegerVector x, IntegerVector y) {
  
  std::vector<long long int> xconv = as<std::vector<long long int> >(x);
  std::vector<long long int> yconv = as<std::vector<long long int> >(y);
  std::vector<long long int> z(x.size() + y.size());
  
  std::copy(xconv.begin(), xconv.end(), z.begin());
  std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
  
  IntegerVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Append StringVector to the End of Another StringVector
//' 
//' Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @param A Any StringVector.
//' @param B Any other StringVector.
//' 
//' @return Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.concat_str)]]
StringVector concat_str(StringVector x, StringVector y) {
  
  std::vector<std::string> xconv = as<std::vector<std::string> >(x);
  std::vector<std::string> yconv = as<std::vector<std::string> >(y);
  std::vector<std::string> z(x.size() + y.size());
  
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  
  StringVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Create Stageframe for Population Matrix Projection Analysis
//' 
//' Function \code{sf_create()} returns a data frame describing each ahistorical
//' life history stage in the life history model. This data frame can be used as 
//' input into MPM creation functions including \code{\link{flefko3}()}, 
//' \code{\link{flefko2}()}, \code{\link{aflefko2}()}, \code{\link{rlefko3}()},
//' and \code{\link{rlefko2}()}, in which it determines how each stage is
//' treated during matrix estimation.
//' 
//' @param sizes A numeric vector of the typical or representative size of each
//' life history stage. If making function-based MPMs, then this should be a
//' vector composed of the midpoints of each size bin. If denoting the boundary
//' of an automated size classification group, then should denote the absolute
//' minimum size of that group, or the absolute size of that group (see Notes).
//' @param stagenames A vector of stage names, in the same order as elements in
//' sizes. Can also be set to \code{ipm} for automated size classification (see
//' Notes section).
//' @param sizesb An optional numeric vector for a second size metric for each
//' life history stage. Only to be used if stages are defined by at least two
//' size metrics in all cases. Same issues apply as in \code{sizes}.
//' @param sizesc An optional numeric vector for a third size metric for each
//' life history stage. Only to be used if stages are defined by at least three
//' size metrics in all cases. Same issues apply as in \code{sizes}.
//' @param repstatus A vector denoting the binomial reproductive status of each
//' life history stage. Defaults to 1.
//' @param obsstatus A vector denoting the binomial observation status of each
//' life history stage. Defaults to 1, but may be changed for unobservable 
//' stages.
//' @param propstatus A vector denoting whether each life history stage is a 
//' propagule. Such stages are generally only used in fecundity estimation. 
//' Defaults to 0.
//' @param matstatus A vector denoting whether each stage is mature. Must be
//' composed of binomial values if given. Defaults to 1 for all stages defined 
//' in \code{sizes}.
//' @param immstatus A vector denoting whether each stage is immature. Must be
//' composed of binomial values if given. Defaults to the complement of vector
//' \code{matstatus}.
//' @param minage An optional vector denoting the minimum age at which a stage
//' can occur. Only used in age x stage matrix development. Defaults to NA.
//' @param maxage An optional vector denoting the maximum age at which a stage
//' should occur. Only used in age x stage matrix development. Defaults to NA.
//' @param indataset A vector designating which stages are found within the 
//' dataset. While \code{\link{rlefko2}()} and \code{\link{rlefko3}()} can use
//' all stages in the input dataset, \code{\link{flefko3}()} and
//' \code{\link{flefko2}()} can only handle size-classified stages with
//' non-overlapping combinations of size and status variables. Stages that do
//' not actually exist within the dataset should be marked as 0 in this vector.
//' @param binhalfwidth A numeric vector giving the half-width of size bins.
//' Required to classify individuals appropriately within size classes.
//' Defaults to 0.5 for all sizes.
//' @param binhalfwidthb A numeric vector giving the half-width of size bins
//' used for the optional second size metric. Required to classify individuals
//' appropriately with two or three size classes. Defaults to 0.5 for all sizes.
//' @param binhalfwidthc A numeric vector giving the half-width of size bins
//' used for the optional third size metric. Required to classify individuals
//' appropriately with three size classes. Defaults to 0.5 for all sizes.
//' @param group An integer vector providing information on each respective
//' stage's size classification group. If used, then function-based MPM creation
//' functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
//' \code{\link{aflefko2}()} will estimate transitions only within these groups
//' and for allowed cross-group transitions noted within the supplement table.
//' Defaults to 0.
//' @param comments An optional vector of text entries holding useful text
//' descriptions of all stages.
//' @param roundsize This parameter sets the precision of size classification,
//' and equals the number of digits used in rounding sizes. Defaults to 5.
//' @param roundsizeb This parameter sets the precision of size classification
//' in the optional second size metric, and equals the number of digits used in
//' rounding sizes. Defaults to 5.
//' @param roundsizec This parameter sets the precision of size classification
//' in the optional third size metric, and equals the number of digits used in
//' rounding sizes. Defaults to 5.
//' @param ipmbins An integer giving the number of size bins to create using the
//' primary size classification variable. This number is in addition to any
//' stages that are not size classified. Defaults to 100, and numbers greater
//' than this yield a warning about the loss of statistical power and increasing
//' chance of matrix over-parameterization resulting from increasing numbers of
//' stages.
//' @param ipmbinsb An optional integer giving the number of size bins to create
//' using the secondary size classification variable. This number is in addition
//' to any stages that are not size classified, as well as in addition to any
//' automated size classification using the primary and tertiary size variables.
//' Defaults to NA, and must be set to a positive integer for automated size
//' classification to progress.
//' @param ipmbinsc An optional integer giving the number of size bins to create
//' using the tertiary size classification variable. This number is in addition
//' to any stages that are not size classified, as well as in addition to any
//' automated size classification using the primary and secondary size
//' variables. Defaults to NA, and must be set to a positive integer for
//' automated size classification to progress.
//' 
//' @return A data frame of class \code{stageframe}, which includes information
//' on the stage name, size, reproductive status, observation status, propagule 
//' status, immaturity status, maturity status, presence within the core dataset, 
//' stage group classification, raw bin half-width, and the minimum, 
//' center, and maximum of each size bin, as well as its width. If minimum and
//' maximum ages were specified, then these are also included. Also includes an 
//' empty string variable that can be used to describe stages meaningfully. This
//' object can be used as the \code{stageframe} input for \code{\link{flefko3}()} 
//' \code{\link{flefko2}()}, \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.
//' 
//' Variables in this data frame include the following:
//' \item{stage}{The unique names of the stages to be analyzed.}
//' \item{size}{The typical or representative size at which each stage occurs.}
//' \item{size_b}{Size at which each stage occurs in terms of a second size
//' variable, if one exists.}
//' \item{size_c}{Size at which each stage occurs in terms of a third size
//' variable, if one exists.}
//' \item{min_age}{The minimum age at which the stage may occur.}
//' \item{max_age}{The maximum age at which the stage may occur.}
//' \item{repstatus}{A binomial variable showing whether each stage is
//' reproductive.}
//' \item{obsstatus}{A binomial variable showing whether each stage is
//' observable.}
//' \item{propstatus}{A binomial variable showing whether each stage is a
//' propagule.}
//' \item{immstatus}{A binomial variable showing whether each stage can occur as
//' immature.}
//' \item{matstatus}{A binomial variable showing whether each stage occurs in
//' maturity.}
//' \item{indataset}{A binomial variable describing whether each stage occurs in
//' the input dataset.}
//' \item{binhalfwidth_raw}{The half-width of the size bin, as input.}
//' \item{sizebin_min}{The minimum size at which the stage may occur.}
//' \item{sizebin_max}{The maximum size at which the stage may occur.}
//' \item{sizebin_center}{The midpoint of the size bin at which the stage may
//' occur.}
//' \item{sizebin_width}{The width of the size bin corresponding to the stage.}
//' \item{binhalfwidthb_raw}{The half-width of the size bin of a second size
//' variable, as input.}
//' \item{sizebinb_min}{The minimum size at which the stage may occur.}
//' \item{sizebinb_max}{The maximum size at which the stage may occur.}
//' \item{sizebinb_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a second size variable.}
//' \item{sizebinb_width}{The width of the size bin corresponding to the stage,
//' in terms of a second size variable.}
//' \item{binhalfwidthc_raw}{The half-width of the size bin of a third size
//' variable, as input.}
//' \item{sizebinc_min}{The minimum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_max}{The maximum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a third size variable.}
//' \item{sizebinc_width}{The width of the size bin corresponding to the stage,
//' in terms of a third size variable.}
//' \item{group}{An integer denoting the size classification group that the
//' stage falls within.}
//' \item{comments}{A text field for stage descriptions.}
//' 
//' @section Notes:
//' If an IPM or function-based matrix with automated size classification is
//' desired, then two stages that occur within the dataset and represent the
//' lower and upper size limits of the IPM must be marked with \code{ipm} in
//' the stagenames vector. These stages should have all characteristics other
//' than size equal, and the size input for whichever size will be classified
//' automatically must include the minimum in one stage and the maximum in the
//' other. The actual characteristics of the first stage encountered in the
//' inputs will be used as the template for the creation of these sizes. Note
//' that \code{ipm} refers to size classification with the primary size
//' variable. To automate size classification with the secondary size variable,
//' use \code{ipmb}, and to automate size classification with the tertiary size
//' variable, use \code{ipmc}. To nest automated size classifications, use 
//' \code{ipmab} for the primary and secondary size variables, \code{ipmac} for
//' the primary and tertiary size variables, \code{ipmbc} for the secondary and
//' tertiary size variables, and \code{ipmabc} for all three size variables.
//' The primary size variable can also be set with \code{ipma}.
//' 
//' If two or more groups of stages, each with its own characteristics, are to
//' be developed for an IPM or function-based MPM, then an even number of stages
//' with two stages marking the minimum and maximum size of each group should be
//' marked with the same code as given above, with all other characteristics
//' equal within each group.
//' 
//' Stage classification groups set with the \code{group} variable create zones
//' within function-based matrices in which survival transitions are estimated.
//' These groups should not be set if transitions are possible between all
//' stages regardless of group. To denote specific transitions as estimable
//' between stage groups, use the \code{\link{supplemental}()} function.
//' 
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
//'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
//'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
//'   yearcol = "year2", indivcol = "individ")
//' 
//' ehrlen3mean <- lmean(ehrlen3)
//' ehrlen3mean$A[[1]]
//' 
//' # Cypripedium example
//' data(cypdata)
//' 
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' 
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//'                        
//' cyp2mean <- lmean(cypmatrix2r)
//' cyp2mean
//' 
//' @export sf_create
// [[Rcpp::export]]
Rcpp::List sf_create (NumericVector sizes, Nullable<StringVector> stagenames = R_NilValue,
  Nullable<NumericVector> sizesb = R_NilValue, Nullable<NumericVector> sizesc = R_NilValue,
  Nullable<IntegerVector> repstatus = R_NilValue, Nullable<IntegerVector> obsstatus = R_NilValue,
  Nullable<IntegerVector> propstatus = R_NilValue, Nullable<IntegerVector> matstatus = R_NilValue,
  Nullable<IntegerVector> immstatus = R_NilValue, Nullable<NumericVector> minage = R_NilValue,
  Nullable<NumericVector> maxage = R_NilValue, Nullable<IntegerVector> indataset = R_NilValue,
  Nullable<NumericVector> binhalfwidth = R_NilValue, Nullable<NumericVector> binhalfwidthb = R_NilValue,
  Nullable<NumericVector> binhalfwidthc = R_NilValue, Nullable<IntegerVector> group = R_NilValue,
  Nullable<StringVector> comments = R_NilValue, int roundsize = 5, int roundsizeb = 5,
  int roundsizec = 5, int ipmbins = 100, int ipmbinsb = NA_INTEGER, int ipmbinsc = NA_INTEGER) {
  
  Rcpp::List output_longlist(29);
  Rcpp::CharacterVector varnames {"stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center", "sizebin_width",
    "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
    "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
    "group", "comments"};
  
  int matsize = sizes.size(); // Core vector size
  int used_sizes {1};
  arma::uvec ipm_calls_a (matsize, fill::zeros);
  arma::uvec ipm_calls_b (matsize, fill::zeros);
  arma::uvec ipm_calls_c (matsize, fill::zeros);
  arma::uvec ipm_calls_ab (matsize, fill::zeros);
  arma::uvec ipm_calls_bc (matsize, fill::zeros);
  arma::uvec ipm_calls_ac (matsize, fill::zeros);
  arma::uvec ipm_calls_abc (matsize, fill::zeros);
  
  StringVector stagenames_true (matsize, NA_STRING);
  NumericVector sizesb_true (matsize, NA_REAL);
  NumericVector sizesc_true (matsize, NA_REAL);
  IntegerVector repstatus_true (matsize, 1);
  IntegerVector obsstatus_true (matsize, 1);
  IntegerVector propstatus_true (matsize, 0);
  IntegerVector matstatus_true (matsize, 1);
  IntegerVector immstatus_true (matsize, 0);
  IntegerVector indataset_true (matsize, 1);
  NumericVector minage_true (matsize, NA_REAL);
  NumericVector maxage_true (matsize, NA_REAL);
  NumericVector binhalfwidth_true (matsize, 0.5);
  NumericVector binhalfwidthb_true (matsize, NA_REAL);
  NumericVector binhalfwidthc_true (matsize, NA_REAL);
  NumericVector sizebin_min (matsize, NA_REAL);
  NumericVector sizebin_max (matsize, NA_REAL);
  NumericVector sizebin_center (matsize, NA_REAL);
  NumericVector sizebin_width (matsize, NA_REAL);
  NumericVector sizebinb_min (matsize, NA_REAL);
  NumericVector sizebinb_max (matsize, NA_REAL);
  NumericVector sizebinb_center (matsize, NA_REAL);
  NumericVector sizebinb_width (matsize, NA_REAL);
  NumericVector sizebinc_min (matsize, NA_REAL);
  NumericVector sizebinc_max (matsize, NA_REAL);
  NumericVector sizebinc_center (matsize, NA_REAL);
  NumericVector sizebinc_width (matsize, NA_REAL);
  IntegerVector group_true (matsize, 0);
  StringVector comments_true (matsize, "No description");
  
  bool no_age = true;
  
  if (stagenames.isNotNull()) {
    Rcpp::StringVector stagenames_thru(stagenames);
    
    if (stagenames_thru.length() == matsize) {
      stagenames_true = stagenames_thru;
      
      StringVector st_t(stagenames_thru.size());
      
      std::transform(stagenames_thru.begin(), stagenames_thru.end(), st_t.begin(), 
        make_string_transformer(tolower));
      
      for (int i = 0; i < stagenames_thru.size(); i++) {
        String check_elem = trimws(st_t(i));
        
        if (check_elem == "ipm" || check_elem == "ipma" || check_elem == "ipm_a" ||
          check_elem == "ipm1" || check_elem == "ipm_1") {
          ipm_calls_a(i) = 1;
        } else if (check_elem == "ipmb" || check_elem == "ipm_b" || check_elem == "ipm2" ||
          check_elem == "ipm_2") {
          ipm_calls_b(i) = 1;
        } else if (check_elem == "ipmc" || check_elem == "ipm_c" || check_elem == "ipm3" ||
          check_elem == "ipm_3") {
          ipm_calls_c(i) = 1;
        } else if (check_elem == "ipmab" || check_elem == "ipm_ab" || check_elem == "ipm12" ||
          check_elem == "ipm_12") {
          ipm_calls_ab(i) = 1;
        } else if (check_elem == "ipmac" || check_elem == "ipm_ac" || check_elem == "ipm13" ||
          check_elem == "ipm_13") {
          ipm_calls_ac(i) = 1;
        } else if (check_elem == "ipmbc" || check_elem == "ipm_bc" || check_elem == "ipm23" ||
          check_elem == "ipm_23") {
          ipm_calls_bc(i) = 1;
        } else if (check_elem == "ipmabc" || check_elem == "ipm_abc" || check_elem == "ipm123" ||
          check_elem == "ipm_123") {
          ipm_calls_abc(i) = 1;
        }
      }
    } else if (stagenames_thru.length() == 1) {
      
      if (!stagenames_thru.is_na(0)) {
        for (int i = 0; i < matsize; i++) {
          Rcpp::String part1(stagenames_thru(0));
          part1 += (static_cast<char>(i + 1));
          
          stagenames_true(i) = part1;
        }
      } else {
        for (int i = 0; i < matsize; i++) {
          Rcpp::String part1("Stage");
          part1 += (static_cast<char>(i + 1));
          
          stagenames_true(i) = part1;
        }
      }
    } else {
      throw Rcpp::exception("Vector stagenames should be the same length as vector sizes.", false);
    }
  } else {
    for (int i = 0; i < matsize; i++) {
      Rcpp::String part1("Stage");
      part1 += (static_cast<char>(i + 1));
      
      stagenames_true(i) = part1;
    }
  }
  
  if (sizesb.isNotNull()) {
    Rcpp::NumericVector sizesb_thru(sizesb);
    used_sizes++;
    
    if (sizesb_thru.length() == matsize) {
      sizesb_true = sizesb_thru;
    } else if (sizesb_thru.length() == 1) {
      NumericVector try_size (matsize, sizesb_thru(0));
      sizesb_true = try_size;
    } else {
      throw Rcpp::exception("Vector sizesb should be the same length as vector sizes.", false);
    }
  }
  if (sizesc.isNotNull()) {
    Rcpp::NumericVector sizesc_thru(sizesc);
    used_sizes++;
    
    if (sizesc_thru.length() == matsize) {
      sizesc_true = sizesc_thru;
    } else if (sizesc_thru.length() == 1) {
      NumericVector try_size (matsize, sizesc_thru(0));
      sizesc_true = try_size;
    } else {
      throw Rcpp::exception("Vector sizesc should be the same length as vector sizes.", false);
    }
    
    if (used_sizes != 3) {
      throw Rcpp::exception("Vector sizesc should only be set if vector sizesb is also set.", false);
    }
  }
  if (minage.isNotNull()) {
    Rcpp::NumericVector minage_thru(minage);
    
    if (minage_thru.length() == matsize) {
      minage_true = minage_thru;
    } else if (minage_thru.length() == 1) {
      NumericVector try_mna (matsize, minage_thru(0));
      minage_true = try_mna;
    } else {
      throw Rcpp::exception("Vector minage should be the same length as vector sizes.", false);
    }
    
    no_age = false;
  }
  if (maxage.isNotNull()) {
    Rcpp::NumericVector maxage_thru(maxage);
    
    if (maxage_thru.length() == matsize) {
      maxage_true = maxage_thru;
    } else if (maxage_thru.length() == 1) {
      NumericVector try_mxa (matsize, maxage_thru(0));
      maxage_true = try_mxa;
    } else {
      throw Rcpp::exception("Vector maxage should be the same length as vector sizes.", false);
    }
    //if (no_age) { // We don't really need this any more.
    //  throw Rcpp::exception("Vector minage is required if vector maxage is provided.", false);
    //}
  } // else if (!no_age) { // We do not really need this bit
    //throw Rcpp::exception("Vector maxage is required if vector minage is provided.", false);
  //}
  
  if (repstatus.isNotNull()) {
    Rcpp::IntegerVector repstatus_thru(repstatus);
    
    if (repstatus_thru.length() == matsize) {
      repstatus_true = repstatus_thru;
    } else if (repstatus_thru.length() == 1) {
      IntegerVector try_rep (matsize, repstatus_thru(0));
      repstatus_true = try_rep;
    } else {
      throw Rcpp::exception("Vector repstatus should be the same length as vector sizes.", false);
    }
    
    if (max(repstatus_true) > 1 || min(repstatus_true) < 0) {
      throw Rcpp::exception("Vector repstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (obsstatus.isNotNull()) {
    Rcpp::IntegerVector obsstatus_thru(obsstatus);
    
    if (obsstatus_thru.length() == matsize) {
      obsstatus_true = obsstatus_thru;
    } else if (obsstatus_thru.length() == 1) {
      IntegerVector try_obs (matsize, obsstatus_thru(0));
      obsstatus_true = try_obs;
    } else {
      throw Rcpp::exception("Vector obsstatus should be the same length as vector sizes.", false);
    }
    
    if (max(obsstatus_true) > 1 || min(obsstatus_true) < 0) {
      throw Rcpp::exception("Vector obsstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (propstatus.isNotNull()) {
    Rcpp::IntegerVector propstatus_thru(propstatus);
    
    if (propstatus_thru.length() == matsize) {
      propstatus_true = propstatus_thru;
    } else if (propstatus_thru.length() == 1) {
      IntegerVector try_prop (matsize, propstatus_thru(0));
      propstatus_true = try_prop;
    } else {
      throw Rcpp::exception("Vector propstatus should be the same length as vector sizes.", false);
    }
    
    if (max(propstatus_true) > 1 || min(propstatus_true) < 0) {
      throw Rcpp::exception("Vector propstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (matstatus.isNotNull()) {
    Rcpp::IntegerVector matstatus_thru(matstatus);
    
    if (matstatus_thru.length() == matsize) {
      matstatus_true = matstatus_thru;
    } else if (matstatus_thru.length() == 1) {
      IntegerVector try_mat (matsize, matstatus_thru(0));
      matstatus_true = try_mat;
    } else {
      throw Rcpp::exception("Vector matstatus should be the same length as vector sizes.", false);
    }
    
    if (max(matstatus_true) > 1 || min(matstatus_true) < 0) {
      throw Rcpp::exception("Vector matstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (immstatus.isNotNull()) {
    Rcpp::IntegerVector immstatus_thru(immstatus);
    
    if (immstatus_thru.length() == matsize) {
      immstatus_true = immstatus_thru;
    } else if (immstatus_thru.length() == 1) {
      IntegerVector try_imm (matsize, immstatus_thru(0));
      immstatus_true = try_imm;
    } else {
      throw Rcpp::exception("Vector immstatus should be the same length as vector sizes.", false);
    }
    
    if (max(immstatus_true) > 1 || min(immstatus_true) < 0) {
      throw Rcpp::exception("Vector immstatus should be composed only of 0s and 1s.", false);
    }
  } else {
    IntegerVector trialones(matsize, 1);
    immstatus_true = trialones - matstatus_true;
  }
  if (indataset.isNotNull()) {
    Rcpp::IntegerVector indataset_thru(indataset);
    
    if (indataset_thru.length() == matsize) {
      indataset_true = indataset_thru;
    } else if (indataset_thru.length() == 1) {
      IntegerVector try_ind (matsize, indataset_thru(0));
      indataset_true = try_ind;
    } else {
      throw Rcpp::exception("Vector indataset should be the same length as vector sizes.", false);
    }
    
    if (max(indataset_true) > 1 || min(indataset_true) < 0) {
      throw Rcpp::exception("Vector indataset should be composed only of 0s and 1s.", false);
    }
  }
  
  if (binhalfwidth.isNotNull()) {
    Rcpp::NumericVector binhalfwidth_thru(binhalfwidth);
    
    if (binhalfwidth_thru.length() == matsize) {
      binhalfwidth_true = binhalfwidth_thru;
    } else if (binhalfwidth_thru.length() == 1) {
      NumericVector try_hwa (matsize, binhalfwidth_thru(0));
      binhalfwidth_true = try_hwa;
    } else {
      throw Rcpp::exception("Vector binhalfwidth should be the same length as vector sizes.", false);
    }
  }
  if (binhalfwidthb.isNotNull()) {
    if (used_sizes == 1) {
      throw Rcpp::exception("Vector binhalfwidthb should only be used if multiple size variables are being used for classification.", false);
    }
    
    Rcpp::NumericVector binhalfwidthb_thru(binhalfwidthb);
    
    if (binhalfwidthb_thru.length() == matsize) {
      binhalfwidthb_true = binhalfwidthb_thru;
    } else if (binhalfwidthb_thru.length() == 1) {
      NumericVector try_hwb (matsize, binhalfwidthb_thru(0));
      binhalfwidthb_true = try_hwb;
    } else {
      throw Rcpp::exception("Vector binhalfwidthb should be the same length as vector sizes.", false);
    }
  } else if (used_sizes > 1) {
    for (int i = 0; i < binhalfwidthb_true.length(); i++) {
      binhalfwidthb_true(i) = 0.5;
    }
  }
  if (binhalfwidthc.isNotNull()) {
    if (used_sizes < 3) {
      throw Rcpp::exception("Vector binhalfwidthc should only be used if three size variables are being used for classification.", false);
    }
    
    Rcpp::NumericVector binhalfwidthc_thru(binhalfwidthc);
    
    if (binhalfwidthc_thru.length() == matsize) {
      binhalfwidthc_true = binhalfwidthc_thru;
    } else if (binhalfwidthc_thru.length() == 1) {
      NumericVector try_hwc (matsize, binhalfwidthc_thru(0));
      binhalfwidthc_true = try_hwc;
    } else {
      throw Rcpp::exception("Vector binhalfwidthc should be the same length as vector sizes.", false);
    }
  } else if (used_sizes > 2) {
    for (int i = 0; i < binhalfwidthc_true.length(); i++) {
      binhalfwidthc_true(i) = 0.5;
    }
  }
  
  if (group.isNotNull()) {
    Rcpp::IntegerVector group_thru(group);
    
    if (group_thru.length() == matsize) {
      group_true = group_thru;
    } else if (group_thru.length() == 1) {
      IntegerVector try_grp (matsize, group_thru(0));
      group_true = try_grp;
    } else {
      throw Rcpp::exception("Vector group should be the same length as vector sizes.", false);
    }
    
    if (min(group_true) < 0) {
      throw Rcpp::exception("Please use only positive integers for group designations.", false);
    }
  }
  
  if (comments.isNotNull()) {
    Rcpp::StringVector comments_thru(comments);
    
    if (comments_thru.length() == matsize) {
      comments_true = comments_thru;
    } else if (comments_thru.length() == 1) {
      StringVector try_com (matsize, comments_thru(0));
      comments_true = try_com;
    } else {
      throw Rcpp::exception("Comments vector should be the same length as vector sizes.", false);
    }
  }
  
  
  if (ipmbins < 2) {
    throw Rcpp::exception("Please enter a valid integer greater than 1 for ipmbins option.", false);
  } else if (ipmbins > 100) {
    Rf_warningcall(R_NilValue, "High ipmbin numbers may lead to dramatic decreases in statistical power and overparameterized matrices.");
  }
  
  // Now the automated size classification processing
  int ipm_a = sum(ipm_calls_a);
  int ipm_b = sum(ipm_calls_b);
  int ipm_c = sum(ipm_calls_c);
  int ipm_ab = sum(ipm_calls_ab);
  int ipm_bc = sum(ipm_calls_bc);
  int ipm_ac = sum(ipm_calls_ac);
  int ipm_abc = sum(ipm_calls_abc);
  
  // This vector will act as an unique marker for most characteristics
  arma::ivec pair_check (matsize, fill::zeros); 
  arma::uvec entries_to_delete (1);
  int no_entries_to_delete {0};
  
  for (int i = 0; i < matsize; i++) {
    pair_check(i) = repstatus_true(i) * 1000000 + obsstatus_true(i) * 100000 + 
      propstatus_true(i) * 10000 + immstatus_true(i) * 1000 + matstatus_true(i) * 100 +
      indataset_true(i) * 10 + group_true(i);
  }

  // Automated size classification with sizea
  if (ipm_a > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_a % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_a); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize = sizes(check_elems(i));
          minsize = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize = sizes(check_elems(i));
          maxsize = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbins, 0.0);
        StringVector newstagenames (ipmbins, "");
        NumericVector newsizesb (ipmbins, 0.0);
        NumericVector newsizesc (ipmbins, 0.0);
        NumericVector newminage (ipmbins, 0.0);
        NumericVector newmaxage (ipmbins, 0.0);
        IntegerVector newrepstatus (ipmbins, 0);
        IntegerVector newobsstatus (ipmbins, 0);
        IntegerVector newpropstatus (ipmbins, 0);
        IntegerVector newmatstatus (ipmbins, 0);
        IntegerVector newimmstatus (ipmbins, 0);
        IntegerVector newindataset (ipmbins, 0);
        NumericVector newbinhalfwidth (ipmbins, 0.0);
        NumericVector newbinhalfwidthb (ipmbins, 0.0);
        NumericVector newbinhalfwidthc (ipmbins, 0.0);
        IntegerVector newgroup (ipmbins, 0);
        StringVector newcomments (ipmbins, "");
        
        for (int j = 0; j < ipmbins; j++) {
          newsizes(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizesb(j) = sizesb_true(check_elems(i));
          newsizesc(j) = sizesc_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizes(j));
          newstagenames(j) = "sza_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // IPM classification with sizeb
  if (ipm_b > 0) {
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the secondary size variable has not been set.", false);
    }
    
    if (ipm_b % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_b); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }

      if (go_ahead == 1) {
        
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize = sizesb_true(check_elems(i));
          minsize = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize = sizesb_true(check_elems(i));
          maxsize = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;

        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbinsb, 0.0);
        StringVector newstagenames (ipmbinsb, "");
        NumericVector newsizesb (ipmbinsb, 0.0);
        NumericVector newsizesc (ipmbinsb, 0.0);
        NumericVector newminage (ipmbinsb, 0.0);
        NumericVector newmaxage (ipmbinsb, 0.0);
        IntegerVector newrepstatus (ipmbinsb, 0);
        IntegerVector newobsstatus (ipmbinsb, 0);
        IntegerVector newpropstatus (ipmbinsb, 0);
        IntegerVector newmatstatus (ipmbinsb, 0);
        IntegerVector newimmstatus (ipmbinsb, 0);
        IntegerVector newindataset (ipmbinsb, 0);
        NumericVector newbinhalfwidth (ipmbinsb, 0.0);
        NumericVector newbinhalfwidthb (ipmbinsb, 0.0);
        NumericVector newbinhalfwidthc (ipmbinsb, 0.0);
        IntegerVector newgroup (ipmbinsb, 0);
        StringVector newcomments (ipmbinsb, "");
        
        for (int j = 0; j < ipmbinsb; j++) {
          newsizesb(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizes(j) = sizes(check_elems(i));
          newsizesc(j) = sizesc_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizesb(j));
          newstagenames(j) = "szb_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // IPM classification with sizec
  if (ipm_c > 0) {
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the tertiary size variable has not been set.", false);
    }
    
    if (ipm_c % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_c); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize = sizesc_true(check_elems(i));
          minsize = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize = sizesc_true(check_elems(i));
          maxsize = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbinsc, 0.0);
        StringVector newstagenames (ipmbinsc, "");
        NumericVector newsizesb (ipmbinsc, 0.0);
        NumericVector newsizesc (ipmbinsc, 0.0);
        NumericVector newminage (ipmbinsc, 0.0);
        NumericVector newmaxage (ipmbinsc, 0.0);
        IntegerVector newrepstatus (ipmbinsc, 0);
        IntegerVector newobsstatus (ipmbinsc, 0);
        IntegerVector newpropstatus (ipmbinsc, 0);
        IntegerVector newmatstatus (ipmbinsc, 0);
        IntegerVector newimmstatus (ipmbinsc, 0);
        IntegerVector newindataset (ipmbinsc, 0);
        NumericVector newbinhalfwidth (ipmbinsc, 0.0);
        NumericVector newbinhalfwidthb (ipmbinsc, 0.0);
        NumericVector newbinhalfwidthc (ipmbinsc, 0.0);
        IntegerVector newgroup (ipmbinsc, 0);
        StringVector newcomments (ipmbinsc, "");
        
        for (int j = 0; j < ipmbinsc; j++) {
          newsizesc(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizesb(j) = sizesb_true(check_elems(i));
          newsizes(j) = sizes(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizesc(j));
          newstagenames(j) = "szc_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea & sizeb together
  if (ipm_ab > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_ab % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_ab); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_b {0};
        double maxsize_b {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsb), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsb), "");
        NumericVector newsizesb ((ipmbins * ipmbinsb), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsb), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsb), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsb), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsb), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsb), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsb), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsb), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsb), 0);
        StringVector newcomments ((ipmbins * ipmbinsb), "");
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsb; k++) {
            newsizes((j * ipmbinsb) + k) = minsize_a + (j*standard_increment_a) + standard_midpoint_a;
            newsizesb((j * ipmbinsb) + k) = minsize_b + (k*standard_increment_b) + standard_midpoint_b;;
            newsizesc((j * ipmbinsb) + k) = sizesc_true(check_elems(i));
            newminage((j * ipmbinsb) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsb) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsb) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsb) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsb) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsb) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsb) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsb) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsb) + k) = standard_midpoint_a;
            newbinhalfwidthb((j * ipmbinsb) + k) = standard_midpoint_b;
            newbinhalfwidthc((j * ipmbinsb) + k) = binhalfwidthc_true(check_elems(i));
            newgroup((j * ipmbinsb) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsb) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizes((j * ipmbinsb) + k));
            std::string sizenums2 = std::to_string(newsizesb((j * ipmbinsb) + k));
            newstagenames((j * ipmbinsb) + k) = "sza_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsb) + k) += "_szb_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsb) + k) += "_";
            newstagenames((j * ipmbinsb) + k) += std::to_string(newgroup((j * ipmbinsb) + k));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea & sizec together
  if (ipm_ac > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_ac % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_ac); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsc), "");
        NumericVector newsizesb ((ipmbins * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsc), 0);
        StringVector newcomments ((ipmbins * ipmbinsc), "");
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsc; k++) {
            newsizes((j * ipmbinsc) + k) = minsize_a + (j*standard_increment_a) + standard_midpoint_a;
            newsizesc((j * ipmbinsc) + k) = minsize_c + (k*standard_increment_c) + standard_midpoint_c;;
            newsizesb((j * ipmbinsc) + k) = sizesb_true(check_elems(i));
            newminage((j * ipmbinsc) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsc) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsc) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsc) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsc) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsc) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsc) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsc) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsc) + k) = standard_midpoint_a;
            newbinhalfwidthb((j * ipmbinsc) + k) = binhalfwidthb_true(check_elems(i));
            newbinhalfwidthc((j * ipmbinsc) + k) = standard_midpoint_c;
            newgroup((j * ipmbinsc) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsc) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizes((j * ipmbinsc) + k));
            std::string sizenums2 = std::to_string(newsizesc((j * ipmbinsc) + k));
            newstagenames((j * ipmbinsc) + k) = "sza_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_szc_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_";
            newstagenames((j * ipmbinsc) + k) += std::to_string(newgroup((j * ipmbinsc)));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizeb & sizec together
  if (ipm_bc > 0) {
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_bc % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_bc); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_b {0};
        double maxsize_b {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbinsb * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbinsb * ipmbinsc), "");
        NumericVector newsizesb ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbinsb * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbinsb * ipmbinsc), 0);
        StringVector newcomments ((ipmbinsb * ipmbinsc), "");
        
        for (int j = 0; j < ipmbinsb; j++) {
          for (int k = 0; k < ipmbinsc; k++) {
            newsizesb((j * ipmbinsc) + k) = minsize_b + (j*standard_increment_b) + standard_midpoint_b;
            newsizesc((j * ipmbinsc) + k) = minsize_c + (k*standard_increment_c) + standard_midpoint_c;;
            newsizes((j * ipmbinsc) + k) = sizes(check_elems(i));
            newminage((j * ipmbinsc) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsc) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsc) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsc) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsc) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsc) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsc) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsc) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsc) + k) = binhalfwidth_true(check_elems(i));
            newbinhalfwidthb((j * ipmbinsc) + k) = standard_midpoint_b;
            newbinhalfwidthc((j * ipmbinsc) + k) = standard_midpoint_c;
            newgroup((j * ipmbinsc) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsc) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizesb((j * ipmbinsc) + k));
            std::string sizenums2 = std::to_string(newsizesc((j * ipmbinsc) + k));
            newstagenames((j * ipmbinsc) + k) = "szb_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_szc_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_";
            newstagenames((j * ipmbinsc) + k) += std::to_string(newgroup((j * ipmbinsc) + k));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea, sizeb, & sizec together
  if (ipm_abc > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_abc % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_abc); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_b {0};
        double maxsize_b {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsb * ipmbinsc), "");
        NumericVector newsizesb ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsb * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsb * ipmbinsc), 0);
        StringVector newcomments ((ipmbins * ipmbinsb * ipmbinsc), "");
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsb; k++) {
            for (int l = 0; l < ipmbinsc; l++) {
              newsizes((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_a + 
                (j*standard_increment_a) + standard_midpoint_a;
              newsizesb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_b + 
                (k*standard_increment_b) + standard_midpoint_b;
              newsizesc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_c + 
                (l*standard_increment_c) + standard_midpoint_c;
              newminage((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minage_true(check_elems(i));
              newmaxage((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = maxage_true(check_elems(i));
              newrepstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = repstatus_true(check_elems(i));
              newobsstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = obsstatus_true(check_elems(i));
              newpropstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = propstatus_true(check_elems(i));
              newmatstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = matstatus_true(check_elems(i));
              newimmstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = immstatus_true(check_elems(i));
              newindataset((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = indataset_true(check_elems(i));
              newbinhalfwidth((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = standard_midpoint_a;
              newbinhalfwidthb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = standard_midpoint_b;
              newbinhalfwidthc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = standard_midpoint_c;
              newgroup((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = group_true(check_elems(i));
              newcomments((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = comments_true(check_elems(i));
              
              std::string sizenums1 = std::to_string(newsizes((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              std::string sizenums2 = std::to_string(newsizesb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              std::string sizenums3 = std::to_string(newsizesc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = "sza_" + sizenums1.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_szb_" + sizenums2.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_szc_" + sizenums3.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_";
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += 
                std::to_string(newgroup((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
            }
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Delete unneeded elements - run for loop backwards to delete elements properly
  arma::uvec unique_delete = unique(entries_to_delete);
  int no_unique_delete = unique_delete.n_elem;
  
  if (no_entries_to_delete > 0) {
    for (int i = 0; i < no_unique_delete; i++) {
      sizes.erase(unique_delete(no_unique_delete - (i + 1)));
      sizesb_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizesc_true.erase(unique_delete(no_unique_delete - (i + 1)));
      stagenames_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      minage_true.erase(unique_delete(no_unique_delete - (i + 1)));
      maxage_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      repstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      obsstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      propstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      matstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      immstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      indataset_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      binhalfwidth_true.erase(unique_delete(no_unique_delete - (i + 1)));
      binhalfwidthb_true.erase(unique_delete(no_unique_delete - (i + 1)));
      binhalfwidthc_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      group_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      comments_true.erase(unique_delete(no_unique_delete - (i + 1)));
    }
  }
  
  // Post-IPM calculations
  matsize = sizes.size(); // Redefined length
  int elems_to_grow = matsize - sizebin_min.length();
  NumericVector zeros_to_append (elems_to_grow, NA_REAL);
  
  sizebin_min = concat_dbl(sizebin_min, zeros_to_append);
  sizebin_max = concat_dbl(sizebin_max, zeros_to_append);
  sizebin_center = concat_dbl(sizebin_center, zeros_to_append);
  sizebin_width = concat_dbl(sizebin_width, zeros_to_append);
  
  sizebinb_min = concat_dbl(sizebinb_min, zeros_to_append);
  sizebinb_max = concat_dbl(sizebinb_max, zeros_to_append);
  sizebinb_center = concat_dbl(sizebinb_center, zeros_to_append);
  sizebinb_width = concat_dbl(sizebinb_width, zeros_to_append);
  
  sizebinc_min = concat_dbl(sizebinc_min, zeros_to_append);
  sizebinc_max = concat_dbl(sizebinc_max, zeros_to_append);
  sizebinc_center = concat_dbl(sizebinc_center, zeros_to_append);
  sizebinc_width = concat_dbl(sizebinc_width, zeros_to_append);
  
  for (int i = 0; i < matsize; i++) {
    sizebin_min(i) = sizes(i) - binhalfwidth_true(i);
    sizebin_max(i) = sizes(i) + binhalfwidth_true(i);
    sizebin_center(i) = sizebin_min(i) + ((sizebin_max(i) - sizebin_min(i))/ 2);
    sizebin_width(i) = sizebin_max(i) - sizebin_min(i);
    
    if (used_sizes > 1) {
      sizebinb_min(i) = sizesb_true(i) - binhalfwidthb_true(i);
      sizebinb_max(i) = sizesb_true(i) + binhalfwidthb_true(i);
      sizebinb_center(i) = sizebinb_min(i) + ((sizebinb_max(i) - sizebinb_min(i))/ 2);
      sizebinb_width(i) = sizebinb_max(i) - sizebinb_min(i);
    }
    if (used_sizes > 2) {
      sizebinc_min(i) = sizesc_true(i) - binhalfwidthc_true(i);
      sizebinc_max(i) = sizesc_true(i) + binhalfwidthc_true(i);
      sizebinc_center(i) = sizebinc_min(i) + ((sizebinc_max(i) - sizebinc_min(i))/ 2);
      sizebinc_width(i) = sizebinc_max(i) - sizebinc_min(i);
    }
  }
  sizebin_min = round(sizebin_min, roundsize);
  sizebin_max = round(sizebin_max, roundsize);
  sizebin_center = round(sizebin_center, roundsize);
  sizebin_width = round(sizebin_width, roundsize);
  
  if (used_sizes > 1) {
    sizebinb_min = round(sizebinb_min, roundsizeb);
    sizebinb_max = round(sizebinb_max, roundsizeb);
    sizebinb_center = round(sizebinb_center, roundsizeb);
    sizebinb_width = round(sizebinb_width, roundsizeb);
  }
  if (used_sizes > 2) {
    sizebinc_min = round(sizebinc_min, roundsizec);
    sizebinc_max = round(sizebinc_max, roundsizec);
    sizebinc_center = round(sizebinc_center, roundsizec);
    sizebinc_width = round(sizebinc_width, roundsizec);
  }
  
  output_longlist(0) = stagenames_true;
  output_longlist(1) = sizes;
  output_longlist(2) = sizesb_true;
  output_longlist(3) = sizesc_true;
  
  output_longlist(4) = minage_true;
  output_longlist(5) = maxage_true;
  output_longlist(6) = repstatus_true;
  output_longlist(7) = obsstatus_true;
  output_longlist(8) = propstatus_true;
  output_longlist(9) = immstatus_true;
  output_longlist(10) = matstatus_true;
  
  output_longlist(11) = indataset_true;
  
  output_longlist(12) = binhalfwidth_true;
  output_longlist(13) = sizebin_min;
  output_longlist(14) = sizebin_max;
  output_longlist(15) = sizebin_center;
  output_longlist(16) = sizebin_width;
  
  output_longlist(17) = binhalfwidthb_true;
  output_longlist(18) = sizebinb_min;
  output_longlist(19) = sizebinb_max;
  output_longlist(20) = sizebinb_center;
  output_longlist(21) = sizebinb_width;
  
  output_longlist(22) = binhalfwidthc_true;
  output_longlist(23) = sizebinc_min;
  output_longlist(24) = sizebinc_max;
  output_longlist(25) = sizebinc_center;
  output_longlist(26) = sizebinc_width;
  
  output_longlist(27) = group_true;
  output_longlist(28) = comments_true;
  
  output_longlist.attr("names") = varnames;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, matsize);
  StringVector needed_classes {"data.frame", "stageframe"};
  output_longlist.attr("class") = needed_classes; // data.frame
  
  return output_longlist;
}

//' Standardize Stageframe For MPM Analysis
//' 
//' Function \code{.sf_reassess()} takes a stageframe as input, and uses
//' information supplied there and through the supplement, reproduction and
//' overwrite tables to rearrange this into a format usable by the matrix
//' creation functions, \code{\link{flefko3}()}, \code{\link{flefko2}()},
//' \code{\link{aflefko2}()}, \code{\link{rlefko3}()}, and
//' \code{\link{rlefko2}()}.
//' 
//' @param stageframe The original stageframe.
//' @param supplement The original supplemental data input
//' (class \code{lefkoSD}). Can also equal NA.
//' @param overwrite An overwrite table.
//' @param repmatrix The original reproduction matrix. Can also equal NA or 0.
//' @param agemat A logical value indicating whether MPM is age-by-stage.
//' @param historical A logical value indicating whether MPM is historical.
//' @param format An integer indicating whether matrices will be in Ehrlen format
//' (if set to 1), or deVries format (if set to 2). Setting to deVries format
//' adds one extra stage to account for the prior status of newborns.
//' 
//' @return This function returns a list with a modified stageframe usable in MPM
//' construction, an associated reproduction matrix, and a general supplement
//' table that takes over the input supplement and overwrite tables. Note that
//' if a supplement is provided and a repmatrix is not, or if repmatrix is set
//' to 0, then it will be assumed that a repmatrix should not be used.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sf_reassess)]]
Rcpp::List sf_reassess(DataFrame stageframe, Nullable<DataFrame> supplement,
  Nullable<DataFrame> overwrite, Nullable<NumericMatrix> repmatrix,
  bool agemat = false, bool historical = false, int format = 1) {
  
  bool skiprepmat = false;
  bool supp_provided = false;
  bool repm_provided = false;
  bool over_provided = false;
  
  Rcpp::DataFrame supplement_true;
  arma::mat repmatrix_true;
  int newsupp_rows {0};
  
  StringVector stagevec = stageframe["stage"];
  NumericVector origsizevec = stageframe["size"];
  NumericVector origsizebvec = stageframe["size_b"];
  NumericVector origsizecvec = stageframe["size_c"];
  NumericVector minagevec = stageframe["min_age"];
  NumericVector maxagevec = stageframe["max_age"];
  arma::uvec repvec = stageframe["repstatus"];
  arma::uvec obsvec = stageframe["obsstatus"];
  arma::uvec propvec = stageframe["propstatus"];
  arma::uvec immvec = stageframe["immstatus"];
  arma::uvec matvec = stageframe["matstatus"];
  arma::uvec indvec = stageframe["indataset"];
  NumericVector binvec = stageframe["binhalfwidth_raw"];
  NumericVector binbvec = stageframe["binhalfwidthb_raw"];
  NumericVector bincvec = stageframe["binhalfwidthc_raw"];
  NumericVector sizeminvec = stageframe["sizebin_min"];
  NumericVector sizemaxvec = stageframe["sizebin_max"];
  NumericVector sizectrvec = stageframe["sizebin_center"];
  NumericVector sizewidthvec = stageframe["sizebin_width"];
  NumericVector sizebminvec = stageframe["sizebinb_min"];
  NumericVector sizebmaxvec = stageframe["sizebinb_max"];
  NumericVector sizebctrvec = stageframe["sizebinb_center"];
  NumericVector sizebwidthvec = stageframe["sizebinb_width"];
  NumericVector sizecminvec = stageframe["sizebinc_min"];
  NumericVector sizecmaxvec = stageframe["sizebinc_max"];
  NumericVector sizecctrvec = stageframe["sizebinc_center"];
  NumericVector sizecwidthvec = stageframe["sizebinc_width"];
  arma::ivec groupvec = stageframe["group"];
  StringVector comvec = stageframe["comments"];
  
  StringVector stage3_supp;
  StringVector stage2_supp;
  StringVector stage1_supp;
  StringVector eststage3_supp;
  StringVector eststage2_supp;
  StringVector eststage1_supp;
  NumericVector givenrate_supp;
  NumericVector multiplier_supp;
  IntegerVector convtype_supp;
  IntegerVector convtype_t12_supp;
  int supp_rows {0};
  
  int stageframe_length = repvec.n_elem;
  arma::uvec repentryvec(stageframe_length, fill::zeros);
  
  // Identify all groups in the stageframe
  arma::ivec all_groups = unique(groupvec);
  int no_groups = all_groups.n_elem;
  StringVector group_text(no_groups);
  
  for (int i = 0; i < no_groups; i++) {
    group_text(i) = "group";
    group_text(i) += std::to_string(all_groups(i));
  }
  
  if (supplement.isNotNull()) {
    supp_provided = true;
    
    Rcpp::DataFrame supplement_thru(supplement);
    supplement_true = supplement_thru;
    
    stage3_supp = supplement_true["stage3"];
    stage2_supp = supplement_true["stage2"];
    stage1_supp = supplement_true["stage1"];
    eststage3_supp = supplement_true["eststage3"];
    eststage2_supp = supplement_true["eststage2"];
    eststage1_supp = supplement_true["eststage1"];
    givenrate_supp = supplement_true["givenrate"];
    multiplier_supp = supplement_true["multiplier"]; // Not in overwrite()
    convtype_supp = supplement_true["convtype"];
    convtype_t12_supp = supplement_true["convtype_t12"];
    supp_rows = stage3_supp.length();
  }
  
  if (repmatrix.isNotNull()) {
    repm_provided = true;
    
    NumericMatrix repmatrix_thru(repmatrix);
    repmatrix_true = as<arma::mat>(repmatrix_thru);
    
    arma::rowvec rep_sums = sum(repmatrix_true, 0);
    arma::vec rep_rowsums = sum(repmatrix_true, 1);
    double repmat_sum = sum(rep_sums);
    
    if (rep_sums.n_elem != stageframe_length) {
      throw Rcpp::exception("Object repmatrix must have the rows and columns equal to the number of rows in the stageframe.", false);
    }
    
    if (repmat_sum > 0) {
      for (int i = 0; i < rep_sums.n_elem; i++) {
        if (rep_sums(i) > 0) {
          repvec(i) = 1;
        } else {
          repvec(i) = 0;
        }
        
        if (rep_rowsums(i) > 0) {
          repentryvec(i) = 1;
        } else {
          repentryvec(i) = 0;
        }
      }
    } else {
      repentryvec = propvec + immvec;
      int rev_checksum = sum(repentryvec);
      
      if (rev_checksum == 0) {
        repentryvec(0) = 1;
        Rf_warningcall(R_NilValue, "No information on reproductive entry stages provided. Assuming the first stage is the entry stage into the life cycle.");
      }
      
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
        
      arma::uvec repentry_calls = find(repentryvec);
      arma::uvec rep_calls = find(repvec);
      
      if (repentry_calls.n_elem > 0 && rep_calls.n_elem > 0) {
        for (int i = 0; i < repentry_calls.n_elem; i++) {
          for (int j = 0; j < rep_calls.n_elem; j++) {
            token_mat(repentry_calls(i), rep_calls(j)) = 1;
          }
        }
      }
      repmatrix_true = token_mat;
    }
  }  else if (supp_provided) {
    arma::ivec cv_supp_arma = as<arma::ivec>(convtype_supp);
    arma::uvec mult_elems = find(cv_supp_arma == 3);
    
    int fec_mults = mult_elems.n_elem;
    // arma::uvec needed_reprods(fec_mults, fill::zeros);
    for (int i = 0; i < fec_mults; i++) {
      for (int j = 0; j < stageframe_length; j++) {
        if (stage3_supp(mult_elems(i)) == stagevec(j)) {
          repentryvec(j) = 1;
        }
      }
      
      if (stage3_supp(mult_elems(i)) == "prop") {
        arma::uvec propvec_1elems = find(propvec);
        
        int propvec_no1s = propvec_1elems.n_elem;
        for (int j = 0; j < propvec_no1s; j++) {
          repentryvec(propvec_1elems(j)) = 1;
        }
      }
      
      if (stage3_supp(mult_elems(i)) == "immat") {
        arma::uvec immvec_1elems = find(immvec);
        
        int immvec_no1s = immvec_1elems.n_elem;
        for (int j = 0; j < immvec_no1s; j++) {
          repentryvec(immvec_1elems(j)) = 1;
        }
      }
      
      for (int j = 0; j < no_groups; j++) {
        if (stage3_supp(mult_elems(i)) == group_text(j)) {
          arma::uvec group_elems = find(groupvec == all_groups(j));
          
          int group_noelems = group_elems.n_elem;
          for (int k = 0; k < group_noelems; k++) {
            repentryvec(group_elems(k)) = 1;
          }
        }
      }
    }
    
    if (agemat) {
      arma::ivec needed_reprods(stageframe_length, fill::zeros);
      arma::vec needed_mults(stageframe_length, fill::zeros);
      
      for (int i = 0; i < fec_mults; i++) {
        for (int j = 0; j < stageframe_length; j++) {
          if (stage2_supp(mult_elems(i)) == stagevec(j)) {
            needed_reprods(j) = 1;
            needed_mults(j) = multiplier_supp(mult_elems(i));
          }
          
          if (stage2_supp(mult_elems(i)) == "rep") {
            arma::uvec repvec_1elems = find(repvec);
            
            int repvec_no1s = repvec_1elems.n_elem;
            for (int k = 0; k < repvec_no1s; k++) {
              needed_reprods(repvec_1elems(k)) = 1;
              needed_mults(repvec_1elems(k)) = multiplier_supp(mult_elems(i));
            }
          }
          
          for (int k = 0; k < no_groups; k++) {
            if (stage2_supp(mult_elems(i)) == group_text(k)) {
              arma::uvec group_elems = find(groupvec == all_groups(k));
              
              int group_noelems = group_elems.n_elem;
              for (int l = 0; l < group_noelems; l++) {
                needed_reprods(group_elems(l)) = 1;
                needed_mults(group_elems(l)) = multiplier_supp(mult_elems(i));
              }
            }
          }
        }
      }
      
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
      
      arma::uvec repentry_calls = find(repentryvec);
      arma::uvec neededrep_calls = find(needed_reprods);
      
      if (repentry_calls.n_elem > 0 && neededrep_calls.n_elem > 0) {
        for (int i = 0; i < repentry_calls.n_elem; i++) {
          for (int j = 0; j < neededrep_calls.n_elem; j++) {
            token_mat(repentry_calls(i), neededrep_calls(j)) = needed_mults(neededrep_calls(j));
          }
        }
      }
      repmatrix_true = token_mat;
      skiprepmat = false;
    
    } else {
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
      repmatrix_true = token_mat;
      skiprepmat = true;
    }
  } else {
    for (int i = 0; i < stageframe_length; i++) {
      if (propvec(i) > 0) {
        repentryvec(i) = 1;
      } else if (immvec(i) > 0) {
        repentryvec(i) = 1;
      }
    }
    
    int rev_checksum = sum(repentryvec);
    if(rev_checksum == 0) {
      repentryvec(0) = 1;
      Rf_warningcall(R_NilValue, "No information on reproductive entry stages provided. Assuming the first stage is the entry stage into the life cycle.");
    }
    
    arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
      
    arma::uvec repentry_calls = find(repentryvec);
    arma::uvec rep_calls = find(repvec);
    
    if (repentry_calls.n_elem > 0 && rep_calls.n_elem > 0) {
      for (int i = 0; i < repentry_calls.n_elem; i++) {
        for (int j = 0; j < rep_calls.n_elem; j++) {
          token_mat(repentry_calls(i), rep_calls(j)) = 1;
        }
      }
    }
    repmatrix_true = token_mat;
  }
  
  // Now we reorder the stageframe
  arma::uvec prop_stages = find(propvec);
  arma::uvec prop0_stages = find(propvec == 0);
  arma::uvec imm_stages = find(immvec);
  arma::uvec imm0_stages = find(immvec == 0);
  arma::uvec p0_im_stages = intersect(prop0_stages, imm_stages);
  
  arma::uvec mat_stages = find(matvec);
  arma::uvec mat_imm0_stages = intersect(mat_stages, imm0_stages);
  arma::uvec rep_stages = find(repvec);
  arma::uvec rep0_stages = find(repvec == 0);
  arma::uvec rep0_mat_imm0_stages = intersect(mat_imm0_stages, rep0_stages);
  arma::uvec mat_rep0_stages = intersect(mat_stages, rep0_stages);
  arma::uvec neworder(stageframe_length, fill::zeros);
  arma::uvec tracked_elems(stageframe_length, fill::ones);
  
  int no_prop_stages = prop_stages.n_elem;
  int no_p0_im_stages = p0_im_stages.n_elem;
  int no_r0_im_mt_stages = rep0_mat_imm0_stages.n_elem;
  
  int counter {0};
  for (int j = 0; j < no_prop_stages; j++) {
    neworder(counter) = prop_stages(j);
    tracked_elems(prop_stages(j)) = 0;
    counter++;
  }
  for (int j = 0; j < no_p0_im_stages; j++) {
    neworder(counter) = p0_im_stages(j);
    tracked_elems(p0_im_stages(j)) = 0;
    counter++;
  }
  for (int j = 0; j < no_r0_im_mt_stages; j++) {
    neworder(counter) = rep0_mat_imm0_stages(j);
    tracked_elems(rep0_mat_imm0_stages(j)) = 0;
    counter++;
  }
  
  arma::uvec remaining_elems = find(tracked_elems);
  int no_remaining_elems = remaining_elems.n_elem;
  for (int j = 0; j < no_remaining_elems; j++) {
    neworder(counter) = remaining_elems(j);
    counter++;
  }
  
  StringVector newstagevec(stageframe_length + format);
  NumericVector neworigsizevec(stageframe_length + format);
  NumericVector neworigsizebvec(stageframe_length + format);
  NumericVector neworigsizecvec(stageframe_length + format);
  NumericVector newminagevec(stageframe_length + format);
  NumericVector newmaxagevec(stageframe_length + format);
  arma::uvec newrepvec((stageframe_length + format), fill::zeros);
  arma::uvec newobsvec((stageframe_length + format), fill::zeros);
  arma::uvec newpropvec((stageframe_length + format), fill::zeros);
  arma::uvec newimmvec((stageframe_length + format), fill::zeros);
  arma::uvec newmatvec((stageframe_length + format), fill::zeros);
  arma::uvec newrepentryvec((stageframe_length + format), fill::zeros);
  arma::uvec newindvec((stageframe_length + format), fill::zeros);
  NumericVector newbinvec((stageframe_length + format));
  NumericVector newbinbvec((stageframe_length + format));
  NumericVector newbincvec((stageframe_length + format));
  NumericVector newsizeminvec((stageframe_length + format));
  NumericVector newsizemaxvec((stageframe_length + format));
  NumericVector newsizectrvec((stageframe_length + format));
  NumericVector newsizewidthvec((stageframe_length + format));
  NumericVector newsizebminvec((stageframe_length + format));
  NumericVector newsizebmaxvec((stageframe_length + format));
  NumericVector newsizebctrvec((stageframe_length + format));
  NumericVector newsizebwidthvec((stageframe_length + format));
  NumericVector newsizecminvec((stageframe_length + format));
  NumericVector newsizecmaxvec((stageframe_length + format));
  NumericVector newsizecctrvec((stageframe_length + format));
  NumericVector newsizecwidthvec((stageframe_length + format));
  arma::ivec newgroupvec((stageframe_length + format), fill::zeros);
  StringVector newcomvec((stageframe_length + format));
  arma::uvec newalive((stageframe_length + format), fill::zeros);
  NumericVector newalmostbornvec((stageframe_length + format));
  
  arma::mat repmat1(stageframe_length, stageframe_length);
  arma::mat repmat2(stageframe_length, stageframe_length);
  
  for (int i = 0; i < stageframe_length; i++) {
    newstagevec(i) = stagevec(neworder(i));
    neworigsizevec(i) = origsizevec(neworder(i));
    neworigsizebvec(i) = origsizebvec(neworder(i));
    neworigsizecvec(i) = origsizecvec(neworder(i));
    newminagevec(i) = minagevec(neworder(i));
    if (NumericVector::is_na(newminagevec(i))) newminagevec(i) = 0.0;
    newmaxagevec(i) = maxagevec(neworder(i));
    newrepvec(i) = repvec(neworder(i));
    newobsvec(i) = obsvec(neworder(i));
    newpropvec(i) = propvec(neworder(i));
    newimmvec(i) = immvec(neworder(i));
    newmatvec(i) = matvec(neworder(i));
    newrepentryvec(i) = repentryvec(neworder(i));
    newindvec(i) = indvec(neworder(i));
    newbinvec(i) = binvec(neworder(i));
    newbinbvec(i) = binbvec(neworder(i));
    newbincvec(i) = bincvec(neworder(i));
    newsizeminvec(i) = sizeminvec(neworder(i));
    newsizemaxvec(i) = sizemaxvec(neworder(i));
    newsizectrvec(i) = sizectrvec(neworder(i));
    newsizewidthvec(i) = sizewidthvec(neworder(i));
    newsizebminvec(i) = sizebminvec(neworder(i));
    newsizebmaxvec(i) = sizebmaxvec(neworder(i));
    newsizebctrvec(i) = sizebctrvec(neworder(i));
    newsizebwidthvec(i) = sizebwidthvec(neworder(i));
    newsizecminvec(i) = sizecminvec(neworder(i));
    newsizecmaxvec(i) = sizecmaxvec(neworder(i));
    newsizecctrvec(i) = sizecctrvec(neworder(i));
    newsizecwidthvec(i) = sizecwidthvec(neworder(i));
    newgroupvec(i) = groupvec(neworder(i));
    newcomvec(i) = comvec(neworder(i));
    newalive(i) = 1;
    
    repmat1.col(i) = repmatrix_true.col(neworder(i));
  }
  
  for (int i = 0; i < stageframe_length; i++) {
    repmat2.row(i) = repmat1.row(neworder(i));
  }
  
  if (format == 2) {
    newstagevec(stageframe_length) = "AlmostBorn";
    neworigsizevec(stageframe_length) = 0.0;
    neworigsizebvec(stageframe_length) = 0.0;
    neworigsizecvec(stageframe_length) = 0.0;
    newrepvec(stageframe_length) = 0;
    newobsvec(stageframe_length) = 1;
    newpropvec(stageframe_length) = 0;
    newimmvec(stageframe_length) = 1;
    newmatvec(stageframe_length) = 0;
    newrepentryvec(stageframe_length) = 0;
    newindvec(stageframe_length) = 1;
    newbinvec(stageframe_length) = 0.0;
    newbinbvec(stageframe_length) = 0.0;
    newbincvec(stageframe_length) = 0.0;
    newsizeminvec(stageframe_length) = 0.0;
    newsizemaxvec(stageframe_length) = 0.0;
    newsizectrvec(stageframe_length) = 0.0;
    newsizewidthvec(stageframe_length) = 0.0;
    newsizebminvec(stageframe_length) = 0.0;
    newsizebmaxvec(stageframe_length) = 0.0;
    newsizebctrvec(stageframe_length) = 0.0;
    newsizebwidthvec(stageframe_length) = 0.0;
    newsizecminvec(stageframe_length) = 0.0;
    newsizecmaxvec(stageframe_length) = 0.0;
    newsizecctrvec(stageframe_length) = 0.0;
    newsizecwidthvec(stageframe_length) = 0.0;
    newgroupvec(stageframe_length) = no_groups + 1;
    newcomvec(stageframe_length) = "Almost born (t-1)";
    newalive(stageframe_length) = 0;
    newalmostbornvec(stageframe_length) = 1.0;
    
    if (agemat) {
      newminagevec(stageframe_length) = 0.0;
      newmaxagevec(stageframe_length) = 0.0;
    } else {
      newminagevec(stageframe_length) = NA_REAL;
      newmaxagevec(stageframe_length) = NA_REAL;
    }
  }
  
  newstagevec(stageframe_length + (format - 1)) = "Dead";
  neworigsizevec(stageframe_length + (format - 1)) = 0.0;
  neworigsizebvec(stageframe_length + (format - 1)) = 0.0;
  neworigsizecvec(stageframe_length + (format - 1)) = 0.0;
  newrepvec(stageframe_length + (format - 1)) = 0;
  newobsvec(stageframe_length + (format - 1)) = 1;
  newpropvec(stageframe_length + (format - 1)) = 0;
  newimmvec(stageframe_length + (format - 1)) = 0;
  newmatvec(stageframe_length + (format - 1)) = 1;
  newrepentryvec(stageframe_length + (format - 1)) = 0;
  newindvec(stageframe_length + (format - 1)) = 1;
  newbinvec(stageframe_length + (format - 1)) = 0.0;
  newbinbvec(stageframe_length + (format - 1)) = 0.0;
  newbincvec(stageframe_length + (format - 1)) = 0.0;
  newsizeminvec(stageframe_length + (format - 1)) = 0.0;
  newsizemaxvec(stageframe_length + (format - 1)) = 0.0;
  newsizectrvec(stageframe_length + (format - 1)) = 0.0;
  newsizewidthvec(stageframe_length + (format - 1)) = 0.0;
  newsizebminvec(stageframe_length + (format - 1)) = 0.0;
  newsizebmaxvec(stageframe_length + (format - 1)) = 0.0;
  newsizebctrvec(stageframe_length + (format - 1)) = 0.0;
  newsizebwidthvec(stageframe_length + (format - 1)) = 0.0;
  newsizecminvec(stageframe_length + (format - 1)) = 0.0;
  newsizecmaxvec(stageframe_length + (format - 1)) = 0.0;
  newsizecctrvec(stageframe_length + (format - 1)) = 0.0;
  newsizecwidthvec(stageframe_length + (format - 1)) = 0.0;
  newgroupvec(stageframe_length + (format - 1)) = no_groups + 1;
  newcomvec(stageframe_length + (format - 1)) = "Dead";
  newalive(stageframe_length + (format - 1)) = 0;
  
  if (agemat) {
    newminagevec(stageframe_length + (format - 1)) = 0.0;
    newmaxagevec(stageframe_length + (format - 1)) = 0.0;
  } else {
    newminagevec(stageframe_length + (format - 1)) = NA_REAL;
    newmaxagevec(stageframe_length + (format - 1)) = NA_REAL;
  }
  
  arma::ivec stage_id;
  if (format == 1) {
    stage_id = linspace<arma::ivec>(1, (stageframe_length + 1), (stageframe_length + 1));
  } else {
    stage_id = linspace<arma::ivec>(1, (stageframe_length + 2), (stageframe_length + 2));
  }
  
  StringVector newstagevec_check = unique(newstagevec);
  if (newstagevec_check.length() < newstagevec.length()) {
    throw Rcpp::exception("All stage names must be unique.", false);
  }
  
  Rcpp::List newstageframe(33);
  
  newstageframe(0) = Rcpp::IntegerVector(stage_id.begin(), stage_id.end());
  newstageframe(1) = newstagevec;
  newstageframe(2) = neworigsizevec;
  newstageframe(3) = neworigsizebvec;
  newstageframe(4) = neworigsizecvec;
  newstageframe(5) = newminagevec;
  newstageframe(6) = newmaxagevec;
  newstageframe(7) = Rcpp::IntegerVector(newrepvec.begin(), newrepvec.end());
  newstageframe(8) = Rcpp::IntegerVector(newobsvec.begin(), newobsvec.end());
  newstageframe(9) = Rcpp::IntegerVector(newpropvec.begin(), newpropvec.end());
  newstageframe(10) = Rcpp::IntegerVector(newimmvec.begin(), newimmvec.end());
  newstageframe(11) = Rcpp::IntegerVector(newmatvec.begin(), newmatvec.end());
  newstageframe(12) = Rcpp::IntegerVector(newrepentryvec.begin(), newrepentryvec.end());
  newstageframe(13) = Rcpp::IntegerVector(newindvec.begin(), newindvec.end());
  newstageframe(14) = newbinvec;
  newstageframe(15) = newsizeminvec;
  newstageframe(16) = newsizemaxvec;
  newstageframe(17) = newsizectrvec;
  newstageframe(18) = newsizewidthvec;
  newstageframe(19) = newbinbvec;
  newstageframe(20) = newsizebminvec;
  newstageframe(21) = newsizebmaxvec;
  newstageframe(22) = newsizebctrvec;
  newstageframe(23) = newsizebwidthvec;
  newstageframe(24) = newbincvec;
  newstageframe(25) = newsizecminvec;
  newstageframe(26) = newsizecmaxvec;
  newstageframe(27) = newsizecctrvec;
  newstageframe(28) = newsizecwidthvec;
  newstageframe(29) = Rcpp::IntegerVector(newgroupvec.begin(), newgroupvec.end());
  newstageframe(30) = newcomvec;
  newstageframe(31) = Rcpp::IntegerVector(newalive.begin(), newalive.end());
  newstageframe(32) = newalmostbornvec; // For use in deVries-format hMPMs
  
  CharacterVector namevec = {"stage_id", "stage", "original_size", "original_size_b",
    "original_size_c", "min_age", "max_age", "repstatus", "obsstatus", "propstatus",
    "immstatus", "matstatus", "entrystage", "indataset", "binhalfwidth_raw", "sizebin_min",
    "sizebin_max", "sizebin_center", "sizebin_width", "binhalfwidthb_raw", "sizebinb_min",
    "sizebinb_max", "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments", "alive",
    "almostborn"};
  CharacterVector newclasses = {"data.frame", "stageframe"};
  newstageframe.attr("names") = namevec;
  newstageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, (stageframe_length + format));
  newstageframe.attr("class") = newclasses;
  
  // From here on, we cover the new version of .overwrite_reassess
  if (overwrite.isNotNull()) {
    if (supp_provided) {
      throw Rcpp::exception("Please input either a supplement table or an overwrite table, not both.", false);
    }
    
    Rcpp::DataFrame supplement_thru(supplement);
    supplement_true = supplement_thru;
    
    stage3_supp = supplement_true["stage3"];
    stage2_supp = supplement_true["stage2"];
    stage1_supp = supplement_true["stage1"];
    eststage3_supp = supplement_true["eststage3"];
    eststage2_supp = supplement_true["eststage2"];
    eststage1_supp = supplement_true["eststage1"];
    givenrate_supp = supplement_true["givenrate"];
    convtype_supp = supplement_true["convtype"];
    convtype_t12_supp = supplement_true["convtype_t12"];
    
    supp_rows = givenrate_supp.length();
    multiplier_supp = Rcpp::NumericVector::create(1.0, supp_rows); // Not in overwrite()
    
    over_provided = true;
  }
  
  Rcpp::List newsupplement(10);
  
  if (over_provided || supp_provided) {
    StringVector unique_stages = unique(newstagevec);
    StringVector extra_terms = {"rep", "nrep", "immat", "mat", "prop", "npr", "all", "obs", "nobs"};
    
    int no_newstages = unique_stages.length();
    int no_extraterms = extra_terms.length();
    
    StringVector all_possible_stage_terms(no_newstages + no_extraterms + no_groups);
    for (int i = 0; i < no_newstages; i++) {
      all_possible_stage_terms(i) = unique_stages(i);
    }
    for (int i = 0; i < no_extraterms; i++) {
      all_possible_stage_terms(i + no_newstages) = extra_terms(i);
    }
    for (int i = 0; i < no_groups; i++) {
      all_possible_stage_terms(i + no_newstages + no_extraterms) = group_text(i);
    }
    
    // Check for good entries in the supplement or overwrite table
    for (int i = 0; i < stage3_supp.length(); i++) {
      int s3supp_count {0};
      int s2supp_count {0};
      int s1supp_count {0};
      
      bool ests3_used = false;
      bool ests2_used = false;
      bool ests1_used = false;
      
      int ests3supp_count {0};
      int ests2supp_count {0};
      int ests1supp_count {0};
      
      for (int j = 0; j < all_possible_stage_terms.length(); j++) {
        if (stage3_supp(i) == all_possible_stage_terms(j)) s3supp_count++;
        if (stage2_supp(i) == all_possible_stage_terms(j)) s2supp_count++;
        
        if (!StringVector::is_na(eststage3_supp(i))) {
          ests3_used = true;
          if (eststage3_supp(i) == all_possible_stage_terms(j)) ests3supp_count++;
        }
        if (!StringVector::is_na(eststage2_supp(i))) {
          ests2_used = true;
          if (eststage2_supp(i) == all_possible_stage_terms(j)) ests2supp_count++;
        }
        
        if (historical) {
          if (stage1_supp(i) == all_possible_stage_terms(j)) s1supp_count++;
          
          if (!StringVector::is_na(eststage1_supp(i))) {
            ests1_used = true;
            if (eststage1_supp(i) == all_possible_stage_terms(j)) ests1supp_count++;
          }
        } 
      }
      
      if (s3supp_count == 0) {
        throw Rcpp::exception("Stage names in supplement or overwrite table (stage3) must match stageframe",
          false);
      }
      if (s2supp_count == 0) {
        throw Rcpp::exception("Stage names in supplement or overwrite table (stage2) must match stageframe",
          false);
      }
      if (ests3_used) {
        if (s3supp_count == 0) {
          throw Rcpp::exception("Stage names in supplement or overwrite table (eststage3) must match stageframe",
            false);
        }
      }
      if (ests2_used) {
        if (s2supp_count == 0) {
          throw Rcpp::exception("Stage names in supplement or overwrite table (eststage2) must match stageframe",
            false);
        }
      }
      if (historical) {
        if (s1supp_count == 0) {
          throw Rcpp::exception("Stage names in supplement or overwrite table (stage1) must match stageframe",
            false);
        }
        if (ests1_used) {
          if (s1supp_count == 0) {
            throw Rcpp::exception("Stage names in supplement or overwrite table (eststage1) must match stageframe",
              false);
          }
        }
      }
    }
    
    IntegerVector s1_calls (supp_rows, 1);
    IntegerVector s2_calls (supp_rows, 1);
    IntegerVector s3_calls (supp_rows, 1);
    IntegerVector ests1_calls (supp_rows, 1);
    IntegerVector ests2_calls (supp_rows, 1);
    IntegerVector ests3_calls (supp_rows, 1);
    IntegerVector s3_planned (supp_rows, 1);
    IntegerVector s2_planned (supp_rows, 1);
    IntegerVector s1_planned (supp_rows, 1);
      
    IntegerVector s123_calls (supp_rows, 1);
    
    // This section creates the new indices for the edited supplement/overwrite table
    arma::uvec newprop_stages = find(newpropvec);
    arma::uvec newprop0_stages = find(newpropvec == 0);
    arma::uvec newimm_stages = find(newimmvec);
    arma::uvec newalive_stages = find(newalive);
    arma::uvec newmat_stages1 = find(newmatvec);
    arma::uvec newmat_stages = intersect(newalive_stages, newmat_stages1);
    arma::uvec newrep_stages = find(newrepvec);
    arma::uvec newrep0_stages = find(newrepvec == 0);
    arma::uvec newmat_rep0_stages = intersect(newmat_stages, newrep0_stages);
    arma::uvec newobs_stages = find(newobsvec);
    arma::uvec newobs0_stages = find(newobsvec == 0);
    arma::uvec all_stages = find(newalive); // 7 "all"
    int no_current_group {0};
    
    // Now we build the expanded and edited supplement table
    for (int i = 0; i < supp_rows; i++) {
      if (stage3_supp(i) == "prop") {
        s3_calls(i) = newprop_stages.n_elem;
      } else if (stage3_supp(i) == "npr") {
        s3_calls(i) = newprop0_stages.n_elem;
      } else if (stage3_supp(i) == "immat") {
        s3_calls(i) = newimm_stages.n_elem;
      } else if (stage3_supp(i) == "mat") {
        s3_calls(i) = newmat_stages.n_elem;
      } else if (stage3_supp(i) == "rep") {
        s3_calls(i) = newrep_stages.n_elem;
      } else if (stage3_supp(i) == "nrep") {
        s3_calls(i) = newmat_rep0_stages.n_elem;
      } else if (stage3_supp(i) == "obs") {
        s3_calls(i) = newobs_stages.n_elem;
      } else if (stage3_supp(i) == "nobs") {
        s3_calls(i) = newobs0_stages.n_elem;
      } else if (stage3_supp(i) == "all") {
        s3_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (stage3_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            s3_calls(i) = no_current_group;
          }
        }
      }
      if (s3_calls(i) == 0) s3_calls(i) = 1;
      
      if (eststage3_supp(i) == "prop") {
        ests3_calls(i) = newprop_stages.n_elem;
      } else if (eststage3_supp(i) == "npr") {
        ests3_calls(i) = newprop0_stages.n_elem;
      } else if (eststage3_supp(i) == "immat") {
        ests3_calls(i) = newimm_stages.n_elem;
      } else if (eststage3_supp(i) == "mat") {
        ests3_calls(i) = newmat_stages.n_elem;
      } else if (eststage3_supp(i) == "rep") {
        ests3_calls(i) = newrep_stages.n_elem;
      } else if (eststage3_supp(i) == "nrep") {
        ests3_calls(i) = newmat_rep0_stages.n_elem;
      } else if (eststage3_supp(i) == "obs") {
        ests3_calls(i) = newobs_stages.n_elem;
      } else if (eststage3_supp(i) == "nobs") {
        ests3_calls(i) = newobs0_stages.n_elem;
      } else if (eststage3_supp(i) == "all") {
        ests3_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (eststage3_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            ests3_calls(i) = no_current_group;
          }
        }
      }
      if (ests3_calls(i) == 0) ests3_calls(i) = 1;
      
      if (stage2_supp(i) == "prop") {
        s2_calls(i) = newprop_stages.n_elem;
      } else if (stage2_supp(i) == "npr") {
        s2_calls(i) = newprop0_stages.n_elem;
      } else if (stage2_supp(i) == "immat") {
        s2_calls(i) = newimm_stages.n_elem;
      } else if (stage2_supp(i) == "mat") {
        s2_calls(i) = newmat_stages.n_elem;
      } else if (stage2_supp(i) == "rep") {
        s2_calls(i) = newrep_stages.n_elem;
      } else if (stage2_supp(i) == "nrep") {
        s2_calls(i) = newmat_rep0_stages.n_elem;
      } else if (stage2_supp(i) == "obs") {
        s2_calls(i) = newobs_stages.n_elem;
      } else if (stage2_supp(i) == "nobs") {
        s2_calls(i) = newobs0_stages.n_elem;
      } else if (stage2_supp(i) == "all") {
        s2_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (stage2_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            s2_calls(i) = no_current_group;
          }
        }
      }
      if (s2_calls(i) == 0) s2_calls(i) = 1;
      
      if (eststage2_supp(i) == "prop") {
        ests2_calls(i) = newprop_stages.n_elem;
      } else if (eststage2_supp(i) == "npr") {
        ests2_calls(i) = newprop0_stages.n_elem;
      } else if (eststage2_supp(i) == "immat") {
        ests2_calls(i) = newimm_stages.n_elem;
      } else if (eststage2_supp(i) == "mat") {
        ests2_calls(i) = newmat_stages.n_elem;
      } else if (eststage2_supp(i) == "rep") {
        ests2_calls(i) = newrep_stages.n_elem;
      } else if (eststage2_supp(i) == "nrep") {
        ests2_calls(i) = newmat_rep0_stages.n_elem;
      } else if (eststage2_supp(i) == "obs") {
        ests2_calls(i) = newobs_stages.n_elem;
      } else if (eststage2_supp(i) == "nobs") {
        ests2_calls(i) = newobs0_stages.n_elem;
      } else if (eststage2_supp(i) == "all") {
        ests2_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (eststage2_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            ests2_calls(i) = no_current_group;
          }
        }
      }
      if (ests2_calls(i) == 0) ests2_calls(i) = 1;
      
      if (stage1_supp(i) == "prop") {
        s1_calls(i) = newprop_stages.n_elem;
      } else if (stage1_supp(i) == "npr") {
        s1_calls(i) = newprop0_stages.n_elem;
      } else if (stage1_supp(i) == "immat") {
        s1_calls(i) = newimm_stages.n_elem;
      } else if (stage1_supp(i) == "mat") {
        s1_calls(i) = newmat_stages.n_elem;
      } else if (stage1_supp(i) == "rep") {
        s1_calls(i) = newrep_stages.n_elem;
      } else if (stage1_supp(i) == "nrep") {
        s1_calls(i) = newmat_rep0_stages.n_elem;
      } else if (stage1_supp(i) == "obs") {
        s1_calls(i) = newobs_stages.n_elem;
      } else if (stage1_supp(i) == "nobs") {
        s1_calls(i) = newobs0_stages.n_elem;
      } else if (stage1_supp(i) == "all") {
        s1_calls(i) = all_stages.n_elem;
      } else if (StringVector::is_na(stage1_supp(i))) {
        s1_calls(i) = 1;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (stage1_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            s1_calls(i) = no_current_group;
          }
        }
      }
      if (s1_calls(i) == 0) s1_calls(i) = 1;
      
      if (eststage1_supp(i) == "prop") {
        ests1_calls(i) = newprop_stages.n_elem;
      } else if (eststage1_supp(i) == "npr") {
        ests1_calls(i) = newprop0_stages.n_elem;
      } else if (eststage1_supp(i) == "immat") {
        ests1_calls(i) = newimm_stages.n_elem;
      } else if (eststage1_supp(i) == "mat") {
        ests1_calls(i) = newmat_stages.n_elem;
      } else if (eststage1_supp(i) == "rep") {
        ests1_calls(i) = newrep_stages.n_elem;
      } else if (eststage1_supp(i) == "nrep") {
        ests1_calls(i) = newmat_rep0_stages.n_elem;
      } else if (eststage1_supp(i) == "obs") {
        ests1_calls(i) = newobs_stages.n_elem;
      } else if (eststage1_supp(i) == "nobs") {
        ests1_calls(i) = newobs0_stages.n_elem;
      } else if (eststage1_supp(i) == "all") {
        ests1_calls(i) = all_stages.n_elem;
      } else if (StringVector::is_na(eststage1_supp(i))) {
        ests1_calls(i) = 1;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (eststage1_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            ests1_calls(i) = no_current_group;
          }
        }
      }
      if (ests1_calls(i) == 0) ests1_calls(i) = 1;
      
      if (!StringVector::is_na(eststage3_supp(i))) {
        if (eststage3_supp(i) != stage3_supp(i)) {
          if (s3_calls(i) == 1 && ests3_calls(i) > 1) {
            s3_planned(i) = ests3_calls(i);
          } else if (s3_calls(i) > 1 && ests3_calls(i) > 1) {
            throw Rcpp::exception("If stage group shorthand is used to designate both a transition and a proxy, then the shorthand group must be the same in both cases.", false);
          }
        } else {
          s3_planned(i) = s3_calls(i);
        }
      } else {
        s3_planned(i) = s3_calls(i);
      }
      
      if (!StringVector::is_na(eststage2_supp(i))) {
        if (eststage2_supp(i) != stage2_supp(i)) {
          if (s2_calls(i) == 1 && ests2_calls(i) > 1) {
            s2_planned(i) = ests2_calls(i);
          } else if (s2_calls(i) > 1 && ests2_calls(i) > 1) {
            throw Rcpp::exception("If stage group shorthand is used to designate both a transition and a proxy, then the shorthand group must be the same in both cases.", false);
          }
        } else {
          s2_planned(i) = s2_calls(i);
        }
      } else {
        s2_planned(i) = s2_calls(i);
      }
      
      if (!StringVector::is_na(eststage1_supp(i))) {
        if (historical && eststage1_supp(i) != stage1_supp(i)) {
          if (s1_calls(i) == 1 && ests1_calls(i) > 1) {
            s1_planned(i) = ests1_calls(i);
          } else if (s1_calls(i) > 1 && ests1_calls(i) > 1) {
            throw Rcpp::exception("If stage group shorthand is used to designate both a transition and a proxy, then the shorthand group must be the same in both cases.", false);
          }
        } else if (historical) {
          s1_planned(i) = s1_calls(i);
        } else if (!historical) {
          s1_planned(i) = 1;
        }
      } else {
        s1_planned(i) = s1_calls(i);
      }
      
      s123_calls(i) = s3_planned(i) * s2_planned(i) * s1_planned(i);
    }
    
    NumericVector basepoints(supp_rows, 0.0);
    for (int i = 0; i < (supp_rows - 1); i++) {
      basepoints(i+1) = basepoints(i) + s123_calls(i);
    }
    
    newsupp_rows = sum(s123_calls);
    
    StringVector stage3_newsupp(newsupp_rows);
    StringVector stage2_newsupp(newsupp_rows);
    StringVector stage1_newsupp(newsupp_rows);
    StringVector eststage3_newsupp(newsupp_rows);
    StringVector eststage2_newsupp(newsupp_rows);
    StringVector eststage1_newsupp(newsupp_rows);
    NumericVector givenrate_newsupp(newsupp_rows);
    IntegerVector convtype_newsupp(newsupp_rows);
    IntegerVector convtype_t12_newsupp(newsupp_rows);
    NumericVector multiplier_newsupp(newsupp_rows);
    
    int overall_counter {0};
    int group_check {0};
    
    int group_baseline3 {0};
    int group_baseline2 {0};
    int group_baseline1 {0};
    int group_baselinee3 {0};
    int group_baselinee2 {0};
    int group_baselinee1 {0};
    
    int group_ratchet3 {0};
    int group_ratchet2 {0};
    int group_ratchet1 {0};
    int group_ratchete3 {0};
    int group_ratchete2 {0};
    int group_ratchete1 {0};
    
    int prevl3 {0};
    int prevl2 {0};
    int prevl1 {0};
    int prevle3 {0};
    int prevle2 {0};
    int prevle1 {0};
    
    for (int i = 0; i < supp_rows; i++) {
      overall_counter = 0;
      for (int j = 0; j < s1_planned(i); j++) {
        for (int k = 0; k < s2_planned(i); k++) {
          for (int l = 0; l < s3_planned(i); l++) {
            if (stage3_supp(i) == "prop") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(l));
            } else if (stage3_supp(i) == "npr") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(l));
            } else if (stage3_supp(i) == "immat") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(l));
            } else if (stage3_supp(i) == "mat") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(l));
            } else if (stage3_supp(i) == "rep") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(l));
            } else if (stage3_supp(i) == "nrep") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(l));
            } else if (stage3_supp(i) == "obs") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs_stages(l));
            } else if (stage3_supp(i) == "nobs") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs0_stages(l));
            } else if (stage3_supp(i) == "all") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(l));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (stage3_supp(i) == group_text(m)) {
                  if (l == 0) group_ratchet3 = 0;
                  if (l != prevl3 && l != 0) group_ratchet3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet3 > (current_group_length - 1)) {
                    group_ratchet3 = 0;
                  }
                  
                  if (group_ratchet3 == 0) {
                    group_baseline3 = l;
                  }
                  
                  stage3_newsupp(basepoints(i) + overall_counter) = 
                    newstagevec(current_group(l - group_baseline3));
                  
                  prevl3 = l;
                }
              }
              
              if (group_check == 0) {
               stage3_newsupp(basepoints(i) + overall_counter) = stage3_supp(i);
              }
              
              group_check = 0;
            }
            givenrate_newsupp(basepoints(i) + overall_counter) = givenrate_supp(i);
            multiplier_newsupp(basepoints(i) + overall_counter) = multiplier_supp(i);
            convtype_newsupp(basepoints(i) + overall_counter) = convtype_supp(i);
            convtype_t12_newsupp(basepoints(i) + overall_counter) = convtype_t12_supp(i);
            
            if (eststage3_supp(i) == "prop") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(l));
            } else if (eststage3_supp(i) == "npr") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(l));
            } else if (eststage3_supp(i) == "immat") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(l));
            } else if (eststage3_supp(i) == "mat") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(l));
            } else if (eststage3_supp(i) == "rep") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(l));
            } else if (eststage3_supp(i) == "nrep") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(l));
            } else if (eststage3_supp(i) == "obs") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs_stages(l));
            } else if (eststage3_supp(i) == "nobs") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs0_stages(l));
            } else if (eststage3_supp(i) == "all") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(l));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (eststage3_supp(i) == group_text(m)) {
                  if (l == 0) group_ratchete3 = 0;
                  if (l != prevle3 && l != 0) group_ratchete3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchete3 > (current_group_length - 1)) {
                    group_ratchete3 = 0;
                  }
                  
                  if (group_ratchete3 == 0) {
                    group_baselinee3 = l;
                  }
                  
                  eststage3_newsupp(basepoints(i) + overall_counter) = 
                    newstagevec(current_group(l - group_baselinee3));
                  
                  prevle3 = l;
                }
              }
              
              if (group_check == 0) {
                eststage3_newsupp(basepoints(i) + overall_counter) = eststage3_supp(i);
              }
              
              group_check = 0;
            }
            
            if (stage2_supp(i) == "prop") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(k));
            } else if (stage2_supp(i) == "npr") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(k));
            } else if (stage2_supp(i) == "immat") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(k));
            } else if (stage2_supp(i) == "mat") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(k));
            } else if (stage2_supp(i) == "rep") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(k));
            } else if (stage2_supp(i) == "nrep") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(k));
            } else if (stage2_supp(i) == "obs") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs_stages(k));
            } else if (stage2_supp(i) == "nobs") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs0_stages(k));
            } else if (stage2_supp(i) == "all") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(k));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (stage2_supp(i) == group_text(m)) {
                  if (k == 0) group_ratchet2 = 0;
                  if (k != prevl2 && k != 0) group_ratchet2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet2 > (current_group_length - 1)) {
                    group_ratchet2 = 0;
                  }
                  
                  if (group_ratchet2 == 0) {
                    group_baseline2 = k;
                  }
                  
                  stage2_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(k - group_baseline2));
                  
                  prevl2 = k;
                }
              }
              
              if (group_check == 0) {
               stage2_newsupp(basepoints(i) + overall_counter) = stage2_supp(i);
              }
              
              group_check = 0;
            }
            
            if (eststage2_supp(i) == "prop") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(k));
            } else if (eststage2_supp(i) == "npr") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(k));
            } else if (eststage2_supp(i) == "immat") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(k));
            } else if (eststage2_supp(i) == "mat") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(k));
            } else if (eststage2_supp(i) == "rep") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(k));
            } else if (eststage2_supp(i) == "nrep") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(k));
            } else if (eststage2_supp(i) == "obs") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs_stages(k));
            } else if (eststage2_supp(i) == "nobs") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs0_stages(k));
            } else if (eststage2_supp(i) == "all") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(k));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (eststage2_supp(i) == group_text(m)) {
                  if (k == 0) group_ratchete2 = 0;
                  if (k != prevle2 && k != 0) group_ratchete2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchete2 > (current_group_length - 1)) {
                    group_ratchete2 = 0;
                  }
                  
                  if (group_ratchete2 == 0) {
                    group_baselinee2 = k;
                  }
                  
                  eststage2_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(k - group_baselinee2));
                  
                  prevle2 = k;
                }
              }
              
              if (group_check == 0) {
                eststage2_newsupp(basepoints(i) + overall_counter) = eststage2_supp(i);
              }
              
              group_check = 0;
            }
            
            if (stage1_supp(i) == "prop") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(j));
            } else if (stage1_supp(i) == "npr") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(j));
            } else if (stage1_supp(i) == "immat") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(j));
            } else if (stage1_supp(i) == "mat") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(j));
            } else if (stage1_supp(i) == "rep") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(j));
            } else if (stage1_supp(i) == "nrep") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(j));
            } else if (stage1_supp(i) == "obs") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs_stages(j));
            } else if (stage1_supp(i) == "nobs") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs0_stages(j));
            } else if (stage1_supp(i) == "all") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(j));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (stage1_supp(i) == group_text(m)) {
                  if (j == 0) group_ratchet1 = 0;
                  if (j != prevl1 && j != 0) group_ratchet1 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet1 > (current_group_length - 1)) {
                    group_ratchet1 = 0;
                  }
                  
                  if (group_ratchet1 == 0) {
                    group_baseline1 = j;
                  }
                  
                  stage1_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(j - group_baseline1));
                  
                  prevl1 = j;
                }
              }
              
              if (group_check == 0) {
               stage1_newsupp(basepoints(i) + overall_counter) = stage1_supp(i);
              }
              
              group_check = 0;
            }
            
            if (eststage1_supp(i) == "prop") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(j));
            } else if (eststage1_supp(i) == "npr") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(j));
            } else if (eststage1_supp(i) == "immat") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(j));
            } else if (eststage1_supp(i) == "mat") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(j));
            } else if (eststage1_supp(i) == "rep") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(j));
            } else if (eststage1_supp(i) == "nrep") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(j));
            } else if (eststage1_supp(i) == "obs") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs_stages(j));
            } else if (eststage1_supp(i) == "nobs") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newobs0_stages(j));
            } else if (eststage1_supp(i) == "all") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(j));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (eststage1_supp(i) == group_text(m)) {
                  if (j == 0) group_ratchete1 = 0;
                  if (j != prevle1 && j != 0) group_ratchete1 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchete1 > (current_group_length - 1)) {
                    group_ratchete1 = 0;
                  }
                  
                  if (group_ratchete1 == 0) {
                    group_baselinee1 = j;
                  }
                  
                  eststage1_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(j - group_baselinee1));
                  
                  prevle1 = j;
                }
              }
              
              if (group_check == 0) {
                eststage1_newsupp(basepoints(i) + overall_counter) = eststage1_supp(i);
              }
              
              group_check = 0;
            }
            
            overall_counter++;
          }
        }
      }
    }
    
    newsupplement(0) = stage3_newsupp;
    newsupplement(1) = stage2_newsupp;
    newsupplement(2) = stage1_newsupp;
    newsupplement(3) = eststage3_newsupp;
    newsupplement(4) = eststage2_newsupp;
    newsupplement(5) = eststage1_newsupp;
    newsupplement(6) = givenrate_newsupp;
    newsupplement(7) = multiplier_newsupp;
    newsupplement(8) = convtype_newsupp;
    newsupplement(9) = convtype_t12_newsupp;
  } else {
    newsupplement(0) = NULL;
    newsupplement(1) = NULL;
    newsupplement(2) = NULL;
    newsupplement(3) = NULL;
    newsupplement(4) = NULL;
    newsupplement(5) = NULL;
    newsupplement(6) = NULL;
    newsupplement(7) = NULL;
    newsupplement(8) = NULL;
    newsupplement(9) = NULL;
  }
  
  CharacterVector su_namevec = {"stage3", "stage2", "stage1", "eststage3", "eststage2",
    "eststage1", "givenrate", "multiplier", "convtype", "convtype_t12"};
  CharacterVector su_newclasses = {"data.frame", "lefkoSD"};
  newsupplement.attr("names") = su_namevec;
  newsupplement.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newsupp_rows);
  newsupplement.attr("class") = su_newclasses;
  
  return Rcpp::List::create(_["stageframe"] = newstageframe, _["repmatrix"] = repmat2,
    _["ovtable"] = newsupplement);
}

//' Create Stageframe for Population Matrix Projection Analysis
//' 
//' Function \code{sf_leslie()} returns a data frame describing each age in a
//' Leslie MPM in terms of ahistorical stage information. This function is
//' internal to \code{\link{rleslie}()} and \code{\link{fleslie}()}.
//' 
//' @param min_age The first age to include in the matrix.
//' @param max_age The maximum age to include in the matrix.
//' @param min_fecage The first age in which reproduction is possible.
//' @param max_fecage The final age in which reproduction is possible.
//' @param cont A logical value indicating whether survival continues past the
//' last described age.
//' 
//' @return A data frame of class \code{stageframe}, which includes information
//' on the stage name, size, reproductive status, observation status, propagule 
//' status, immaturity status, maturity status, presence within the core dataset, 
//' stage group classification, raw bin half-width, and the minimum, 
//' center, and maximum of each size bin, as well as its width. If minimum and
//' maximum ages were specified, then these are also included. Also includes an 
//' empty string variable that can be used to describe stages meaningfully.
//' 
//' Variables in this data frame include the following:
//' \item{stage}{The unique names of the stages to be analyzed.}
//' \item{size}{The typical or representative size at which each stage occurs.}
//' \item{size_b}{Size at which each stage occurs in terms of a second size
//' variable, if one exists.}
//' \item{size_c}{Size at which each stage occurs in terms of a third size
//' variable, if one exists.}
//' \item{min_age}{The minimum age at which the stage may occur.}
//' \item{max_age}{The maximum age at which the stage may occur.}
//' \item{repstatus}{A binomial variable showing whether each stage is
//' reproductive.}
//' \item{obsstatus}{A binomial variable showing whether each stage is
//' observable.}
//' \item{propstatus}{A binomial variable showing whether each stage is a
//' propagule.}
//' \item{immstatus}{A binomial variable showing whether each stage can occur as
//' immature.}
//' \item{matstatus}{A binomial variable showing whether each stage occurs in
//' maturity.}
//' \item{indataset}{A binomial variable describing whether each stage occurs in
//' the input dataset.}
//' \item{binhalfwidth_raw}{The half-width of the size bin, as input.}
//' \item{sizebin_min}{The minimum size at which the stage may occur.}
//' \item{sizebin_max}{The maximum size at which the stage may occur.}
//' \item{sizebin_center}{The midpoint of the size bin at which the stage may
//' occur.}
//' \item{sizebin_width}{The width of the size bin corresponding to the stage.}
//' \item{binhalfwidthb_raw}{The half-width of the size bin of a second size
//' variable, as input.}
//' \item{sizebinb_min}{The minimum size at which the stage may occur.}
//' \item{sizebinb_max}{The maximum size at which the stage may occur.}
//' \item{sizebinb_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a second size variable.}
//' \item{sizebinb_width}{The width of the size bin corresponding to the stage,
//' in terms of a second size variable.}
//' \item{binhalfwidthc_raw}{The half-width of the size bin of a third size
//' variable, as input.}
//' \item{sizebinc_min}{The minimum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_max}{The maximum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a third size variable.}
//' \item{sizebinc_width}{The width of the size bin corresponding to the stage,
//' in terms of a third size variable.}
//' \item{group}{An integer denoting the size classification group that the
//' stage falls within.}
//' \item{comments}{A text field for stage descriptions.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sf_leslie)]]
Rcpp::List sf_leslie (int min_age, int max_age, int min_fecage,int max_fecage,
  bool cont) {
  
  int age_range = max_age - min_age;
  int total_ages = age_range + 1;
  if (age_range < 1) {
    throw Rcpp::exception("There must be at least two ages to create a Leslie MPM.", false);
  }
  
  Rcpp::List output_longlist(29);
  Rcpp::CharacterVector varnames {"stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center", "sizebin_width",
    "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
    "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
    "group", "comments"};
  
  int matsize = total_ages;
  
  StringVector agenames_true (matsize, NA_STRING);
  NumericVector size_true (matsize, NA_REAL);
  NumericVector sizesb_true (matsize, NA_REAL);
  NumericVector sizesc_true (matsize, NA_REAL);
  IntegerVector repstatus_true (matsize, 0);
  IntegerVector obsstatus_true (matsize, 1);
  IntegerVector propstatus_true (matsize, 0);
  IntegerVector matstatus_true (matsize, 0);
  IntegerVector immstatus_true (matsize, 1);
  IntegerVector indataset_true (matsize, 1);
  NumericVector minage_true (matsize, 0.0);
  NumericVector maxage_true (matsize, 0.0);
  NumericVector binhalfwidth_true (matsize, NA_REAL);
  NumericVector binhalfwidthb_true (matsize, NA_REAL);
  NumericVector binhalfwidthc_true (matsize, NA_REAL);
  NumericVector sizebin_min (matsize, NA_REAL);
  NumericVector sizebin_max (matsize, NA_REAL);
  NumericVector sizebin_center (matsize, NA_REAL);
  NumericVector sizebin_width (matsize, NA_REAL);
  NumericVector sizebinb_min (matsize, NA_REAL);
  NumericVector sizebinb_max (matsize, NA_REAL);
  NumericVector sizebinb_center (matsize, NA_REAL);
  NumericVector sizebinb_width (matsize, NA_REAL);
  NumericVector sizebinc_min (matsize, NA_REAL);
  NumericVector sizebinc_max (matsize, NA_REAL);
  NumericVector sizebinc_center (matsize, NA_REAL);
  NumericVector sizebinc_width (matsize, NA_REAL);
  IntegerVector group_true (matsize, 0);
  StringVector comments_true (matsize, "No description");
  
  for (int i = 0; i < total_ages; i++) {
    Rcpp::String part1("Age");
    part1 += (static_cast<char>(i + min_age));
    agenames_true(i) = part1;
    
    if ((i + min_age) >= min_fecage) {
      repstatus_true(i) = 1;
      matstatus_true(i) = 1;
      immstatus_true(i) = 0;
    }
    
    minage_true(i) = static_cast<double>(i + min_age);
    if ((i + min_age < max_age) || !cont) {
      maxage_true(i) = static_cast<double>(i + min_age);
    } else {
      maxage_true(i) = NA_REAL;
    }
  }
  
  output_longlist(0) = agenames_true;
  output_longlist(1) = size_true;
  output_longlist(2) = sizesb_true;
  output_longlist(3) = sizesc_true;
  
  output_longlist(4) = minage_true;
  output_longlist(5) = maxage_true;
  output_longlist(6) = repstatus_true;
  output_longlist(7) = obsstatus_true;
  output_longlist(8) = propstatus_true;
  output_longlist(9) = immstatus_true;
  output_longlist(10) = matstatus_true;
  
  output_longlist(11) = indataset_true;
  
  output_longlist(12) = binhalfwidth_true;
  output_longlist(13) = sizebin_min;
  output_longlist(14) = sizebin_max;
  output_longlist(15) = sizebin_center;
  output_longlist(16) = sizebin_width;
  
  output_longlist(17) = binhalfwidthb_true;
  output_longlist(18) = sizebinb_min;
  output_longlist(19) = sizebinb_max;
  output_longlist(20) = sizebinb_center;
  output_longlist(21) = sizebinb_width;
  
  output_longlist(22) = binhalfwidthc_true;
  output_longlist(23) = sizebinc_min;
  output_longlist(24) = sizebinc_max;
  output_longlist(25) = sizebinc_center;
  output_longlist(26) = sizebinc_width;
  
  output_longlist(27) = group_true;
  output_longlist(28) = comments_true;
  
  output_longlist.attr("names") = varnames;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, matsize);
  StringVector needed_classes {"data.frame", "stageframe"};
  output_longlist.attr("class") = needed_classes; // data.frame
  
  return output_longlist;
}

// Matrix Extimation

//' Compares Two Strings Literally
//' 
//' This function compares two strings element by element. Returns \code{FALSE}
//' in case of any differences whatsoever.
//' 
//' @name stringcompare_hard
//' 
//' @param str1 The first string
//' @param str2 The second string
//' 
//' @return A logical value. In case of any difference at all, it will return
//' \code{FALSE}.
//' 
//' @keywords internal
//' @noRd
bool stringcompare_hard(std::string str1, std::string str2) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  bool same = true;
  
  if (str1_length == str2_length && str1_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (str1[i] != str2[i]) {
        same = false;
      }
    }
  } else if (str1_length != str2_length) {
    same = false;
  }
  
  return same;
}

//' Compares Two Strings, Assessing Inclusion
//' 
//' This function compares two strings, and will assess whether \code{str2} is
//' contained within \code{str1}.
//' 
//' @name stringcompare_soft
//' 
//' @param str1 The first string
//' @param str2 The second string
//' 
//' @return A list of two values. The first is a logical value indicating
//' whether \code{str2} occurs within \code{str1}. The second element is an
//' integer indicating at what element of \code{str1} \code{str2} begins.
//' \code{FALSE}.
//' 
//' @keywords internal
//' @noRd
List stringcompare_soft(std::string str1, std::string str2) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  int rem_check {0};
  bool same = false;
  unsigned int start_index {0};
  
  //int rem_check = str1_length;
  
  if (str1_length >= str2_length && str2_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (str1[i] != str2[rem_check]) {
        rem_check = 0;
        // same = false;
      } else {
        if (rem_check == 0) start_index = i;
        rem_check += 1;
        if (rem_check >= str2_length) break;
      }
    }
    
    if (rem_check == str2_length) {
      same = true;
    }
  }
  
  List output = List::create(_["contains"] = same, _["start_index"] = start_index);
  
  return output;
}

//' Compares Two Strings, Assessing Inclusion
//' 
//' This function compares two strings, and will assess whether \code{str2} is
//' contained within \code{str1}. It is a simpler version of 
//' \code{stringcompare_soft()} that yields only the logical result.
//' 
//' @name stringcompare_simple
//' 
//' @param str1 The first string
//' @param str2 The second string
//' @param lower A logical value indicating whether to change all inputs to
//' lower case before checking.
//' 
//' @return A logical value indicating whether \code{str2} occurs within
//' \code{str1}.
//' 
//' @keywords internal
//' @noRd
bool stringcompare_simple(std::string str1, std::string str2, bool lower = false) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  int rem_check {0};
  bool same = false;
  unsigned int start_index {0};
  
  //int rem_check = str1_length;
  
  if (str1_length >= str2_length && str2_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (!lower) {
        if (str1[i] != str2[rem_check]) {
          rem_check = 0;
          // same = false;
        } else {
          if (rem_check == 0) start_index = i;
          rem_check += 1;
          if (rem_check >= str2_length) break;
        }
      } else {
        if (tolower(str1[i]) != tolower(str2[rem_check])) {
          rem_check = 0;
          // same = false;
        } else {
          if (rem_check == 0) start_index = i;
          rem_check += 1;
          if (rem_check >= str2_length) break;
        }
      }
    }
    
    if (rem_check == str2_length) {
      same = true;
    }
  }
  
  return same;
}

//' Compares Three Strings for Interaction Notation
//' 
//' This function compares checks to see if one string is composed of the other
//' two strings in R's interaction notation.
//' 
//' @name stringcompare_x
//' 
//' @param str1 The first string. Used for comparison.
//' @param str2 The second string. Will be incorporated into interaction format.
//' @param str3 The third string. Will be incorporated into interaction format.
//' 
//' @return A logical value. In case of any difference at all, it will return
//' \code{FALSE}.
//' 
//' @keywords internal
//' @noRd
bool stringcompare_x(std::string str1, std::string str2, std::string str3) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  int str3_length = str3.size();
  int combined_length = str2_length + str3_length + 1;
  bool same = false;
  bool same1 = true;
  bool same2 = true;
  
  if (str1_length == combined_length && str1_length > 0) {
    std::string x1 = str2;
    x1 += ":";
    x1 += str3;
    
    std::string x2 = str3;
    x2 += ":";
    x2 += str2;
    
    for (int i = 0; i < str1_length; i++) {
      if (str1[i] != x1[i]) {
        same1 = false;
      }
      if (str1[i] != x2[i]) {
        same2 = false;
      }
    }
  } else {
    same1 = false;
    same2 = false;
  }
  
  if (same1 || same2) same = true;
  
  return same;
}

//' Sort String Elements
//' 
//' This function is based on code obtained from R Bloggers
//' (see https://www.r-bloggers.com/2013/01/handling-strings-with-rcpp/). It
//' sorts the elements of a string vector in alphabetical order.
//' 
//' @name stringsort
//' 
//' @param string_input A string vector.
//' 
//' @return The sorted string vector.
//' 
//' @keywords internal
//' @noRd
CharacterVector stringsort(CharacterVector string_input ) {
  int len = string_input.size();
  
  std::vector<std::string> converted(len);
  for (int i=0; i < len; i++) converted[i] = as<std::string>(string_input(i));
  std::sort( converted.begin(), converted.end() );
  
  CharacterVector new_converted(len);
  new_converted = converted;
  
  return new_converted;
}

//' Sort Integer Elements
//' 
//' This function is based on code obtained from the Rcpp Gallery by Ross
//' Bennett (see https://gallery.rcpp.org/articles/sorting/). It sorts the
//' elements of an integer vector.
//' 
//' @name int_sort
//' 
//' @param int_input An integer vector.
//' 
//' @return The sorted integer vector.
//' 
//' @keywords internal
//' @noRd
IntegerVector int_sort(IntegerVector x) {
   IntegerVector y = clone(x);
   std::sort(y.begin(), y.end());
   
   return y;
}

//' Function to Index a Numeric Vector According to a Reference Vector
//' 
//' Function \code{refsort_num()} takes a numeric matrix and replaces it with an
//' integer vector showing the position of each element in the input vector
//' within the reference vector.
//' 
//' @name refsort_num
//' 
//' @param vec The matrix to index
//' @param ref The vector to use as a reference
//' 
//' @return An integer vector with integers referring to elements in vector
//' \code{ref}.
//' 
//' @keywords internal
//' @noRd
IntegerMatrix refsort_num(NumericMatrix vec, NumericVector ref) {
  int vec_length = vec.length();
  int ref_length = ref.length();
  
  IntegerMatrix output(vec.nrow(), vec.ncol());
  
  for (int i = 0; i < vec_length; i++) {
    for (int j = 0; j < ref_length; j++) {
      if (vec[i] == ref[j]) output[i] = j + 1;
    }
  }
  
  return output;
}

//' Function to Index a Numeric Vector According to a Reference Vector
//' 
//' Function \code{refsort_str()} takes a string vector and replaces it with an
//' integer vector showing the position of each element in the input vector
//' within the reference vector.
//' 
//' @name refsort_str
//' 
//' @param vec The vector to index
//' @param ref The vector to use as a reference
//' 
//' @return An integer vector with integers referring to elements in vector
//' \code{ref}.
//' 
//' @keywords internal
//' @noRd
IntegerVector refsort_str(CharacterVector vec, CharacterVector ref) {
  int vec_length = vec.length();
  int ref_length = ref.length();
  
  IntegerVector output(vec_length);
  
  for (int i = 0; i < vec_length; i++) {
    for (int j = 0; j < ref_length; j++) {
      if (stringcompare_hard(as<std::string>(vec[i]), as<std::string>(ref[j]))) output[i] = j + 1;
    }
  }
  
  return output;
}

//' Create hstages Index Object
//' 
//' Function \code{hst_maker()} creates \code{hstages} index data frames from
//' \code{stageframe} inputs.
//' 
//' @name hst_maker
//' 
//' @param sframe The ahistorical stageframe used in MPM development.
//' 
//' @return A data frame with the following columns:
//' \item{stage_id_2}{Integer index of stage in time \emph{t}+1.}
//' \item{stage_id_1}{Integer index of stage in time \emph{t}.}
//' \item{stage_2}{String name of stage in time \emph{t}+1.}
//' \item{stage_1}{String name of stage in time \emph{t}.}
//' 
//' @keywords internal
//' @noRd
DataFrame hst_maker (DataFrame sframe) {
  StringVector stage_name = as<StringVector>(sframe["stage"]);
  int true_stages = stage_name.length();
  
  IntegerVector stage_id = seq(1, true_stages);
  int h_stages = true_stages * true_stages;
  
  IntegerVector stage_id_2 (h_stages);
  IntegerVector stage_id_1 (h_stages);
  StringVector stage_2 (h_stages);
  StringVector stage_1 (h_stages);
  
  for (int s1 = 0; s1 < true_stages; s1++) {
    for (int s2 = 0; s2 < true_stages; s2++) {
      int current_elem = (s1 * true_stages) + s2;
      
      stage_id_2[current_elem] = stage_id[s2];
      stage_id_1[current_elem] = stage_id[s1];
      stage_2[current_elem] = stage_name[s2];
      stage_1[current_elem] = stage_name[s1];
    }
  }
  
  DataFrame output = DataFrame::create(_["stage_id_2"] = stage_id_2,
    _["stage_id_1"] = stage_id_1, _["stage_2"] = stage_2, _["stage_1"] = stage_1);
  
  return output;
}

//' Create agestages Index Object
//' 
//' Function \code{age_maker()} creates \code{agestages} index data frames from
//' \code{stageframe} inputs.
//' 
//' @name age_maker
//' 
//' @param sframe The ahistorical stageframe used in MPM development.
//' 
//' @return A data frame with the following columns:
//' \item{stage_id}{Integer index of stage.}
//' \item{stage}{String name of stage.}
//' \item{age}{The age of stage in current time.}
//' 
//' @keywords internal
//' @noRd
DataFrame age_maker (DataFrame sframe, int start_age, int last_age) {
  StringVector stage_name = as<StringVector>(sframe["stage"]);
  int true_stages = stage_name.length();
  
  IntegerVector stage_id = seq(1, true_stages);
  IntegerVector all_ages = seq(start_age, last_age);
  int num_ages = all_ages.length();
  int age_stages = true_stages * num_ages;
  
  IntegerVector stage_id_new (age_stages);
  StringVector stage_new (age_stages);
  IntegerVector age_new (age_stages);
  
  for (int s1 = 0; s1 < num_ages; s1++) {
    for (int s2 = 0; s2 < true_stages; s2++) {
      int current_elem = (s1 * true_stages) + s2;
      
      stage_id_new[current_elem] = stage_id[s2];
      stage_new[current_elem] = stage_name[s2];
      age_new[current_elem] = all_ages[s1];
    }
  }
  
  DataFrame output = DataFrame::create(_["stage_id"] = stage_id_new,
    _["stage"] = stage_new, _["age"] = age_new);
  
  return output;
}

//' Re-index Projection Matrix On Basis of Overwrite Table
//' 
//' Function \code{ovreplace()} takes matrix indices provided by functions
//' \code{\link{rlefko3}()}, \code{\link{rlefko2}()}, \code{\link{flefko3}()},
//' \code{\link{flefko2}()}, and \code{\link{aflefko2}()} and updates them with
//' information provided in the overwrite table used as input in that function.
//' 
//' @name ovreplace
//' 
//' @param allst321 Vector containing the original element-by-element matrix
//' index.
//' @param idx321old Vector containing the indices of matrix elements to be
//' updated.
//' @param idx321new Vector containing the replacement matrix element indices.
//' @param convtype Vector denoting survival transition (1), fecundity (2), or
//' fecundity multiplier (3).
//' @param eststag3 Vector of new stages in time \emph{t}+1.
//' @param gvnrate Vector of replacement transition values.
//' @param multipl Vector of fecundity multipliers.
//' 
//' @return A matrix. Column 1 is the given rate for a survival transitions,
//' Column 2 is the proxy transition to be used to estimate that transition.
//' Column 3 is the given rate for a fecundity transitions. Column 4 is the
//' proxy transition to be used to estimate that transition. Column 5 is a
//' vector of fecundity multipliers, in cases where no given rate or proxy is to
//' be used but fecundity is to be multiplied by some value. Column 6 is a
//' vector of survival transition multipliers. Column 7 is a vector of fecundity
//' transition multipliers.
//' 
//' @keywords internal
//' @noRd
arma::mat ovreplace(arma::vec allst321, arma::vec idx321old,
  arma::vec idx321new, arma::vec convtype, arma::vec eststag3, 
  arma::vec gvnrate, arma::vec multipl) {
  
  int n = idx321new.n_elem;
  
  arma::mat replacements(allst321.n_elem, 7);
  replacements.fill(-1.0);
  
  for (int i = 0; i < n; i++) {
    arma::uvec correctplace = find(allst321 == idx321old[i]);
    
    int m = correctplace.n_elem; 
    
    for (int j = 0; j < m; j++) {
      if (convtype[i] == 1.0) {
        if (gvnrate[i] >= 0) {replacements(correctplace[j], 0) = gvnrate[i];}
        if (eststag3[i] != -1 && idx321new[i] >= 0) {replacements(correctplace[j], 1) = idx321new[i];}
      }
      
      if (convtype[i] == 2.0) {
        if (gvnrate[i] >= 0) {replacements(correctplace[j], 2) = gvnrate[i];}
        if (eststag3[i] != -1 && idx321new[i] >= 0) {replacements(correctplace[j], 3) = idx321new[i];}
      }
      
      if (convtype[i] == 3.0) {
        replacements(correctplace[j], 4) = multipl[i];
      } else if (convtype[i] == 1) {
        replacements(correctplace[j], 5) = multipl[i];
      } else if (convtype[i] == 2) {
        replacements(correctplace[j], 6) = multipl[i];
      }
    }
  }
  
  return replacements;
}

//' Create Element Index for Matrix Estimation
//' 
//' Function \code{theoldpizzle()} creates a data frame object used by 
//' functions \code{\link{specialpatrolgroup}()},
//' \code{\link{normalpatrolgroup}()}, and \code{jerzeibalowski()} to estimate
//' raw and function-derived matrices.
//' 
//' @name theoldpizzle
//'
//' @param StageFrame The stageframe object identifying the life history model
//' being operationalized.
//' @param OverWrite The overwrite table used in analysis, as modified by 
//' \code{.overwrite_reassess}. Must be processed via \code{.overwrite_reassess}
//' rather than being a raw overwrite or supplement table.
//' @param repmatrix The reproductive matrix used in analysis.
//' @param firstage The first age to be used in the analysis. Should typically
//' be \code{0} for pre-breeding and \code{1} for post-breeding life history
//' models. If not building age-by-stage MPMs, then should be set to \code{0}.
//' @param finalage The final age to be used in analysis. If not building
//' age-by-stage MPMs, then should be set to \code{0}.
//' @param format Indicates whether historical matrices should be in (1) Ehrlen
//' or (2) deVries format.
//' @param style The style of analysis, where 0 is historical, 1 is ahistorical,
//' and 2 is age-by-stage.
//' @param cont Denotes whether age-by-stage matrix continues past the final
//' age.
//' @param filter An integer denoting whether to filter the DataFrame to
//' eliminate unusable rows, and if so, how to do so. Possible values: \code{0}:
//' no filtering, \code{1}: filter out rows with \code{index321 == -1}, and
//' \code{2}: filter out rows with \code{aliveandequal == -1}.
//' 
//' @return The output is a large data frame describing every element to be
//' estimated in matrices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.theoldpizzle)]]
Rcpp::List theoldpizzle(DataFrame StageFrame, DataFrame OverWrite,
  arma::mat repmatrix, int firstage, int finalage, int format, int style,
  int cont, int filter) {
  
  StringVector ovstage3 = as<StringVector>(OverWrite["stage3"]);
  StringVector ovstage2 = as<StringVector>(OverWrite["stage2"]);
  StringVector ovstage1 = as<StringVector>(OverWrite["stage1"]);
  StringVector oveststage3 = as<StringVector>(OverWrite["eststage3"]);
  StringVector oveststage2 = as<StringVector>(OverWrite["eststage2"]);
  StringVector oveststage1 = as<StringVector>(OverWrite["eststage1"]);
  arma::vec ovgivenrate = as<arma::vec>(OverWrite["givenrate"]);
  arma::vec ovmultiplier = as<arma::vec>(OverWrite["multiplier"]);
  arma::vec ovconvtype = as<arma::vec>(OverWrite["convtype"]);
  arma::vec ovconvt12 = as<arma::vec>(OverWrite["convtype_t12"]);
  int ovrows = ovconvtype.n_elem;
  
  int totalages = (finalage - firstage) + 1;
  
  arma::vec ovindex3(ovrows * totalages);
  arma::vec ovindex2(ovrows * totalages);
  arma::vec ovindex1(ovrows * totalages);
  arma::vec ovnew3(ovrows * totalages);
  arma::vec ovnew2(ovrows * totalages);
  arma::vec ovnew1(ovrows * totalages);
  arma::vec ovindexold321(ovrows * totalages);
  arma::vec ovindexnew321(ovrows * totalages);
  arma::vec ovnewgivenrate(ovrows * totalages);
  arma::vec ovnewmultiplier(ovrows * totalages);
  arma::vec ovconvtypeage(ovrows * totalages);
  ovindex3.fill(-1.0);
  ovindex2.fill(-1.0);
  ovindex1.fill(-1.0);
  ovnew3.fill(-1.0);
  ovnew2.fill(-1.0);
  ovnew1.fill(-1.0);
  ovindexold321.fill(-1.0);
  ovindexnew321.fill(-1.0);
  ovnewgivenrate.fill(-1.0);
  ovnewmultiplier.zeros();
  ovconvtypeage.fill(-1.0);
  
  arma::vec newstageid = as<arma::vec>(StageFrame["stage_id"]);
  StringVector origstageid = as<StringVector>(StageFrame["stage"]);
  arma::vec binsizectr = as<arma::vec>(StageFrame["sizebin_center"]);
  arma::vec repstatus = as<arma::vec>(StageFrame["repstatus"]);
  arma::vec obsstatus = as<arma::vec>(StageFrame["obsstatus"]);
  arma::vec immstatus = as<arma::vec>(StageFrame["immstatus"]);
  arma::vec matstatus = as<arma::vec>(StageFrame["matstatus"]);
  arma::vec indata = as<arma::vec>(StageFrame["indataset"]);
  arma::vec binsizewidth = as<arma::vec>(StageFrame["sizebin_width"]);
  arma::vec alive = as<arma::vec>(StageFrame["alive"]);
  arma::vec minage = as<arma::vec>(StageFrame["min_age"]);
  arma::vec maxage = as<arma::vec>(StageFrame["max_age"]);
  arma::vec group = as<arma::vec>(StageFrame["group"]);
  arma::vec almostborn = as<arma::vec>(StageFrame["almostborn"]);
  
  arma::vec binsizebctr = as<arma::vec>(StageFrame["sizebinb_center"]);
  arma::vec binsizecctr = as<arma::vec>(StageFrame["sizebinc_center"]);
  arma::vec binsizebwidth = as<arma::vec>(StageFrame["sizebinb_width"]);
  arma::vec binsizecwidth = as<arma::vec>(StageFrame["sizebinc_width"]);
  
  // This section determines the length of the matrix map data frame
  int nostages = newstageid.n_elem;
  int nostages_nodead = nostages - 1;
  int nostages_nounborn = nostages;
  int nostages_nodead_nounborn = nostages_nodead;
  int prior_stage = -1;
  arma::vec ovrepentry_prior(nostages, fill::zeros);
  
  int totallength {0};
  
  if (style == 2) {
    totallength = (nostages * nostages * totalages * totalages);
  } else if (style == 1) {
    totallength = (nostages * nostages_nodead);
  } else {
    if (format == 2) {
      nostages_nodead_nounborn = nostages - 2;
      prior_stage = nostages_nodead_nounborn;
      nostages_nounborn = nostages - 1;
      totallength = (2 * nostages_nodead_nounborn * nostages_nounborn * nostages_nounborn);
    } else {
      totallength = (nostages * (nostages_nodead * nostages_nodead));
    }
  }
  
  // This section sets up the repmatrix. First we will determine whether the
  // repmatrix has been entered in historical or ahistorical format, since this
  // does not necessarily match the MPM type
  int reprows = repmatrix.n_rows;
  int repmattype = 0;
  
  if (reprows == (nostages - 1) || reprows == (nostages - 2)) {
    repmattype = 1; // The repmatrix is ahistorical in dimensions
  } else if (reprows == ((nostages - 1) * (nostages - 1)) || 
      reprows == ((nostages - 2) * (nostages - 2))) {
    repmattype = 2; // The repmatrix is historical in dimensions
  }
  
  // Here we set up the vectors that will be put together into the matrix map data frame
  arma::vec stage3(totallength, fill::zeros);
  arma::vec stage2n(totallength, fill::zeros);
  arma::vec stage2o(totallength, fill::zeros);
  arma::vec stage1(totallength, fill::zeros);
  
  arma::vec size3(totallength, fill::zeros);
  arma::vec size2n(totallength, fill::zeros);
  arma::vec size2o(totallength, fill::zeros);
  arma::vec size1(totallength, fill::zeros);
  
  arma::vec sizeb3(totallength, fill::zeros);
  arma::vec sizeb2n(totallength, fill::zeros);
  arma::vec sizeb2o(totallength, fill::zeros);
  arma::vec sizeb1(totallength, fill::zeros);
  
  arma::vec sizec3(totallength, fill::zeros);
  arma::vec sizec2n(totallength, fill::zeros);
  arma::vec sizec2o(totallength, fill::zeros);
  arma::vec sizec1(totallength, fill::zeros);
  
  arma::vec obs3(totallength, fill::zeros);
  arma::vec obs2n(totallength, fill::zeros);
  arma::vec obs2o(totallength, fill::zeros);
  arma::vec obs1(totallength, fill::zeros);
  
  arma::vec rep3(totallength, fill::zeros);
  arma::vec rep2n(totallength, fill::zeros);
  arma::vec rep2o(totallength, fill::zeros);
  arma::vec rep1(totallength, fill::zeros);
  
  arma::vec mat3(totallength, fill::zeros);
  arma::vec mat2n(totallength, fill::zeros);
  arma::vec mat2o(totallength, fill::zeros);
  arma::vec mat1(totallength, fill::zeros);
  
  arma::vec imm3(totallength, fill::zeros);
  arma::vec imm2n(totallength, fill::zeros);
  arma::vec imm2o(totallength, fill::zeros);
  arma::vec imm1(totallength, fill::zeros);
  
  arma::vec repentry3(totallength, fill::zeros);
  arma::vec repentry2o(totallength, fill::zeros);
  arma::vec almostborn1(totallength, fill::zeros);
  
  arma::vec binwidth(totallength, fill::zeros);
  arma::vec binbwidth(totallength, fill::zeros);
  arma::vec bincwidth(totallength, fill::zeros);
  
  arma::vec indata3(totallength, fill::zeros);
  arma::vec indata2n(totallength, fill::zeros);
  arma::vec indata2o(totallength, fill::zeros);
  arma::vec indata1(totallength, fill::zeros);
  
  arma::vec minage3(totallength, fill::zeros);
  arma::vec minage2(totallength, fill::zeros);
  arma::vec maxage3(totallength, fill::zeros);
  arma::vec maxage2(totallength, fill::zeros);
  
  arma::vec grp3(totallength, fill::zeros);
  arma::vec grp2n(totallength, fill::zeros);
  arma::vec grp2o(totallength, fill::zeros);
  arma::vec grp1(totallength, fill::zeros);
  
  arma::vec actualage(totallength, fill::zeros);
  arma::vec index321(totallength);
  arma::vec index21(totallength);
  arma::vec indatalong(totallength, fill::zeros);
  arma::vec aliveequal(totallength);
  arma::vec included(totallength, fill::zeros);
  index321.fill(-1.0);
  index21.fill(-1.0);
  aliveequal.fill(-1.0);
  
  arma::mat asadditions(totallength, 5, fill::zeros);
  
  arma::vec ovgivent(totallength);
  arma::vec ovestt(totallength);
  arma::vec ovgivenf(totallength);
  arma::vec ovestf(totallength);
  arma::vec ovrepentry(totallength, fill::zeros);
  arma::vec ovsurvmult(totallength, fill::ones);
  arma::vec ovfecmult(totallength, fill::ones);
  ovgivent.fill(-1.0);
  ovestt.fill(-1.0);
  ovgivenf.fill(-1.0);
  ovestf.fill(-1.0);
  
  int repm_elem {-1};
  double deadandnasty {0};
  long long int currentindex {0};
  
  // This step changes the stage names to stage numbers per the input stageframe for styles 0 and 1
  if (style < 2) {
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
        for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
          if (ovstage3(i) == origstageid(j)) {
            ovindex3(i) = newstageid(j);
          }
          
          if (ovstage2(i) == origstageid(j)) {
            ovindex2(i) = newstageid(j);
          }
          
          if (ovstage1(i) == origstageid(j)) {
            ovindex1(i) = newstageid(j);
          }
          
          if (oveststage3(i) == origstageid(j)) {
            ovnew3(i) = newstageid(j);
          }
          
          if (oveststage2(i) == origstageid(j)) {
            ovnew2(i) = newstageid(j);
          }
          
          if (oveststage1(i) == origstageid(j)) {
            ovnew1(i) = newstageid(j);
          }
        } // j for loop
      } // i for loop
    } // ovrows if statement
  } // style if statement
  
  // Now we cover the main data frame creation loops
  // When style = 0, this will create AllStages for the historical case
  // When format = 2, the historical MPM will be in deVries format
  if (style == 0 && format == 2) {
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows\
        
        // Here are the new versions
        if (ovconvtype(i) > 1.0) { // Catches all changes to fecundity and reproductive multipliers
          ovindexold321(i) = (ovindex3(i) - 1) + (prior_stage * nostages) + 
            ((ovindex2(i) - 1) * nostages * nostages) + 
            ((ovindex1(i) - 1) * nostages * nostages * nostages);
            
          ovindexnew321(i) = (ovnew3(i) - 1) + (prior_stage * nostages) + 
            ((ovnew2(i) - 1) * nostages * nostages) + 
            ((ovnew1(i) - 1) * nostages * nostages * nostages);
        } else if (ovconvt12(i) == 2.0) { // Catches all survival terms with historical reproduction events
          ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages) + 
            ((ovindex2(i) - 1) * nostages * nostages) + 
            (prior_stage * nostages * nostages * nostages);
            
          ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages) + 
            ((ovnew2(i) - 1) * nostages * nostages) + 
            (prior_stage * nostages * nostages * nostages);
        } else { // This refers to full survival transitions
          ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages) + 
            ((ovindex2(i) - 1) * nostages * nostages) + 
            ((ovindex1(i) - 1) * nostages * nostages * nostages);
            
          ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages) + 
            ((ovnew2(i) - 1) * nostages * nostages) + 
            ((ovnew1(i) - 1) * nostages * nostages * nostages);
        }
        if (ovindexold321(i) < 0.0) ovindexold321(i) = -1.0;
        if (ovindexnew321(i) < 0.0) ovindexnew321(i) = -1.0;
        
        if (!NumericVector::is_na(ovgivenrate(i))) {
          ovnewgivenrate(i) = ovgivenrate(i);
        }
        if (NumericVector::is_na(ovmultiplier(i))) {
          ovmultiplier(i) = 1;
        }
        ovnewmultiplier(i) = ovmultiplier(i);
        
        if (ovconvtype(i) == 3.0) {
          for (int j = 0; j < nostages; j++) {
            if (origstageid(j) == ovstage3(i)) ovrepentry_prior(j) = 1.0;
          }
        }
      } // i for loop
    } // ovrows if statement
    
    arma::uvec marked_for_repentry (nostages, fill::zeros); // Only used in deVries format
    
    for (int time1 = 0; time1 < nostages_nodead; time1++) {
      for (int time2o = 0; time2o < nostages_nodead_nounborn; time2o++) {
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            if (time3 != prior_stage) {
              if (time2n == time2o || time2n == prior_stage){
                
                included(currentindex) = 1.0;
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2o);
                stage1(currentindex) = newstageid(time1);
                
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2o);
                size1(currentindex) = binsizectr(time1);
                
                sizeb3(currentindex) = binsizebctr(time3);
                sizeb2n(currentindex) = binsizebctr(time2n);
                sizeb2o(currentindex) = binsizebctr(time2o);
                sizeb1(currentindex) = binsizebctr(time1);
                
                if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0.0;
                
                sizec3(currentindex) = binsizecctr(time3);
                sizec2n(currentindex) = binsizecctr(time2n);
                sizec2o(currentindex) = binsizecctr(time2o);
                sizec1(currentindex) = binsizecctr(time1);
                
                if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
                if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
                if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0.0;
                
                obs3(currentindex) = obsstatus(time3);
                obs2n(currentindex) = obsstatus(time2n);
                obs2o(currentindex) = obsstatus(time2o);
                obs1(currentindex) = obsstatus(time1);
                
                rep3(currentindex) = repstatus(time3);
                rep2n(currentindex) = repstatus(time2n);
                rep2o(currentindex) = repstatus(time2o);
                rep1(currentindex) = repstatus(time1);
                
                mat3(currentindex) = matstatus(time3);
                mat2n(currentindex) = matstatus(time2n);
                mat2o(currentindex) = matstatus(time2o);
                mat1(currentindex) = matstatus(time1);
                
                imm3(currentindex) = immstatus(time3);
                imm2n(currentindex) = immstatus(time2n);
                imm2o(currentindex) = immstatus(time2o);
                imm1(currentindex) = immstatus(time1);
                
                // This statement fills in the repentry info from the repmatrix
                if (time2n == prior_stage && time3 < prior_stage && time2o < prior_stage) {
                  if (repmattype == 1) { 
                    repm_elem = time3 + (time2o * nostages_nodead_nounborn);
                  } else if (repmattype == 2) {
                    repm_elem = time3 + (time2o * nostages_nodead_nounborn) + 
                      (time2o * nostages_nodead_nounborn * nostages_nodead_nounborn) +
                      (time1 * nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn);
                  } else repm_elem = -1;
                  
                  if (repmatrix(repm_elem) > 0.0) {
                    repentry3(currentindex) = repmatrix(repm_elem);
                    if (repentry3(currentindex) == 0.0 && ovrepentry_prior(time3) != 0.0) {
                      repentry3(currentindex) = 1.0;
                      marked_for_repentry(stage3(currentindex)) = 1;
                    } 
                  }
                } else repentry3(currentindex) = 0.0;
                
                almostborn1(currentindex) = almostborn(time1);
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2o);
                indata1(currentindex) = indata(time1);
                
                binwidth(currentindex) = binsizewidth(time3);
                binbwidth(currentindex) = binsizebwidth(time3);
                bincwidth(currentindex) = binsizecwidth(time3);
                
                if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
                if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2o);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2o);
                actualage(currentindex) = 0.0;
                
                grp3(currentindex) = group(time3);
                grp2n(currentindex) = group(time2n);
                grp2o(currentindex) = group(time2o);
                grp1(currentindex) = group(time1);
                
                if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
                  deadandnasty = 1.0;
                } else if (stage2o(currentindex) == nostages || stage1(currentindex) == nostages) {
                  deadandnasty = 1.0;
                } else {
                  deadandnasty = 0.0;
                }
                
                if (deadandnasty == 0.0) {
                  // The next index variable gives the element in the final matrix
                  aliveequal(currentindex) = (stage3(currentindex) - 1) + 
                    ((stage2n(currentindex) - 1) * nostages_nodead_nounborn) + 
                    ((stage2o(currentindex) - 1) * nostages_nodead * nostages_nodead_nounborn) + 
                    ((stage1(currentindex) - 1) * nostages_nodead_nounborn * 
                      nostages_nodead * nostages_nodead_nounborn);
                  
                  // The next two index variables are used by ovreplace
                  index321(currentindex) = (stage3(currentindex) - 1) + 
                    ((stage2n(currentindex) - 1) * nostages) + 
                    ((stage2o(currentindex) - 1) * nostages * nostages) + 
                    ((stage1(currentindex) - 1) * nostages * nostages * nostages);
                    
                  index21(currentindex) = (stage2o(currentindex) - 1) + 
                    ((stage1(currentindex) - 1) * nostages);
                }
                
                indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
                  indata2o(currentindex) * indata1(currentindex);
                
                currentindex += 1;
              } // if (time2n == tim2o || time2n == prior_stage) statement
            } // if (time3n != dead_stage) statement
          } // time3 loop
        } // time2n loop
      } // time2o loop
    } // time1 loop 
    
    // Here we edit the data frame to make sure that almostborn situations in time 1
    // lead to estimated elements only if a repentry stage occurs in time 2
    arma::uvec marked_only = find(marked_for_repentry);
    if (marked_only.n_elem > 0) {
      for (int i = 0; i < marked_only.n_elem; i++) {
        arma::uvec total_indices_to_change = find(stage2o == marked_only(i));
        
        if (total_indices_to_change.n_elem > 0) {
          for (int j = 0; j < total_indices_to_change.n_elem; j++) {
            repentry2o(total_indices_to_change(j)) = 1;
          }
        }
      }
    }
    
    arma::uvec alm_only = find(almostborn1);
    if (alm_only.n_elem > 0) {
      for (int i = 0; i < alm_only.n_elem; i++) {
        if (repentry2o(alm_only(i)) < 1.0) {
          index321(alm_only(i)) = -1.0;
        }
      }
    }
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtype, 
        ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      
      ovrepentry = asadditions.col(4);
      ovsurvmult = asadditions.col(5);
      ovfecmult = asadditions.col(6);
      
      arma::uvec workedupindex = find(ovrepentry > 0.0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } else if (style == 0 && format == 1) { // Historical MPM in Ehrlen format
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
        
        ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages_nodead_nounborn) + 
          ((ovindex2(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
          ((ovindex1(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
            nostages_nodead_nounborn);
          
        ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages_nodead) + 
          ((ovnew2(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
          ((ovnew1(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
            nostages_nodead_nounborn);
        
        if (ovindexold321(i) < 0) ovindexold321(i) = -1.0;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1.0;
        
        if (!NumericVector::is_na(ovgivenrate(i))) {
          ovnewgivenrate(i) = ovgivenrate(i);
        }
        if (NumericVector::is_na(ovmultiplier(i))) {
          ovmultiplier(i) = 1.0;
        }
        ovnewmultiplier(i) = ovmultiplier(i);
      } // i for loop
    } // ovrows if statement
    
    for (int time1 = 0; time1 < nostages_nodead; time1++) {
      for (int time2o = 0; time2o < nostages_nodead; time2o++) {
        for (int time3 = 0; time3 < nostages; time3++) {
          
          included(currentindex) = 1.0;
          
          stage3(currentindex) = newstageid(time3);
          stage2n(currentindex) = newstageid(time2o);
          stage2o(currentindex) = newstageid(time2o);
          stage1(currentindex) = newstageid(time1);
          
          size3(currentindex) = binsizectr(time3);
          size2n(currentindex) = binsizectr(time2o);
          size2o(currentindex) = binsizectr(time2o);
          size1(currentindex) = binsizectr(time1);
          
          sizeb3(currentindex) = binsizebctr(time3);
          sizeb2n(currentindex) = binsizebctr(time2o);
          sizeb2o(currentindex) = binsizebctr(time2o);
          sizeb1(currentindex) = binsizebctr(time1);
          
          if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
          if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
          if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
          if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0.0;
                
          sizec3(currentindex) = binsizecctr(time3);
          sizec2n(currentindex) = binsizecctr(time2o);
          sizec2o(currentindex) = binsizecctr(time2o);
          sizec1(currentindex) = binsizecctr(time1);
          
          if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
          if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
          if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
          if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0.0;
          
          obs3(currentindex) = obsstatus(time3);
          obs2n(currentindex) = obsstatus(time2o);
          obs2o(currentindex) = obsstatus(time2o);
          obs1(currentindex) = obsstatus(time1);
          
          rep3(currentindex) = repstatus(time3);
          rep2n(currentindex) = repstatus(time2o);
          rep2o(currentindex) = repstatus(time2o);
          rep1(currentindex) = repstatus(time1);
          
          mat3(currentindex) = matstatus(time3);
          mat2n(currentindex) = matstatus(time2o);
          mat2o(currentindex) = matstatus(time2o);
          mat1(currentindex) = matstatus(time1);
          
          imm3(currentindex) = immstatus(time3);
          imm2n(currentindex) = immstatus(time2o);
          imm2o(currentindex) = immstatus(time2o);
          imm1(currentindex) = immstatus(time1);
          
          //This next section determines repentry3 on the basis of the input repmatrix
          if (time3 < nostages_nodead_nounborn) {
            if (repmattype == 1) {
              repm_elem = time3 + (time2o * nostages_nodead_nounborn);
            } else if (repmattype == 2) {
              repm_elem = time3 + (time2o * nostages_nodead_nounborn) + 
                (time2o * nostages_nodead_nounborn * nostages_nodead_nounborn) +
                (time1 * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                  nostages_nodead_nounborn);
            } else {
              repm_elem = -1;
            }
          }
          
          if(repm_elem > -1) {
            if (repmatrix(repm_elem) > 0.0) {
              repentry3(currentindex) = repmatrix(repm_elem);
            }
          }
          
          if (time3 < nostages_nodead_nounborn) {
            if (repmattype == 1) { // Ahistorical repmatrix
              repentry3(currentindex) = repmatrix((time3 + (nostages_nodead_nounborn * time2o)));
            } else if (repmattype == 2) {  // Historical repmatrix
              repentry3(currentindex) = repmatrix((time3 + (nostages_nodead_nounborn * time2o)) + 
                ((nostages_nodead_nounborn * nostages_nodead_nounborn * time2o)) +
                (nostages_nodead_nounborn * nostages_nodead_nounborn * 
                  nostages_nodead_nounborn * time1));
            }
          } else {
            repentry3(currentindex) = 0.0;
          }
          
          indata3(currentindex) = indata(time3);
          indata2n(currentindex) = indata(time2o);
          indata2o(currentindex) = indata(time2o);
          indata1(currentindex) = indata(time1);
          
          binwidth(currentindex) = binsizewidth(time3);
          binbwidth(currentindex) = binsizebwidth(time3);
          bincwidth(currentindex) = binsizecwidth(time3);
          
          if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
          if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
          
          minage3(currentindex) = minage(time3);
          minage2(currentindex) = minage(time2o);
          maxage3(currentindex) = maxage(time3);
          maxage2(currentindex) = maxage(time2o);
          actualage(currentindex) = 0.0;
          
          grp3(currentindex) = group(time3);
          grp2n(currentindex) = group(time2o);
          grp2o(currentindex) = group(time2o);
          grp1(currentindex) = group(time1);
          
          if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
            deadandnasty = 1.0;
          } else if (stage2o(currentindex) == nostages || stage1(currentindex) == nostages) {
            deadandnasty = 1.0;
          } else {
            deadandnasty = 0.0;
          }
          
          if (deadandnasty == 0.0) {
            aliveequal(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * 
                (nostages - 1)) + ((stage2o(currentindex) - 1) * (nostages - 1) * (nostages - 1)) + 
              ((stage1(currentindex) - 1) * (nostages - 1) * (nostages - 1) * (nostages - 1));
            
            index321(currentindex) = (stage3(currentindex) - 1) + 
              ((stage2n(currentindex) - 1) * nostages_nodead_nounborn) + 
              ((stage2n(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
              ((stage1(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                nostages_nodead_nounborn);
            index21(currentindex) = (stage2n(currentindex) - 1) + ((stage1(currentindex) - 1) * nostages);
          }
          
          indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
            indata2o(currentindex) * indata1(currentindex);
          
          currentindex += 1;
        } // time3 loop
      } // time2o loop
    } // time1 loop 
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtype, 
        ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      ovsurvmult = asadditions.col(5);
      ovfecmult = asadditions.col(6);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0.0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } else if (style == 1) { // Takes care of the ahistorical case
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
      
        ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages);
        ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages);
        
        if (ovindexold321(i) < 0) ovindexold321(i) = -1.0;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1.0;
        
        if (!NumericVector::is_na(ovgivenrate(i))) {
          ovnewgivenrate(i) = ovgivenrate(i);
        }
        if (NumericVector::is_na(ovmultiplier(i))) {
          ovmultiplier(i) = 1;
        }
        ovnewmultiplier(i) = ovmultiplier(i);
      } // i for loop
    } // ovrows if statement
    
    for (int time2n = 0; time2n < nostages_nodead; time2n++) {
      for (int time3 = 0; time3 < nostages; time3++) {
        
        stage3(currentindex) = newstageid(time3);
        stage2n(currentindex) = newstageid(time2n);
        stage2o(currentindex) = newstageid(time2n);
        stage1(currentindex) = 0.0;
        
        size3(currentindex) = binsizectr(time3);
        size2n(currentindex) = binsizectr(time2n);
        size2o(currentindex) = binsizectr(time2n);
        size1(currentindex) = 0.0;
        
        sizeb3(currentindex) = binsizebctr(time3);
        sizeb2n(currentindex) = binsizebctr(time2n);
        sizeb2o(currentindex) = binsizebctr(time2n);
        sizeb1(currentindex) = 0.0;
        
        if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
        if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
        if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
        
        sizec3(currentindex) = binsizecctr(time3);
        sizec2n(currentindex) = binsizecctr(time2n);
        sizec2o(currentindex) = binsizecctr(time2n);
        sizec1(currentindex) = 0.0;
        
        if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
        if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
        if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
        
        obs3(currentindex) = obsstatus(time3);
        obs2n(currentindex) = obsstatus(time2n);
        obs2o(currentindex) = obsstatus(time2n);
        obs1(currentindex) = 0.0;
        
        rep3(currentindex) = repstatus(time3);
        rep2n(currentindex) = repstatus(time2n);
        rep2o(currentindex) = repstatus(time2n);
        rep1(currentindex) = 0.0;
        
        mat3(currentindex) = matstatus(time3);
        mat2n(currentindex) = matstatus(time2n);
        mat2o(currentindex) = matstatus(time2n);
        mat1(currentindex) = 0.0;
        
        imm3(currentindex) = immstatus(time3);
        imm2n(currentindex) = immstatus(time2n);
        imm2o(currentindex) = immstatus(time2n);
        imm1(currentindex) = 0.0;
        
        if (time3 < nostages_nodead) {
          repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
        } else {
          repentry3(currentindex) = 0.0;
        }
        
        indata3(currentindex) = indata(time3);
        indata2n(currentindex) = indata(time2n);
        indata2o(currentindex) = indata(time2n);
        indata1(currentindex) = 1.0;
        
        binwidth(currentindex) = binsizewidth(time3);
        binbwidth(currentindex) = binsizebwidth(time3);
        bincwidth(currentindex) = binsizecwidth(time3);
        
        if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
        if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
        
        minage3(currentindex) = minage(time3);
        minage2(currentindex) = minage(time2n);
        maxage3(currentindex) = maxage(time3);
        maxage2(currentindex) = maxage(time2n);
        actualage(currentindex) = 0.0;
        
        grp3(currentindex) = group(time3);
        grp2n(currentindex) = group(time2n);
        grp2o(currentindex) = group(time2n);
        grp1(currentindex) = 0.0;
        
        if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
          deadandnasty = 1.0;
        } else {
          deadandnasty = 0.0;
        }
        
        if (deadandnasty == 0.0) {
          aliveequal(currentindex) = (stage3(currentindex) - 1) + 
            ((stage2n(currentindex) - 1) * nostages_nodead);

          index321(currentindex) = (stage3(currentindex) - 1) + 
            ((stage2n(currentindex) - 1) * nostages);
          index21(currentindex) = (stage2n(currentindex) - 1);
        }
        
        indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
          indata2o(currentindex);
          
        currentindex += 1;
        
      } // time3 loop
    } // time2n loop
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtype,
        ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      ovsurvmult = asadditions.col(5);
      ovfecmult = asadditions.col(6);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0.0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } else if (style == 2) { // Takes care of the age x stage case
    int age3 {firstage};
    
    for (int time3 = 0; time3 < nostages; time3++) {
      if (NumericVector::is_na(maxage(time3))) {
        maxage(time3) = finalage + cont;
      }
    }
    
    // This sets up the overwrite tables
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      // This first set of loops establishes a number of indices
      for (int age2 = firstage; age2 < (totalages + 1); age2++) {
        for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
          for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
            ovconvtypeage(i + (ovrows * (age2 - firstage))) = ovconvtype(i);
              
            if (age2 < totalages) {
              if (ovconvtype(i) == 1.0) {
                age3 = age2 + 1;
              } else {
                age3 = firstage;
              }
              
              if (ovstage3(i) == origstageid(j)) {
                ovindex3(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (ovstage2(i) == origstageid(j)) {
                ovindex2(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (oveststage3(i) == origstageid(j)) {
                ovnew3(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (oveststage2(i) == origstageid(j)) {
                ovnew2(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (ovindex3(i + (ovrows * (age2 - firstage))) != -1.0 && 
                ovindex2(i + (ovrows * (age2 - firstage))) != -1.0) {
                ovindexold321(i + (ovrows * (age2 - firstage))) = 
                  ovindex3(i + (ovrows * (age2 - firstage))) +
                  ((age3 - firstage) * nostages) +
                  (ovindex2(i + (ovrows * (age2 - firstage))) * nostages * totalages) + 
                  ((age2 - firstage) * nostages * nostages * totalages);
              }
              
              if (ovnew3(i + (ovrows * (age2 - firstage))) != -1.0 &&
                ovnew2(i + (ovrows * (age2 - firstage))) != -1.0) {
                ovindexnew321(i + (ovrows * (age2 - firstage))) =
                  ovnew3(i + (ovrows * (age2 - firstage))) +
                  ((age3 - firstage) * nostages) +
                  (ovnew2(i + (ovrows * (age2 - firstage))) * nostages * totalages) +
                  ((age2 - firstage) * nostages * nostages * totalages);
              }
              
              if (!NumericVector::is_na(ovgivenrate(i))) {
                ovnewgivenrate(i + (ovrows * (age2 - firstage))) = ovgivenrate(i);
              }
              if (NumericVector::is_na(ovmultiplier(i))) {
                ovmultiplier(i) = 1.0;
              }
              ovnewmultiplier(i + (ovrows * (age2 - firstage))) = ovmultiplier(i);
            } else {
              if (ovconvtype(i) == 1.0) {
                age3 = age2;
              } else {
                age3 = firstage;
              }
              
              if (ovstage3(i) == origstageid(j)) {
                ovindex3(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (ovstage2(i) == origstageid(j)) {
                ovindex2(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (oveststage3(i) == origstageid(j)) {
                ovnew3(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (oveststage2(i) == origstageid(j)) {
                ovnew2(i + (ovrows * (age2 - firstage))) = newstageid(j) - 1.0;
              }
              
              if (ovindex3(i + (ovrows * (age2 - firstage))) != -1.0 &&
                ovindex2(i + (ovrows * (age2 - firstage))) != -1.0) {
                ovindexold321(i + (ovrows * (age2 - firstage))) =
                  ovindex3(i + (ovrows * (age2 - firstage))) +
                  ((age3 - firstage) * nostages) +
                  (ovindex2(i + (ovrows * (age2 - firstage))) * nostages * totalages) +
                  ((age2 - firstage) * nostages * nostages * totalages);
              }
              
              if (ovnew3(i + (ovrows * (age2 - firstage))) != -1.0 &&
                ovnew2(i + (ovrows * (age2 - firstage))) != -1.0) {
                ovindexnew321(i + (ovrows * (age2 - firstage))) =
                  ovnew3(i + (ovrows * (age2 - firstage))) +
                  ((age3 - firstage) * nostages) +
                  (ovnew2(i + (ovrows * (age2 - firstage))) * nostages * totalages) +
                  ((age2 - firstage) * nostages * nostages * totalages);
              }
              if (!NumericVector::is_na(ovgivenrate(i))) {
                ovnewgivenrate(i + (ovrows * (age2 - firstage))) = ovgivenrate(i);
              }
              if (NumericVector::is_na(ovmultiplier(i))) {
                ovmultiplier(i) = 1.0;
              }
              ovnewmultiplier(i + (ovrows * (age2 - firstage))) = ovmultiplier(i);
            }
          } // j for loop
          
        if (ovindexold321(i) < 0) ovindexold321(i) = -1.0;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1.0;
          
        } // i for loop
      } // age loop
    } // ovrows if statement
    
    for (int age2 = firstage; age2 <= finalage; age2++) {
      if (age2 < finalage) { // This first loop takes care of transitions from one age to the next
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            // First survival
            age3 = age2 + 1;
            currentindex = time3 + ((age3 - firstage) * nostages) + 
              (time2n * nostages * totalages) +
              ((age2 - firstage) * nostages * nostages * totalages);
            
            stage3(currentindex) = newstageid(time3);
            stage2n(currentindex) = newstageid(time2n);
            stage2o(currentindex) = newstageid(time2n);
            stage1(currentindex) = 0.0;
            
            size3(currentindex) = binsizectr(time3);
            size2n(currentindex) = binsizectr(time2n);
            size2o(currentindex) = binsizectr(time2n);
            size1(currentindex) = 0.0;
            
            sizeb3(currentindex) = binsizebctr(time3);
            sizeb2n(currentindex) = binsizebctr(time2n);
            sizeb2o(currentindex) = binsizebctr(time2n);
            sizeb1(currentindex) = 0.0;
            
            if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
            
            sizec3(currentindex) = binsizecctr(time3);
            sizec2n(currentindex) = binsizecctr(time2n);
            sizec2o(currentindex) = binsizecctr(time2n);
            sizec1(currentindex) = 0.0;
            
            if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
            if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
            if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
            
            obs3(currentindex) = obsstatus(time3);
            obs2n(currentindex) = obsstatus(time2n);
            obs2o(currentindex) = obsstatus(time2n);
            obs1(currentindex) = 0.0;
            
            rep3(currentindex) = repstatus(time3);
            rep2n(currentindex) = repstatus(time2n);
            rep2o(currentindex) = repstatus(time2n);
            rep1(currentindex) = 0.0;
            
            mat3(currentindex) = matstatus(time3);
            mat2n(currentindex) = matstatus(time2n);
            mat2o(currentindex) = matstatus(time2n);
            mat1(currentindex) = 0.0;
            
            imm3(currentindex) = immstatus(time3);
            imm2n(currentindex) = immstatus(time2n);
            imm2o(currentindex) = immstatus(time2n);
            imm1(currentindex) = 0.0;
            
            repentry3(currentindex) = 0.0;
            
            indata3(currentindex) = indata(time3);
            indata2n(currentindex) = indata(time2n);
            indata2o(currentindex) = indata(time2n);
            indata1(currentindex) = 0.0;
            
            binwidth(currentindex) = binsizewidth(time3);
            binbwidth(currentindex) = binsizebwidth(time3);
            bincwidth(currentindex) = binsizecwidth(time3);
            
            if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
            if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
            
            minage3(currentindex) = minage(time3);
            minage2(currentindex) = minage(time2n);
            maxage3(currentindex) = maxage(time3);
            maxage2(currentindex) = maxage(time2n);
            actualage(currentindex) = age2;
            
            grp3(currentindex) = group(time3);
            grp2n(currentindex) = group(time2n);
            grp2o(currentindex) = group(time2n);
            grp1(currentindex) = 0.0;
            
            // The next indexer includes the following order: (1st # of age blocks) + 
            // (1st # of stage cols) + (1st # of age rows) + stage in time 3
            index321(currentindex) = currentindex;
            index21(currentindex) = time2n + ((age2 - firstage) * nostages);
            indatalong(currentindex) = 1.0;
            
            // This section identifies elements with non-zero entries by their
            // element number in the final matrix
            if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
              if (age2 >= minage2(currentindex) && age2 < maxage3(currentindex)) { 
                
                // Survival transitions
                aliveequal(currentindex) =
                  ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                  (time2n * (nostages - 1) * totalages) +
                  ((age3 - firstage) * (nostages - 1)) + time3;
              }
            }
            
            if (time3 < nostages_nodead && time2n < nostages_nodead) {
              
              if (repmatrix((time3 + (nostages_nodead * time2n))) > 0.0) {
                
                // Now fecundity
                age3 = firstage;
                currentindex = time3 + ((age3 - firstage) * nostages) + 
                  (time2n * nostages * totalages) +
                  ((age2 - firstage) * nostages * nostages * totalages);
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2n);
                stage1(currentindex) = 0.0;
                
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2n);
                size1(currentindex) = 0.0;
                
                sizeb3(currentindex) = binsizebctr(time3);
                sizeb2n(currentindex) = binsizebctr(time2n);
                sizeb2o(currentindex) = binsizebctr(time2n);
                sizeb1(currentindex) = 0.0;
                
                if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
                
                sizec3(currentindex) = binsizecctr(time3);
                sizec2n(currentindex) = binsizecctr(time2n);
                sizec2o(currentindex) = binsizecctr(time2n);
                sizec1(currentindex) = 0.0;
                
                if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
                if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
                if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                
                obs3(currentindex) = obsstatus(time3);
                obs2n(currentindex) = obsstatus(time2n);
                obs2o(currentindex) = obsstatus(time2n);
                obs1(currentindex) = 0.0;
                
                rep3(currentindex) = repstatus(time3);
                rep2n(currentindex) = repstatus(time2n);
                rep2o(currentindex) = repstatus(time2n);
                rep1(currentindex) = 0.0;
                
                mat3(currentindex) = matstatus(time3);
                mat2n(currentindex) = matstatus(time2n);
                mat2o(currentindex) = matstatus(time2n);
                mat1(currentindex) = 0.0;
                
                imm3(currentindex) = immstatus(time3);
                imm2n(currentindex) = immstatus(time2n);
                imm2o(currentindex) = immstatus(time2n);
                imm1(currentindex) = 0.0;
                
                if (rep2n(currentindex) > 0.0 && time3 < nostages_nodead && time2n < nostages_nodead) {
                  repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
                } else repentry3(currentindex) = 0.0;
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2n);
                indata1(currentindex) = 0.0;
                
                binwidth(currentindex) = binsizewidth(time3);
                binbwidth(currentindex) = binsizebwidth(time3);
                bincwidth(currentindex) = binsizecwidth(time3);
                
                if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
                if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2n);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2n);
                actualage(currentindex) = age2;
                
                grp3(currentindex) = group(time3);
                grp2n(currentindex) = group(time2n);
                grp2o(currentindex) = group(time2n);
                grp1(currentindex) = 0.0;
                
                // The next indexer includes the following order: (1st # of age blocks) + 
                // (1st # of stage cols) + (1st # of age rows) + stage in time 3
                index321(currentindex) = currentindex;
                index21(currentindex) = time2n + ((age2 - firstage) * nostages);
                indatalong(currentindex) = 1.0;
                
                // This section identifies elements with non-zero entries by their
                // element number in the final matrix
                if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
                  if (age2 >= minage2(currentindex) && age2 <= maxage2(currentindex)) { 
                    
                    // Fecundity transitions
                    aliveequal(currentindex) = 
                      ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                      (time2n * (nostages - 1) * totalages) +
                      ((age3 - firstage) * (nostages - 1)) + time3;
                  }
                } // if statement leading to aliveequal assignment
              } // if statement yielding fecundity estimation
            } // if statement checking time3 and time2n
          } // time3 loop
        } // time2n loop
      } else if (cont == 1) { // Self-loop on final age, if the organism can live past final age
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            // First survival
            age3 = age2;
            currentindex = time3 + ((age3 - firstage) * nostages) + 
              (time2n * nostages * totalages) +
              ((age2 - firstage) * nostages * nostages * totalages);
            
            stage3(currentindex) = newstageid(time3);
            stage2n(currentindex) = newstageid(time2n);
            stage2o(currentindex) = newstageid(time2n);
            stage1(currentindex) = 0.0;
            
            size3(currentindex) = binsizectr(time3);
            size2n(currentindex) = binsizectr(time2n);
            size2o(currentindex) = binsizectr(time2n);
            size1(currentindex) = 0.0;
            
            sizeb3(currentindex) = binsizebctr(time3);
            sizeb2n(currentindex) = binsizebctr(time2n);
            sizeb2o(currentindex) = binsizebctr(time2n);
            sizeb1(currentindex) = 0.0;
            
            if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
            if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
            
            sizec3(currentindex) = binsizecctr(time3);
            sizec2n(currentindex) = binsizecctr(time2n);
            sizec2o(currentindex) = binsizecctr(time2n);
            sizec1(currentindex) = 0.0;
            
            if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
            if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
            if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                
            obs3(currentindex) = obsstatus(time3);
            obs2n(currentindex) = obsstatus(time2n);
            obs2o(currentindex) = obsstatus(time2n);
            obs1(currentindex) = 0.0;
            
            rep3(currentindex) = repstatus(time3);
            rep2n(currentindex) = repstatus(time2n);
            rep2o(currentindex) = repstatus(time2n);
            rep1(currentindex) = 0.0;
            
            mat3(currentindex) = matstatus(time3);
            mat2n(currentindex) = matstatus(time2n);
            mat2o(currentindex) = matstatus(time2n);
            mat1(currentindex) = 0.0;
            
            imm3(currentindex) = immstatus(time3);
            imm2n(currentindex) = immstatus(time2n);
            imm2o(currentindex) = immstatus(time2n);
            imm1(currentindex) = 0.0;
            
            repentry3(currentindex) = 0.0;
            
            indata3(currentindex) = indata(time3);
            indata2n(currentindex) = indata(time2n);
            indata2o(currentindex) = indata(time2n);
            indata1(currentindex) = 0.0;
            
            binwidth(currentindex) = binsizewidth(time3);
            binbwidth(currentindex) = binsizebwidth(time3);
            bincwidth(currentindex) = binsizecwidth(time3);
            
            if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
            if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                
            minage3(currentindex) = minage(time3);
            minage2(currentindex) = minage(time2n);
            maxage3(currentindex) = maxage(time3);
            maxage2(currentindex) = maxage(time2n);
            actualage(currentindex) = age2;
            
            grp3(currentindex) = group(time3);
            grp2n(currentindex) = group(time2n);
            grp2o(currentindex) = group(time2n);
            grp1(currentindex) = 0.0;
            
            // The next indexer includes the following order: (1st # of age blocks) + 
            // (1st # of stage cols) + (1st # of age rows) + stage in time 3
            index321(currentindex) = currentindex;
            index21(currentindex) = time2n + ((age2 - firstage) * nostages);
            indatalong(currentindex) = 1;
            
            // This section identifies elements with non-zero entries by their
            // element number in the final matrix
            if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
              if (age2 >= minage2(currentindex) && age2 < maxage3(currentindex)) { 

                // Survival transitions
                aliveequal(currentindex) = 
                  ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                  (time2n * (nostages - 1) * totalages) +
                  ((age3 - firstage) * (nostages - 1)) + time3;
              }
            }
            
            if (time3 < nostages_nodead && time2n < nostages_nodead) {
              if (repmatrix((time3 + (nostages_nodead * time2n))) > 0.0) {
                
                // Now fecundity
                age3 = firstage;
                currentindex = time3 + ((age3 - firstage) * nostages) + 
                  (time2n * nostages * totalages) +
                  ((age2 - firstage) * nostages * nostages * totalages);
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2n);
                stage1(currentindex) = 0.0;
                
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2n);
                size1(currentindex) = 0.0;
                
                sizeb3(currentindex) = binsizebctr(time3);
                sizeb2n(currentindex) = binsizebctr(time2n);
                sizeb2o(currentindex) = binsizebctr(time2n);
                sizeb1(currentindex) = 0.0;
                
                if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0.0;
                if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0.0;
                
                sizec3(currentindex) = binsizecctr(time3);
                sizec2n(currentindex) = binsizecctr(time2n);
                sizec2o(currentindex) = binsizecctr(time2n);
                sizec1(currentindex) = 0.0;
                
                if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0.0;
                if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0.0;
                if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0.0;
                
                obs3(currentindex) = obsstatus(time3);
                obs2n(currentindex) = obsstatus(time2n);
                obs2o(currentindex) = obsstatus(time2n);
                obs1(currentindex) = 0.0;
                
                rep3(currentindex) = repstatus(time3);
                rep2n(currentindex) = repstatus(time2n);
                rep2o(currentindex) = repstatus(time2n);
                rep1(currentindex) = 0.0;
                
                mat3(currentindex) = matstatus(time3);
                mat2n(currentindex) = matstatus(time2n);
                mat2o(currentindex) = matstatus(time2n);
                mat1(currentindex) = 0.0;
                
                imm3(currentindex) = immstatus(time3);
                imm2n(currentindex) = immstatus(time2n);
                imm2o(currentindex) = immstatus(time2n);
                imm1(currentindex) = 0.0;
                
                if (rep2n(currentindex) == 1) {
                  repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
                } else repentry3(currentindex) = 0.0;
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2n);
                indata1(currentindex) = 0.0;
                
                binwidth(currentindex) = binsizewidth(time3);
                binbwidth(currentindex) = binsizebwidth(time3);
                bincwidth(currentindex) = binsizecwidth(time3);
                
                if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0.0;
                if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0.0;
                
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2n);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2n);
                actualage(currentindex) = age2;
                grp3(currentindex) = group(time3);
                grp2n(currentindex) = group(time2n);
                grp2o(currentindex) = group(time2n);
                grp1(currentindex) = 0.0;
                
                // The next indexer includes the following order: (1st # of age blocks) +
                // (1st # of stage cols) + (1st # of age rows) + stage in time 3
                index321(currentindex) = currentindex;
                index21(currentindex) = time2n + ((age2 - firstage) * nostages);
                indatalong(currentindex) = 1.0;
                
                // This section identifies elements with non-zero entries by their
                // element number in the final matrix
                if (alive(time2n) == 1.0 && alive(time3) == 1.0) {
                  if (age2 >= minage2(currentindex) && age2 <= maxage2(currentindex)) { 
                    
                    // Fecundity transitions
                    aliveequal(currentindex) =
                      ((age2 - firstage) * (nostages - 1) * (nostages - 1) * totalages) + 
                      (time2n * (nostages - 1) * totalages) +
                      ((age3 - firstage) * (nostages - 1)) + time3;
                  }
                } // if statement leading to aliveequal assignment
              } // if statement yielding fecundity estimation
            } // if statement checking time3 and time2n
          } // time3 loop
        } // time2n loop
      }// if-else statement
    } // age2 loop
    
    if (ovrows > 1 || ovconvtype(0) != -1.0) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtypeage,
        ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      ovsurvmult = asadditions.col(5);
      ovfecmult = asadditions.col(6);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0.0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } // Age-by-stage loop (style = 2)
  
  // Now the final output formatting
  Rcpp::List output_longlist(59);
  int stage3_length = 0;
  
  arma::uvec used_indices;
  
  if (filter == 1) {
    used_indices = find(index321 != -1.0);
  } else if (filter == 2) {
    used_indices = find(aliveequal != -1.0);
  }
  
  if (filter > 0) {
    int new_length = used_indices.n_elem;
    stage3_length = new_length;
    
    NumericVector stage3_new(new_length);
    NumericVector stage2n_new(new_length);
    NumericVector stage2o_new(new_length);
    NumericVector stage1_new(new_length);
    
    NumericVector size3_new(new_length);
    NumericVector size2n_new(new_length);
    NumericVector size2o_new(new_length);
    NumericVector size1_new(new_length);
    
    NumericVector sizeb3_new(new_length);
    NumericVector sizeb2n_new(new_length);
    NumericVector sizeb2o_new(new_length);
    NumericVector sizeb1_new(new_length);
    
    NumericVector sizec3_new(new_length);
    NumericVector sizec2n_new(new_length);
    NumericVector sizec2o_new(new_length);
    NumericVector sizec1_new(new_length);
    
    NumericVector obs3_new(new_length);
    NumericVector obs2n_new(new_length);
    NumericVector obs2o_new(new_length);
    NumericVector obs1_new(new_length);
    
    NumericVector rep3_new(new_length);
    NumericVector rep2n_new(new_length);
    NumericVector rep2o_new(new_length);
    NumericVector rep1_new(new_length);
    
    NumericVector mat3_new(new_length);
    NumericVector mat2n_new(new_length);
    NumericVector mat2o_new(new_length);
    NumericVector mat1_new(new_length);
    
    NumericVector imm3_new(new_length);
    NumericVector imm2n_new(new_length);
    NumericVector imm2o_new(new_length);
    NumericVector imm1_new(new_length);
    
    NumericVector repentry3_new(new_length);
    NumericVector indata3_new(new_length);
    NumericVector indata2n_new(new_length);
    NumericVector indata2o_new(new_length);
  
    NumericVector indata1_new(new_length);
    NumericVector binwidth_new(new_length);
    NumericVector binbwidth_new(new_length);
    NumericVector bincwidth_new(new_length);
    
    NumericVector minage3_new(new_length);
    NumericVector minage2_new(new_length);
    NumericVector maxage3_new(new_length);
    NumericVector maxage2_new(new_length);
    NumericVector actualage_new(new_length);
    
    NumericVector grp3_new(new_length);
    NumericVector grp2n_new(new_length);
    NumericVector grp2o_new(new_length);
    NumericVector grp1_new(new_length);
    
    NumericVector indatalong_new(new_length);
    NumericVector ovgivent_new(new_length);
    NumericVector ovestt_new(new_length);
    NumericVector ovgivenf_new(new_length);
  
    NumericVector ovestf_new(new_length);
    NumericVector ovsurvmult_new(new_length);
    NumericVector ovfecmult_new(new_length);
    
    NumericVector aliveequal_new(new_length);
    NumericVector index321_new(new_length);
    NumericVector index21_new(new_length);
    
    for (int i = 0; i < new_length; i++) {
      stage3_new(i) = stage3(used_indices(i));
      stage2n_new(i) = stage2n(used_indices(i));
      stage2o_new(i) = stage2o(used_indices(i));
      stage1_new(i) = stage1(used_indices(i));
      
      size3_new(i) = size3(used_indices(i));
      size2n_new(i) = size2n(used_indices(i));
      size2o_new(i) = size2o(used_indices(i));
      size1_new(i) = size1(used_indices(i));
      
      sizeb3_new(i) = sizeb3(used_indices(i));
      sizeb2n_new(i) = sizeb2n(used_indices(i));
      sizeb2o_new(i) = sizeb2o(used_indices(i));
      sizeb1_new(i) = sizeb1(used_indices(i));
      
      sizec3_new(i) = sizec3(used_indices(i));
      sizec2n_new(i) = sizec2n(used_indices(i));
      sizec2o_new(i) = sizec2o(used_indices(i));
      sizec1_new(i) = sizec1(used_indices(i));
      
      obs3_new(i) = obs3(used_indices(i));
      obs2n_new(i) = obs2n(used_indices(i));
      obs2o_new(i) = obs2o(used_indices(i));
      obs1_new(i) = obs1(used_indices(i));
      
      rep3_new(i) = rep3(used_indices(i));
      rep2n_new(i) = rep2n(used_indices(i));
      rep2o_new(i) = rep2o(used_indices(i));
      rep1_new(i) = rep1(used_indices(i));
      
      mat3_new(i) = mat3(used_indices(i));
      mat2n_new(i) = mat2n(used_indices(i));
      mat2o_new(i) = mat2o(used_indices(i));
      mat1_new(i) = mat1(used_indices(i));
      
      imm3_new(i) = imm3(used_indices(i));
      imm2n_new(i) = imm2n(used_indices(i));
      imm2o_new(i) = imm2o(used_indices(i));
      imm1_new(i) = imm1(used_indices(i));
      
      repentry3_new(i) = repentry3(used_indices(i));
      indata3_new(i) = indata3(used_indices(i));
      indata2n_new(i) = indata2n(used_indices(i));
      indata2o_new(i) = indata2o(used_indices(i));
    
      indata1_new(i) = indata1(used_indices(i));
      binwidth_new(i) = binwidth(used_indices(i));
      binbwidth_new(i) = binbwidth(used_indices(i));
      bincwidth_new(i) = bincwidth(used_indices(i));
      
      minage3_new(i) = minage3(used_indices(i));
      minage2_new(i) = minage2(used_indices(i));
      maxage3_new(i) = maxage3(used_indices(i));
      maxage2_new(i) = maxage2(used_indices(i));
      actualage_new(i) = actualage(used_indices(i));
      
      grp3_new(i) = grp3(used_indices(i));
      grp2n_new(i) = grp2n(used_indices(i));
      grp2o_new(i) = grp2o(used_indices(i));
      grp1_new(i) = grp1(used_indices(i));
      
      indatalong_new(i) = indatalong(used_indices(i));
      ovgivent_new(i) = ovgivent(used_indices(i));
      ovestt_new(i) = ovestt(used_indices(i));
      ovgivenf_new(i) = ovgivenf(used_indices(i));
    
      ovestf_new(i) = ovestf(used_indices(i));
      ovsurvmult_new(i) = ovsurvmult(used_indices(i));
      ovfecmult_new(i) = ovfecmult(used_indices(i));
      
      aliveequal_new(i) = aliveequal(used_indices(i));
      index321_new(i) = index321(used_indices(i));
      index21_new(i) = index21(used_indices(i));
    }
    
    output_longlist(0) = stage3_new;
    output_longlist(1) = stage2n_new;
    output_longlist(2) = stage2o_new;
    output_longlist(3) = stage1_new;
    output_longlist(4) = size3_new;
    output_longlist(5) = size2n_new;
    output_longlist(6) = size2o_new;
    output_longlist(7) = size1_new;
    output_longlist(8) = sizeb3_new;
    output_longlist(9) = sizeb2n_new;
    
    output_longlist(10) = sizeb2o_new;
    output_longlist(11) = sizeb1_new;
    output_longlist(12) = sizec3_new;
    output_longlist(13) = sizec2n_new;
    output_longlist(14) = sizec2o_new;
    output_longlist(15) = sizec1_new;
    output_longlist(16) = obs3_new;
    output_longlist(17) = obs2n_new;
    output_longlist(18) = obs2o_new;
    output_longlist(19) = obs1_new;
    
    output_longlist(20) = rep3_new;
    output_longlist(21) = rep2n_new;
    output_longlist(22) = rep2o_new;
    output_longlist(23) = rep1_new;
    output_longlist(24) = mat3_new;
    output_longlist(25) = mat2n_new;
    output_longlist(26) = mat2o_new;
    output_longlist(27) = mat1_new;
    output_longlist(28) = imm3_new;
    output_longlist(29) = imm2n_new;
    
    output_longlist(30) = imm2o_new;
    output_longlist(31) = imm1_new;
    output_longlist(32) = repentry3_new;
    output_longlist(33) = indata3_new;
    output_longlist(34) = indata2n_new;
    output_longlist(35) = indata2o_new;
    output_longlist(36) = indata1_new;
    output_longlist(37) = binwidth_new;
    output_longlist(38) = binbwidth_new;
    output_longlist(39) = bincwidth_new;
    
    output_longlist(40) = minage3_new;
    output_longlist(41) = minage2_new;
    output_longlist(42) = maxage3_new;
    output_longlist(43) = maxage2_new;
    output_longlist(44) = actualage_new;
    
    output_longlist(45) = grp3_new;
    output_longlist(46) = grp2n_new;
    output_longlist(47) = grp2o_new;
    output_longlist(48) = grp1_new;
    
    output_longlist(49) = indatalong_new;
    output_longlist(50) = ovgivent_new;
    output_longlist(51) = ovestt_new;
    output_longlist(52) = ovgivenf_new;
    output_longlist(53) = ovestf_new;
    output_longlist(54) = ovsurvmult_new;
    output_longlist(55) = ovfecmult_new;
    
    output_longlist(56) = aliveequal_new;
    output_longlist(57) = index321_new;
    output_longlist(58) = index21_new;
    
  } else {
    stage3_length = stage3.n_elem;
    
    output_longlist(0) = Rcpp::NumericVector(stage3.begin(), stage3.end());
    output_longlist(1) = Rcpp::NumericVector(stage2n.begin(), stage2n.end());
    output_longlist(2) = Rcpp::NumericVector(stage2o.begin(), stage2o.end());
    output_longlist(3) = Rcpp::NumericVector(stage1.begin(), stage1.end());
    output_longlist(4) = Rcpp::NumericVector(size3.begin(), size3.end());
    output_longlist(5) = Rcpp::NumericVector(size2n.begin(), size2n.end());
    output_longlist(6) = Rcpp::NumericVector(size2o.begin(), size2o.end());
    output_longlist(7) = Rcpp::NumericVector(size1.begin(), size1.end());
    output_longlist(8) = Rcpp::NumericVector(sizeb3.begin(), sizeb3.end());
    output_longlist(9) = Rcpp::NumericVector(sizeb2n.begin(), sizeb2n.end());
    
    output_longlist(10) = Rcpp::NumericVector(sizeb2o.begin(), sizeb2o.end());
    output_longlist(11) = Rcpp::NumericVector(sizeb1.begin(), sizeb1.end());
    output_longlist(12) = Rcpp::NumericVector(sizec3.begin(), sizec3.end());
    output_longlist(13) = Rcpp::NumericVector(sizec2n.begin(), sizec2n.end());
    output_longlist(14) = Rcpp::NumericVector(sizec2o.begin(), sizec2o.end());
    output_longlist(15) = Rcpp::NumericVector(sizec1.begin(), sizec1.end());
    output_longlist(16) = Rcpp::NumericVector(obs3.begin(), obs3.end());
    output_longlist(17) = Rcpp::NumericVector(obs2n.begin(), obs2n.end());
    output_longlist(18) = Rcpp::NumericVector(obs2o.begin(), obs2o.end());
    output_longlist(19) = Rcpp::NumericVector(obs1.begin(), obs1.end());
    
    output_longlist(20) = Rcpp::NumericVector(rep3.begin(), rep3.end());
    output_longlist(21) = Rcpp::NumericVector(rep2n.begin(), rep2n.end());
    output_longlist(22) = Rcpp::NumericVector(rep2o.begin(), rep2o.end());
    output_longlist(23) = Rcpp::NumericVector(rep1.begin(), rep1.end());
    output_longlist(24) = Rcpp::NumericVector(mat3.begin(), mat3.end());
    output_longlist(25) = Rcpp::NumericVector(mat2n.begin(), mat2n.end());
    output_longlist(26) = Rcpp::NumericVector(mat2o.begin(), mat2o.end());
    output_longlist(27) = Rcpp::NumericVector(mat1.begin(), mat1.end());
    output_longlist(28) = Rcpp::NumericVector(imm3.begin(), imm3.end());
    output_longlist(29) = Rcpp::NumericVector(imm2n.begin(), imm2n.end());
    
    output_longlist(30) = Rcpp::NumericVector(imm2o.begin(), imm2o.end());
    output_longlist(31) = Rcpp::NumericVector(imm1.begin(), imm1.end());
    output_longlist(32) = Rcpp::NumericVector(repentry3.begin(), repentry3.end());
    output_longlist(33) = Rcpp::NumericVector(indata3.begin(), indata3.end());
    output_longlist(34) = Rcpp::NumericVector(indata2n.begin(), indata2n.end());
    output_longlist(35) = Rcpp::NumericVector(indata2o.begin(), indata2o.end());
    output_longlist(36) = Rcpp::NumericVector(indata1.begin(), indata1.end());
    output_longlist(37) = Rcpp::NumericVector(binwidth.begin(), binwidth.end());
    output_longlist(38) = Rcpp::NumericVector(binbwidth.begin(), binbwidth.end());
    output_longlist(39) = Rcpp::NumericVector(bincwidth.begin(), bincwidth.end());
    
    output_longlist(40) = Rcpp::NumericVector(minage3.begin(), minage3.end());
    output_longlist(41) = Rcpp::NumericVector(minage2.begin(), minage2.end());
    output_longlist(42) = Rcpp::NumericVector(maxage3.begin(), maxage3.end());
    output_longlist(43) = Rcpp::NumericVector(maxage2.begin(), maxage2.end());
    output_longlist(44) = Rcpp::NumericVector(actualage.begin(), actualage.end());
    
    output_longlist(45) = Rcpp::NumericVector(grp3.begin(), grp3.end());
    output_longlist(46) = Rcpp::NumericVector(grp2n.begin(), grp2n.end());
    output_longlist(47) = Rcpp::NumericVector(grp2o.begin(), grp2o.end());
    output_longlist(48) = Rcpp::NumericVector(grp1.begin(), grp1.end());
    
    output_longlist(49) = Rcpp::NumericVector(indatalong.begin(), indatalong.end());
    output_longlist(50) = Rcpp::NumericVector(ovgivent.begin(), ovgivent.end());
    output_longlist(51) = Rcpp::NumericVector(ovestt.begin(), ovestt.end());
    output_longlist(52) = Rcpp::NumericVector(ovgivenf.begin(), ovgivenf.end());
    output_longlist(53) = Rcpp::NumericVector(ovestf.begin(), ovestf.end());
    output_longlist(54) = Rcpp::NumericVector(ovsurvmult.begin(), ovsurvmult.end());
    output_longlist(55) = Rcpp::NumericVector(ovfecmult.begin(), ovfecmult.end());
    
    output_longlist(56) = Rcpp::NumericVector(aliveequal.begin(), aliveequal.end());
    output_longlist(57) = Rcpp::NumericVector(index321.begin(), index321.end());
    output_longlist(58) = Rcpp::NumericVector(index21.begin(), index21.end());
  }
  
  CharacterVector namevec = {"stage3", "stage2n", "stage2o", "stage1", "size3",
    "size2n", "size2o", "size1", "sizeb3", "sizeb2n", "sizeb2o", "sizeb1", 
    "sizec3", "sizec2n", "sizec2o", "sizec1", "obs3", "obs2n", "obs2o", "obs1",
    "rep3", "rep2n", "rep2o", "rep1", "mat3", "mat2n", "mat2o", "mat1", "imm3",
    "imm2n", "imm2o", "imm1", "repentry3", "indata3", "indata2n", "indata2o",
    "indata1", "binwidth", "binbwidth", "bincwidth", "minage3", "minage2",
    "maxage3", "maxage2", "actualage", "group3", "group2n", "group2o", "group1",
    "indata", "ovgiven_t", "ovest_t", "ovgiven_f", "ovest_f", "ovsurvmult",
    "ovfecmult", "aliveandequal", "index321", "index21"};
  output_longlist.attr("names") = namevec;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stage3_length);
  output_longlist.attr("class") = "data.frame";
  
  return output_longlist;
}

//' Estimate All Elements of Raw Historical Matrix
//' 
//' Function \code{.specialpatrolgroup()} swiftly calculates matrix transitions
//' in raw historical matrices, and serves as the core workhorse function behind
//' \code{\link{rlefko3}()}.
//' 
//' @name specialpatrolgroup
//' 
//' @param sge9l The Allstages data frame developed for \code{rlefko3()}
//' covering stage pairs across times \emph{t}+1, \emph{t} and \emph{t}-1.
//' Generally termed \code{stageexpansion9}.
//' @param sge3 The data frame covering all stages in times \emph{t} and
//' \emph{t}-1. Generally termed \code{stageexpansion3}.
//' @param MainData The demographic dataset modified to hold \code{usedfec}
//' columns.
//' @param StageFrame The full stageframe for the analysis.
//' @param repmatrix The modified repmatrix used in the course of computation.
//' This is used particularly when deVries-format hMPMs are desired.
//' @param format Indicates whether to output Ehrlen-format hMPMs (1) or
//' deVries-format hMPMs (2).
//' @param err_switch If set to 1, then will also output probsrates and
//' stage2fec.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.specialpatrolgroup)]]
List specialpatrolgroup(DataFrame sge9l, DataFrame sge3, DataFrame MainData,
  DataFrame StageFrame, int format, int err_switch) {
  
  arma::vec sge9stage3 = sge9l["stage3"];
  arma::vec sge9fec32 = sge9l["repentry3"];
  arma::vec sge9rep2 = sge9l["rep2o"];
  arma::vec sge9indata32 = sge9l["indata"];
  arma::vec sge9ovgivent = sge9l["ovgiven_t"];
  arma::vec sge9ovgivenf = sge9l["ovgiven_f"];
  arma::vec sge9ovestt = sge9l["ovest_t"];
  arma::vec sge9ovestf = sge9l["ovest_f"];
  arma::vec sge9ovsurvmult = sge9l["ovsurvmult"];
  arma::vec sge9ovfecmult = sge9l["ovfecmult"];
  arma::vec sge9index321 = sge9l["index321"];
  arma::vec sge9index21 = sge9l["index21"]; // This is sge92index - not sure if this is needed
  arma::vec aliveandequal = sge9l["aliveandequal"];
  
  arma::vec sge3rep2 = sge3["rep2n"];
  arma::vec sge3fec32 = sge3["fec32n"];
  arma::vec sge3index21 = sge3["index21"];
  arma::vec sge3stage2n = sge3["stage2n"];
  arma::vec sge3stage3 = sge3["stage3"];
  
  arma::vec dataindex321 = MainData["index321"];
  arma::vec dataindex21 = MainData["pairindex21"];
  arma::vec dataalive3 = MainData["alive3"];
  arma::vec datausedfec2 = MainData["usedfec2"];
  arma::vec dataindex3 = MainData["index3"];
  arma::vec dataindex2 = MainData["index2"];
  arma::vec dataindex1 = MainData["index1"];

  arma::vec sfsizes = StageFrame["sizebin_center"];
  int nostages = sfsizes.n_elem;
  
  int n = dataindex321.n_elem;
  int no2stages = nostages - 1;
  int noelems = sge9index321.n_elem;
  
  if (format == 2) {
    no2stages = no2stages - 1;
  }
  
  int matrixdim = (nostages - 1) * no2stages;
  
  arma::vec probsrates0(noelems); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1(noelems); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2(noelems); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3(noelems); // 4th vec = total fec for pair stage
  probsrates0.zeros();
  probsrates1.zeros();
  probsrates2.zeros();
  probsrates3.zeros();
  
  arma::mat stage2fec(sge3index21.n_elem, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  // These next structures develop the prior stage
  arma::vec probsrates0p(noelems, fill::zeros); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1p(noelems, fill::zeros); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2p(noelems, fill::zeros); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3p(noelems, fill::zeros); // 4th vec = total fec for pair stage
  
  arma::mat stage2fecp(sge3index21.n_elem, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  // The final matrices, though empty
  arma::mat tmatrix(matrixdim, matrixdim, fill::zeros); // Main output U matrix
  arma::mat fmatrix(matrixdim, matrixdim, fill::zeros); // Main output F matrix
  
  arma::uvec all_repentries = find(sge9fec32 > 0);
  arma::vec all_entry_stages = arma::unique(sge9stage3(all_repentries));
  int aes_count = all_entry_stages.n_elem;
  
  arma::mat dataindex321_prior(n, aes_count);
  dataindex321_prior.fill(-1);
  
  // This section creates an alternative index for use in fecundity calculations under deVries format
  if (format == 2) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < aes_count; j++) {
        dataindex321_prior(i, j) = (all_entry_stages(j) - 1) + ((nostages - 2) * nostages) + 
          ((dataindex2(i) - 1) * nostages * nostages) + 
          ((dataindex1(i) - 1) * nostages * nostages * nostages);
      }
    }
  }
  
  // This main loop counts individuals going through transitions and sums their
  // fecundities, and then adds that info to the 3-trans and 2-trans tables
  for (int i = 0; i < n; i++) { 
    arma::uvec choiceelement = find(sge9index321 == dataindex321(i)); // Added this now
    
    stage2fec((dataindex21(i)), 0) = stage2fec((dataindex21(i)), 0) + 1; // Yields sum of all individuals with particular transition
    
    if (choiceelement.n_elem > 0) {
      probsrates0(choiceelement(0)) = probsrates0(choiceelement(0)) + 1; // Yields sum of all individuals with particular transition
      
      if (dataalive3(i) > 0) {
        stage2fec((dataindex21(i)), 1) = stage2fec((dataindex21(i)), 1) + 1;
      }
      
      stage2fec((dataindex21(i)), 2) = stage2fec((dataindex21(i)), 2) + datausedfec2(i);
    }
    
    if (format == 2) {
      for (int j = 0; j < aes_count; j++) {
        arma::uvec choiceelementp = find(sge9index321 == dataindex321_prior(i, j));
        
        stage2fecp((dataindex21(i)), 0) = stage2fecp((dataindex21(i)), 0) + 1;
      
        if (choiceelementp.n_elem > 0) {
          probsrates0p(choiceelementp(0)) = probsrates0p(choiceelementp(0)) + 1;
          
          if (dataalive3(i) > 0) {
            stage2fecp((dataindex21(i)), 1) = stage2fecp((dataindex21(i)), 1) + 1;
          }
          
          stage2fecp((dataindex21(i)), 2) = stage2fecp((dataindex21(i)), 2) + datausedfec2(i);
        }
      }
    }
  }
  
  // The next bit puts together core data to be used to estimate matrix elements
  for (int i = 0; i < noelems; i++) {
    int baseindex21 = sge9index21(i);
    
    if (baseindex21 > -1) {
      arma::uvec coreelementsforchoice = find(sge3index21 == baseindex21);
      unsigned int thechosenone = coreelementsforchoice(0);
      
      probsrates1(i) = stage2fec(thechosenone, 0);
      probsrates2(i) = stage2fec(thechosenone, 1);
      probsrates3(i) = stage2fec(thechosenone, 2);
      
      if (format == 2) {
        arma::uvec coreelementsforchoicep = find(sge3index21 == baseindex21);
        unsigned int thechosenonep = coreelementsforchoicep(0);
        
        probsrates1p(i) = stage2fecp(thechosenonep, 0);
        probsrates2p(i) = stage2fecp(thechosenonep, 1);
        probsrates3p(i) = stage2fecp(thechosenonep, 2);
      }
    }
  }
  
  // Here we create the matrices
  for (int elem3 = 0; elem3 < noelems; elem3++) {
    
    if (aliveandequal(elem3) != -1) {
      if (sge9ovsurvmult(elem3) < 0) sge9ovsurvmult(elem3) = 1.0;
      
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3)* sge9ovsurvmult(elem3) /
        probsrates1(elem3); // Survival
      
      // Fecundity
      if (sge9ovfecmult(elem3) < 0) sge9ovfecmult(elem3) = 1.0;
      if (format == 2) {
        fmatrix(aliveandequal(elem3)) = sge9fec32(elem3) * sge9rep2(elem3) *
          probsrates3p(elem3) * sge9ovfecmult(elem3) / probsrates1p(elem3);
      } else {
        fmatrix(aliveandequal(elem3)) = sge9fec32(elem3) * sge9rep2(elem3) *
          probsrates3(elem3) * sge9ovfecmult(elem3) / probsrates1(elem3);
      }
    }
  }
  
  // Now we will correct transitions and rates for given stuff
  arma::uvec ovgiventind = find(sge9ovgivent != -1);
  arma::uvec ovgivenfind = find(sge9ovgivenf != -1);
  int ovgtn = ovgiventind.n_elem;
  int ovgfn = ovgivenfind.n_elem;
  
  if (ovgtn > 0) {
    for (int i = 0; i < ovgtn; i++) {
      int matrixelement2 = aliveandequal(ovgiventind(i));
      
      tmatrix(matrixelement2) = sge9ovgivent(ovgiventind(i));
    }
  }
  
  if (ovgfn > 0) {
    for (int i = 0; i < ovgfn; i++) {
      int matrixelement2 = aliveandequal(ovgivenfind(i));
      
      fmatrix(matrixelement2) = sge9ovgivenf(ovgivenfind(i));
    }
  }
  
  // This section replaces transitions for proxy values as given in the overwrite table  
  arma::uvec ovesttind = find(sge9ovestt != -1);
  arma::uvec ovestfind = find(sge9ovestf != -1);
  int ovestn = ovesttind.n_elem;
  int ovesfn = ovestfind.n_elem;
  
  if (ovestn > 0) {
    for (int i = 0; i < ovestn; i++) {
      arma::uvec replacement = find(sge9index321 == sge9ovestt(ovesttind(i)));
      
      if (replacement.n_elem > 0) {
        tmatrix(aliveandequal(ovesttind(i))) = tmatrix(aliveandequal(replacement(0)));
      }
      
    }
  }
  
  if (ovesfn > 0) {
    for (int i = 0; i < ovesfn; i++) {
      arma::uvec replacement = find(sge9index321 == sge9ovestf(ovestfind(i)));
      
      if (replacement.n_elem > 0) {
        fmatrix(aliveandequal(ovestfind(i))) = fmatrix(aliveandequal(replacement(0)));
      }
    }
  }
  
  // The next bit changes NAs to 0
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  if (err_switch == 1) {
    arma::mat concatenated_crap = arma::join_horiz(sge9index321, aliveandequal);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates0);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates1);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates2);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates3);
    concatenated_crap = arma::join_horiz(concatenated_crap, sge9fec32);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates0p);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates1p);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates2p);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates3p);
    
    arma::mat s2f = arma::join_horiz(stage2fec, stage2fecp);

    return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix,
      _["concrp"] = concatenated_crap, _["s2f"] = s2f, _["dataprior"] = dataindex321_prior);
  } else {
    return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
  }
}

//' Estimate All Elements of Raw Ahistorical Population Projection Matrix
//' 
//' Function \code{.normalpatrolgroup()} swiftly calculates matrix transitions
//' in raw ahistorical matrices, and serves as the core workhorse function
//' behind \code{\link{rlefko2}()}.
//' 
//' @name normalpatrolgroup
//' 
//' @param sge3 The Allstages data frame developed for \code{rlefko2()} covering
//' stage pairs across times \emph{t}+1 and \emph{t}. Generally termed
//' \code{stageexpansion3}.
//' @param sge2 The data frame covering all stages in time \emph{t}. Generally
//' termed \code{stageexpansion2}.
//' @param MainData The demographic dataset modified to hold \code{usedfec} and
//' \code{usedstage} columns.
//' @param StageFrame The full stageframe for the analysis.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.normalpatrolgroup)]]
List normalpatrolgroup(DataFrame sge3, DataFrame sge2, DataFrame MainData,
  DataFrame StageFrame) {
  
  arma::vec sge3fec32 = sge3["repentry3"];
  arma::vec sge3rep2 = sge3["rep2n"];
  arma::vec sge3indata32 = sge3["indata"];
  arma::vec sge3ovgivent = sge3["ovgiven_t"];
  arma::vec sge3ovgivenf = sge3["ovgiven_f"];
  arma::vec sge3ovestt = sge3["ovest_t"];
  arma::vec sge3ovestf = sge3["ovest_f"];
  arma::vec sge3ovsurvmult = sge3["ovsurvmult"];
  arma::vec sge3ovfecmult = sge3["ovfecmult"];
  arma::vec sge3index32 = sge3["index321"];
  arma::vec sge3index2 = sge3["stage2n"];
  arma::vec aliveandequal = sge3["aliveandequal"];
  
  arma::vec sge2rep2 = sge2["rep2"];
  arma::vec sge2fec3 = sge2["fec3"];
  arma::vec sge2index2 = sge2["index2"];
  arma::vec sge2stage2 = sge2["stage2"];
  
  arma::vec dataindex32 = MainData["index32"];
  arma::vec dataindex2 = MainData["index2"];
  arma::vec dataalive3 = MainData["alive3"];
  arma::vec datausedfec2 = MainData["usedfec2"];
  
  arma::vec sfsizes = StageFrame["sizebin_center"];
  int nostages = sfsizes.n_elem;
  
  int n = dataindex32.n_elem;
  int no2stages = sge2index2.n_elem - 1; // The -1 removes the dead stage, which is still within sge2
  int noelems = sge3index32.n_elem;
  
  arma::vec probsrates0(noelems, fill::zeros); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1(noelems, fill::zeros); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2(noelems, fill::zeros); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3(noelems, fill::zeros); // 4th vec = total fec for pair stage
  
  arma::mat stage2fec(no2stages, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  arma::mat tmatrix((nostages-1), (nostages-1), fill::zeros); // Main output U matrix
  arma::mat fmatrix((nostages-1), (nostages-1), fill::zeros); // Main output F matrix
  
  // This main loop counts individuals going through transitions and sums their
  // fecundities, and then adds that info to the 3-trans and 2-trans tables
  for (int i = 0; i < n; i++) { 
    
    // The next line yields sum of all individuals with particular transition
    probsrates0(dataindex32(i)) = probsrates0(dataindex32(i)) + 1; 
    
    // The next line yields sum of all individuals with particular transition
    stage2fec((dataindex2(i)), 0) = stage2fec((dataindex2(i)), 0) + 1; 
    if (dataalive3(i) > 0) {
      stage2fec((dataindex2(i)), 1) = stage2fec((dataindex2(i)), 1) + 1;
    }
    
    stage2fec((dataindex2(i)), 2) = stage2fec((dataindex2(i)), 2) + datausedfec2(i);
    
  }
  
  // This next loop populates vectors of individuals according to stage in time t
  for (int i = 0; i < no2stages; i++) {
    unsigned int foradding = ((sge2stage2(i) - 1) * nostages);
    
    for (int j = 0; j < nostages; j++) {
      unsigned int entry = foradding + j;
      
      probsrates1(entry) = stage2fec(i, 0);
      probsrates2(entry) = stage2fec(i, 1);
      probsrates3(entry) = stage2fec(i, 2);
    }
  }
  
  // Here we populate the main U and F matrices
  for (int elem3 = 0; elem3 < noelems; elem3++) {
    
    if (aliveandequal(elem3) != -1) {
      
      // The next lines DO leave NaNs in the matrices when 0 individuals are summed through in probsrates1
      if (sge3ovsurvmult(elem3) < 0) sge3ovsurvmult(elem3) = 1.0;
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
        probsrates1(elem3);
        
      if (sge3ovfecmult(elem3) < 0) sge3ovfecmult(elem3) = 1.0;
      fmatrix(aliveandequal(elem3)) = sge3fec32(elem3) * sge3rep2(elem3) * 
        probsrates3(elem3) * sge3ovfecmult(elem3) / probsrates1(elem3);
    }
  }
  
  // This section corrects for transitions given in the overwrite table
  arma::uvec ovgiventind = find(sge3ovgivent != -1);
  arma::uvec ovgivenfind = find(sge3ovgivenf != -1);
  int ovgtn = ovgiventind.n_elem;
  int ovgfn = ovgivenfind.n_elem;
  
  if (ovgtn > 0) {
    for (int i = 0; i < ovgtn; i++) {
      int matrixelement2 = aliveandequal(ovgiventind(i));
      
      tmatrix(matrixelement2) = sge3ovgivent(ovgiventind(i));
    }
  }
  
  if (ovgfn > 0) {
    for (int i = 0; i < ovgfn; i++) {
      int matrixelement2 = aliveandequal(ovgivenfind(i));
      
      fmatrix(matrixelement2) = sge3ovgivenf(ovgivenfind(i));
    }
  }
  
  // This section replaces transitions with proxies as given in the overwrite table
  arma::uvec ovesttind = find(sge3ovestt != -1);
  arma::uvec ovestfind = find(sge3ovestf != -1);
  int ovestn = ovesttind.n_elem;
  int ovesfn = ovestfind.n_elem;
  
  if (ovestn > 0) {
    for (int i = 0; i < ovestn; i++) {
      arma::uvec replacement = find(sge3index32 == sge3ovestt(ovesttind(i)));
      
      tmatrix(aliveandequal(ovesttind(i))) = tmatrix(aliveandequal(replacement(0)));
    }
  }
  
  if (ovesfn > 0) {
    for (int i = 0; i < ovesfn; i++) {
      arma::uvec replacement = find(sge3index32 == sge3ovestf(ovestfind(i)));
      
      fmatrix(aliveandequal(ovestfind(i))) = fmatrix(aliveandequal(replacement(0)));
    }
  }
  
  // The next bit changes NAs to 0
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
}

//' Estimate All Elements of Raw Ahistorical Population Projection Matrix
//' 
//' Function \code{.minorpatrolgroup()} swiftly calculates matrix transitions
//' in raw Leslie MPMs, and is used internally in \code{\link{rleslie}()}.
//' 
//' @name minorpatrolgroup
//' 
//' @param MainData The demographic dataset modified internally to have needed
//' variables for living status, reproduction status, and fecundity.
//' @param StageFrame The full stageframe for the analysis.
//' @param fectime An integer coding to estimate fecundity using time \emph{t}
//' (\code{2}) or time \emph{t}+1 \code{(3)}.
//' @param cont Should a self-loop transition be estimated for the final age.
//' @param lastage An integer coding for the last age to use in matrix
//' construction.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.minorpatrolgroup)]]
Rcpp::List minorpatrolgroup(DataFrame MainData, DataFrame StageFrame,
  int fectime, bool cont, double fec_mod) {
  
  arma::ivec data_age = MainData["usedobsage"];
  arma::vec data_alive2 = MainData["usedalive2"];
  arma::vec data_alive3 = MainData["usedalive3"];
  arma::vec data_usedfec = MainData["usedfec2"];
  arma::vec data_usedfec3 = MainData["usedfec3"];
  arma::vec data_usedrepst2 = MainData["usedrepst2"];
  arma::vec data_usedrepst3 = MainData["usedrepst3"];

  if (fectime == 3) {
    data_usedfec = data_usedfec3;
  }
  
  IntegerVector sf_minage = StageFrame["min_age"];
  IntegerVector sf_maxage = StageFrame["max_age"];
  IntegerVector sf_repstatus = StageFrame["repstatus"];
  int noages = sf_minage.length();
  
  arma::vec probsrates0(noages, fill::zeros); // 1st vec = total indivs in t
  arma::vec probsrates1(noages, fill::zeros); // 2nd vec = total indivs alive in t+1
  arma::vec probsrates2(noages, fill::zeros); // 3rd vec = total fec
  
  arma::mat tmatrix(noages, noages, fill::zeros); // Main output U matrix
  arma::mat fmatrix(noages, noages, fill::zeros); // Main output F matrix
  
  // This main loop counts individuals going through transitions and calculates
  // survival and fecundity
  arma::uvec data_allalive = find(data_alive2);
  int survsum {0};
  double fecsum {0};
  
  for (int i = 0; i < noages; i++) { 
    
    arma::uvec data_indices = find(data_age == sf_minage(i));
    arma::uvec aget_alive = intersect(data_allalive, data_indices);
    int num_aget_alive = aget_alive.n_elem;
    
    if (num_aget_alive > 0) {
      for (int j = 0; j < num_aget_alive; j++) {
        if (data_alive3(aget_alive(j)) > 0) survsum++;
        fecsum = fecsum + data_usedfec(aget_alive(j));
      }
    }
    
    probsrates0(i) = num_aget_alive;
    if (num_aget_alive > 0) {
      probsrates1(i) = static_cast<double>(survsum) / static_cast<double>(num_aget_alive);
      if (sf_repstatus(i) > 0) probsrates2(i) = fec_mod * fecsum / static_cast<double>(num_aget_alive);
    } else {
      probsrates1(i) = 0;
      probsrates2(i) = 0;
    }
    
    if (i < (noages - 1)) tmatrix(i+1, i) = probsrates1(i);
    fmatrix(0, i) = probsrates2(i);
    
    survsum = 0;
    fecsum = 0;
  }
  
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
}

//' Estimate All Elements of Raw Age-By-Stage Population Projection Matrix
//' 
//' Function \code{.subvertedpatrolgroup()} swiftly calculates matrix
//' transitions in raw ahistorical matrices, and serves as the core workhorse
//' function behind \code{\link{arlefko2}()}.
//' 
//' @name subvertedpatrolgroup
//' 
//' @param sge3 The Allstages data frame developed for \code{rlefko2()} covering
//' stage pairs across times \emph{t}+1 and \emph{t}. Generally termed
//' \code{stageexpansion3}.
//' @param sge2 The data frame covering all stages in time \emph{t}. Generally
//' termed \code{stageexpansion2}.
//' @param MainData The demographic dataset modified to hold \code{usedfec} and
//' \code{usedstage} columns.
//' @param StageFrame The full stageframe for the analysis.
//' @param firstage The first true age to start the matrix with.
//' @param finalage The last true age to estimate.
//' @param cont A logical value indicating whether to lump survival past the
//' last age into a final age transition set on the supermatrix diagonal.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.subvertedpatrolgroup)]]
List subvertedpatrolgroup(DataFrame sge3, DataFrame sge2, DataFrame MainData,
  DataFrame StageFrame, int firstage, int finalage, bool cont) {
  
  arma::vec sge3fec32 = sge3["repentry3"];
  arma::vec sge3rep2 = sge3["rep2n"];
  arma::vec sge3indata32 = sge3["indata"];
  arma::vec sge3ovgivent = sge3["ovgiven_t"];
  arma::vec sge3ovgivenf = sge3["ovgiven_f"];
  arma::vec sge3ovestt = sge3["ovest_t"];
  arma::vec sge3ovestf = sge3["ovest_f"];
  arma::vec sge3ovsurvmult = sge3["ovsurvmult"];
  arma::vec sge3ovfecmult = sge3["ovfecmult"];
  arma::vec sge3index321 = sge3["index321"];
  arma::vec sge3index21 = sge3["index21"];
  arma::vec sge3index2 = sge3["stage2n"];
  arma::vec aliveandequal = sge3["aliveandequal"];
  
  arma::vec sge2rep2 = sge2["rep2"];
  arma::vec sge2fec3 = sge2["fec3"];
  arma::vec sge2index21 = sge2["index21"];
  arma::vec sge2stage2 = sge2["stage2"];
  
  arma::vec dataindex321 = MainData["index321"];
  arma::vec dataindex21 = MainData["index21"];
  arma::vec dataalive3 = MainData["alive3"];
  arma::vec datausedfec2 = MainData["usedfec2"];
  
  int totalages = finalage - firstage + 1;
  
  arma::vec sfsizes = StageFrame["sizebin_center"];
  int nostages = sfsizes.n_elem;
  
  int n = dataindex321.n_elem;
  int no21stages = sge2index21.n_elem; // This includes the dead stage in every age
  int noelems = sge3index321.n_elem;
  unsigned int the_chosen_bun {0};
  
  arma::vec probsrates0(noelems, fill::zeros); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1(noelems, fill::zeros); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2(noelems, fill::zeros); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3(noelems, fill::zeros); // 4th vec = total fec for pair stage
  
  arma::mat stage21fec(no21stages, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  arma::mat tmatrix(((nostages-1) * totalages), ((nostages-1) * totalages), fill::zeros); // Main output U matrix
  arma::mat fmatrix(((nostages-1) * totalages), ((nostages-1) * totalages), fill::zeros); // Main output F matrix
  
  // This main loop counts individuals going through transitions and sums their
  // fecundities, and then adds that info to the 3-trans and 2-trans tables
  for (int i = 0; i < n; i++) { 
    
    // The next line yields sum of all individuals with particular transition
    arma::uvec chosen_index321_vec = find(sge3index321 == dataindex321(i));
    
    if (chosen_index321_vec.n_elem > 0) {
      int chosen_index321 = chosen_index321_vec(0);
      probsrates0(chosen_index321) = probsrates0(chosen_index321) + 1; 
    }
    
    // The next line yields sum of all individuals with particular transition
    arma::uvec chosen_index21_vec = find(sge2index21 == dataindex21(i));
    int chosen_index21 = chosen_index21_vec(0);
    
    stage21fec(chosen_index21, 0) = stage21fec(chosen_index21, 0) + 1; 
    if (dataalive3(i) > 0) {
      stage21fec(chosen_index21, 1) = stage21fec(chosen_index21, 1) + 1;
    }
    
    stage21fec(chosen_index21, 2) = stage21fec(chosen_index21, 2) + datausedfec2(i);
  }
  
  // This next loop populates vectors of individuals according to stage in time t
  for (int i = 0; i < noelems; i++) {
    arma::uvec classy_aks = find(sge2index21 == sge3index21(i));
    

    if (classy_aks.n_elem > 0) {
      the_chosen_bun = classy_aks(0);
      
      probsrates1(i) = stage21fec(the_chosen_bun, 0);
      probsrates2(i) = stage21fec(the_chosen_bun, 1);
      probsrates3(i) = stage21fec(the_chosen_bun, 2);
    }
  }
  
  // Here we populate the main U and F matrices
  for (int elem3 = 0; elem3 < noelems; elem3++) {
    
    if (aliveandequal(elem3) != -1) {
      
      // The next lines DO leave NaNs in the matrices when 0 individuals are summed through in probsrates1
      if (sge3ovsurvmult(elem3) < 0) sge3ovsurvmult(elem3) = 1.0;
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
        probsrates1(elem3);
        
      if (sge3ovfecmult(elem3) < 0) sge3ovfecmult(elem3) = 1.0;
      fmatrix(aliveandequal(elem3)) = sge3fec32(elem3) * sge3rep2(elem3) * 
        probsrates3(elem3) * sge3ovfecmult(elem3) / probsrates1(elem3);
    }
  }
  
  // This section corrects for transitions given in the overwrite table
  arma::uvec ovgiventind = find(sge3ovgivent != -1);
  arma::uvec ovgivenfind = find(sge3ovgivenf != -1);
  int ovgtn = ovgiventind.n_elem;
  int ovgfn = ovgivenfind.n_elem;
  
  if (ovgtn > 0) {
    for (int i = 0; i < ovgtn; i++) {
      int matrixelement2 = aliveandequal(ovgiventind(i));
      
      if (matrixelement2 != -1) tmatrix(matrixelement2) = sge3ovgivent(ovgiventind(i));
    }
  }
  
  if (ovgfn > 0) {
    for (int i = 0; i < ovgfn; i++) {
      int matrixelement2 = aliveandequal(ovgivenfind(i));
      
      if (matrixelement2 != -1) fmatrix(matrixelement2) = sge3ovgivenf(ovgivenfind(i));
    }
  }
  
  // This section replaces transitions with proxies as given in the overwrite table
  arma::uvec ovesttind = find(sge3ovestt != -1);
  arma::uvec ovestfind = find(sge3ovestf != -1);
  int ovestn = ovesttind.n_elem;
  int ovesfn = ovestfind.n_elem;
  
  if (ovestn > 0) {
    for (int i = 0; i < ovestn; i++) {
      arma::uvec replacement = find(sge3index321 == sge3ovestt(ovesttind(i)));
      
      if (aliveandequal(ovesttind(i)) != -1 && aliveandequal(replacement(0)) != -1) {
        tmatrix(aliveandequal(ovesttind(i))) = tmatrix(aliveandequal(replacement(0)));
      }
    }
  }
  
  if (ovesfn > 0) {
    for (int i = 0; i < ovesfn; i++) {
      arma::uvec replacement = find(sge3index321 == sge3ovestf(ovestfind(i)));
      
      if (aliveandequal(ovestfind(i)) != -1 && aliveandequal(replacement(0)) != -1) {
        fmatrix(aliveandequal(ovestfind(i))) = fmatrix(aliveandequal(replacement(0)));
      }
    }
  }
  
  // The next bit changes NAs to 0
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
}

//' Creates Matrices of Year and Patch Terms in Models
//' 
//' Function \code{revelations()} creates a matrix holding either the year or
//' patch coefficients from all vital rate models. This reduces memory load in
//' functions \code{\link{jerzeibalowski}()}, which may be important in some
//' systems or compilers.
//' 
//' @name revelations
//' 
//' @param survproxy The proxy vital rate model covering survival from the main
//' matrix estimator function.
//' @param obsproxy The proxy vital rate model covering observation status from
//' the main matrix estimator function.
//' @param sizeproxy The proxy vital rate model covering primary size from the
//' main matrix estimator function.
//' @param sizebproxy The proxy vital rate model covering secondary size from
//' the main matrix estimator function.
//' @param sizecproxy The proxy vital rate model covering tertiary size from the
//' main matrix estimator function.
//' @param repstproxy The proxy vital rate model covering reproductive status
//' from the main matrix estimator function.
//' @param fecproxy The proxy vital rate model covering fecundity from the main
//' matrix estimator function.
//' @param jsurvproxy The proxy vital rate model covering juvenile survival from
//' the main matrix estimator function.
//' @param jobsproxy The proxy vital rate model covering juvenile observation
//' status from the main matrix estimator function.
//' @param jsizeproxy The proxy vital rate model covering juvenile primary size
//' from the main matrix estimator function.
//' @param jsizebproxy The proxy vital rate model covering juvenile secondary
//' size from the main matrix estimator function.
//' @param jsizecproxy The proxy vital rate model covering juvenile tertiary
//' size from the main matrix estimator function.
//' @param jrepstproxy The proxy vital rate model covering juvenile reproductive
//' status from the main matrix estimator function.
//' @param jmatstproxy The proxy vital rate model covering juvenile probability
//' of becoming mature from the main matrix estimator function.
//' @param mat_switch An integer coding for year (\code{1}) or patch (\code{2}).
//' 
//' @return A matrix with 14 columns corresponding to the number of vital rates
//' and number of columns equal to the number of year or patches.
//' 
//' @keywords internal
//' @noRd
NumericMatrix revelations(List survproxy, List obsproxy, List sizeproxy,
  List sizebproxy, List sizecproxy, List repstproxy, List fecproxy,
  List jsurvproxy, List jobsproxy, List jsizeproxy, List jsizebproxy,
  List jsizecproxy, List jrepstproxy,List jmatstproxy, int mat_switch) {
  
  NumericMatrix final_mat;
  
  if (mat_switch == 1) {
    NumericVector survyear = as<NumericVector>(survproxy["years"]);
    NumericVector obsyear = as<NumericVector>(obsproxy["years"]);
    NumericVector sizeyear = as<NumericVector>(sizeproxy["years"]);
    NumericVector sizebyear = as<NumericVector>(sizebproxy["years"]);
    NumericVector sizecyear = as<NumericVector>(sizecproxy["years"]);
    NumericVector repstyear = as<NumericVector>(repstproxy["years"]);
    NumericVector fecyear = as<NumericVector>(fecproxy["years"]);
    NumericVector jsurvyear = as<NumericVector>(jsurvproxy["years"]);
    NumericVector jobsyear = as<NumericVector>(jobsproxy["years"]);
    NumericVector jsizeyear = as<NumericVector>(jsizeproxy["years"]);
    NumericVector jsizebyear = as<NumericVector>(jsizebproxy["years"]);
    NumericVector jsizecyear = as<NumericVector>(jsizecproxy["years"]);
    NumericVector jrepstyear = as<NumericVector>(jrepstproxy["years"]);
    NumericVector jmatstyear = as<NumericVector>(jmatstproxy["years"]);
    
    int matrows = survyear.length();
    
    NumericMatrix year_mat(matrows, 14);
    year_mat(_, 0) = survyear;
    year_mat(_, 1) = obsyear;
    year_mat(_, 2) = sizeyear;
    year_mat(_, 3) = sizebyear;
    year_mat(_, 4) = sizecyear;
    year_mat(_, 5) = repstyear;
    year_mat(_, 6) = fecyear;
    year_mat(_, 7) = jsurvyear;
    year_mat(_, 8) = jobsyear;
    year_mat(_, 9) = jsizeyear;
    year_mat(_, 10) = jsizebyear;
    year_mat(_, 11) = jsizecyear;
    year_mat(_, 12) = jrepstyear;
    year_mat(_, 13) = jmatstyear;
    
    final_mat = year_mat;
    
  } else if (mat_switch == 2) {
    
    NumericVector survpatch = as<NumericVector>(survproxy["patches"]);
    NumericVector obspatch = as<NumericVector>(obsproxy["patches"]);
    NumericVector sizepatch = as<NumericVector>(sizeproxy["patches"]);
    NumericVector sizebpatch = as<NumericVector>(sizebproxy["patches"]);
    NumericVector sizecpatch = as<NumericVector>(sizecproxy["patches"]);
    NumericVector repstpatch = as<NumericVector>(repstproxy["patches"]);
    NumericVector fecpatch = as<NumericVector>(fecproxy["patches"]);
    NumericVector jsurvpatch = as<NumericVector>(jsurvproxy["patches"]);
    NumericVector jobspatch = as<NumericVector>(jobsproxy["patches"]);
    NumericVector jsizepatch = as<NumericVector>(jsizeproxy["patches"]);
    NumericVector jsizebpatch = as<NumericVector>(jsizebproxy["patches"]);
    NumericVector jsizecpatch = as<NumericVector>(jsizecproxy["patches"]);
    NumericVector jrepstpatch = as<NumericVector>(jrepstproxy["patches"]);
    NumericVector jmatstpatch = as<NumericVector>(jmatstproxy["patches"]);
    
    int matrows = survpatch.length();
    
    NumericMatrix patch_mat(matrows, 14);
    patch_mat(_, 0) = survpatch;
    patch_mat(_, 1) = obspatch;
    patch_mat(_, 2) = sizepatch;
    patch_mat(_, 3) = sizebpatch;
    patch_mat(_, 4) = sizecpatch;
    patch_mat(_, 5) = repstpatch;
    patch_mat(_, 6) = fecpatch;
    patch_mat(_, 7) = jsurvpatch;
    patch_mat(_, 8) = jobspatch;
    patch_mat(_, 9) = jsizepatch;
    patch_mat(_, 10) = jsizebpatch;
    patch_mat(_, 11) = jsizecpatch;
    patch_mat(_, 12) = jrepstpatch;
    patch_mat(_, 13) = jmatstpatch;
    
    final_mat = patch_mat;
  }
  
  return final_mat;
}

//' Creates a Summation of Most Terms Needed in Vital Rate Calculation
//' 
//' Function \code{rimeotam()} provides the majority of the work in creating
//' the linear model sum to be used in vital rate estimation in the MPM. Works
//' specifically with functions \code{\link{jerzeibalowski}()} and
//' \code{\link{motherbalowski}()}.
//' 
//' @name rimeotam
//' 
//' @param maincoefs The coefficients portion of the vital rate model proxy.
//' @param fl1_i Reproductive status in time \emph{t}*-1.
//' @param fl2n_i Reproductive status in time \emph{t}.
//' @param sz1_i Primary size in time \emph{t}-1.
//' @param sz2o_i Primary size in time \emph{t}.
//' @param szb1_i Secondary size in time \emph{t}-1.
//' @param szb2o_i Secondary size in time \emph{t}.
//' @param szc1_i Tertiary size in time \emph{t}-1.
//' @param szc2o_i Tertiary size in time \emph{t}.
//' @param aage2_i Used age in time \emph{t}.
//' @param inda_1 Value of numeric individual covariate a in time \emph{t}-1.
//' @param inda_2 Value of numeric individual covariate a in time \emph{t}.
//' @param indb_1 Value of numeric individual covariate b in time \emph{t}-1.
//' @param indb_2 Value of numeric individual covariate b in time \emph{t}.
//' @param indc_1 Value of numeric individual covariate c in time \emph{t}-1.
//' @param indc_2 Value of numeric individual covariate c in time \emph{t}.
//' @param used_dens Density value used.
//' @param zi A logical value indicating whether model coefficients refer to the
//' zero inflation portion of a model.
//' 
//' @return A single numeric value giving the sum of the products of the linear
//' coefficients and the used status values.
//' 
//' @keywords internal
//' @noRd
double rimeotam(NumericVector maincoefs, double fl1_i, double fl2n_i, double sz1_i,
  double sz2o_i, double szb1_i, double szb2o_i, double szc1_i, double szc2o_i,
  double aage2_i, double inda_1, double inda_2, double indb_1, double indb_2,
  double indc_1, double indc_2, double used_dens, bool zi) {
  
  int add1 {0};
  int add2 {0};
  
  if (zi) {
    add1 = 46;
    add2 = 100;
  }
  
  double parti = maincoefs(0 + add1) + (maincoefs(1 + add1) * fl1_i) + (maincoefs(2 + add1) * fl2n_i) +
    (maincoefs(3 + add1) * sz1_i) + (maincoefs(4 + add1) * sz2o_i) + (maincoefs(5 + add1) * fl2n_i * fl1_i) + 
    (maincoefs(6 + add1) * sz2o_i * sz1_i) + (maincoefs(7 + add1) * sz1_i * fl1_i) +
    (maincoefs(8 + add1) * sz2o_i * fl2n_i) + (maincoefs(9 + add1) * sz2o_i * fl1_i) + 
    (maincoefs(10 + add1) * sz1_i * fl2n_i) + (maincoefs(11 + add1) * aage2_i) + 
    (maincoefs(12 + add1) * aage2_i * sz1_i) + (maincoefs(13 + add1) * aage2_i * sz2o_i) + 
    (maincoefs(14 + add1) * aage2_i * fl1_i) + (maincoefs(15 + add1) * aage2_i * fl2n_i) + 
    (maincoefs(16 + add1) * inda_2) + (maincoefs(17 + add1) * indb_2) + (maincoefs(18 + add1) * indc_2) + 
    (maincoefs(19 + add1) * inda_1) + (maincoefs(20 + add1) * indb_1) + (maincoefs(21 + add1) * indc_1) + 
    (maincoefs(22 + add1) * inda_2 * sz2o_i) + (maincoefs(23 + add1) * indb_2 * sz2o_i) + 
    (maincoefs(24 + add1) * indc_2 * sz2o_i) + (maincoefs(25 + add1) * inda_2 * fl2n_i) + 
    (maincoefs(26 + add1) * indb_2 * fl2n_i) + (maincoefs(27 + add1) * indc_2 * fl2n_i) + 
    (maincoefs(28 + add1) * inda_1 * sz1_i) + (maincoefs(29 + add1) * indb_1 * sz1_i) +
    (maincoefs(30 + add1) * indc_1 * sz1_i) + (maincoefs(31 + add1) * inda_1 * fl1_i) + 
    (maincoefs(32 + add1) * indb_1 * fl1_i) + (maincoefs(33 + add1) * indc_1 * fl1_i) + 
    (maincoefs(34 + add1) * inda_2 * indb_2) + (maincoefs(35 + add1) * inda_2 * indc_2) + 
    (maincoefs(36 + add1) * indb_2 * indc_2) + (maincoefs(37 + add1) * inda_1 * indb_1) + 
    (maincoefs(38 + add1) * inda_1 * indc_1) + (maincoefs(39 + add1) * indb_1 * indc_1) + 
    (maincoefs(40 + add1) * inda_2 * indb_1) + (maincoefs(41 + add1) * inda_1 * indb_2) +
    (maincoefs(42 + add1) * inda_2 * indc_1) + (maincoefs(43 + add1) * inda_1 * indc_2) + 
    (maincoefs(44 + add1) * indb_2 * indc_1) + (maincoefs(45 + add1) * indb_1 * indc_2);
    
  double partii = (maincoefs(100 + add2) * szb2o_i) + (maincoefs(101 + add2) * szb1_i) + 
    (maincoefs(102 + add2) * szc2o_i) + (maincoefs(103 + add2) * szc1_i) + 
    (maincoefs(104 + add2) * used_dens) + (maincoefs(105 + add2) * szb1_i * szb2o_i);
    
  double partiii = (maincoefs(106 + add2) * szc1_i * szc2o_i) + (maincoefs(107 + add2) * sz1_i * szb1_i) + 
    (maincoefs(108 + add2) * sz1_i * szc1_i) + (maincoefs(109 + add2) * szb1_i * szc1_i) + 
    (maincoefs(110 + add2) * sz2o_i * szb2o_i) + (maincoefs(111 + add2) * sz2o_i * szc2o_i) + 
    (maincoefs(112 + add2) * szb2o_i * szc2o_i) + (maincoefs(113 + add2) * sz1_i * szb2o_i) + 
    (maincoefs(114 + add2) * sz1_i * szc2o_i) + (maincoefs(115 + add2) * szb1_i * szc2o_i) + 
    (maincoefs(116 + add2) * sz2o_i * szb1_i) + (maincoefs(117 + add2) * sz2o_i * szc1_i) + 
    (maincoefs(118 + add2) * szb2o_i * szc1_i) + (maincoefs(119 + add2) * sz2o_i * used_dens) + 
    (maincoefs(120 + add2) * szb2o_i * used_dens) + (maincoefs(121 + add2) * szc2o_i * used_dens) + 
    (maincoefs(122 + add2) * sz1_i * used_dens) + (maincoefs(123 + add2) * szb1_i * used_dens) + 
    (maincoefs(124 + add2) * szc1_i * used_dens) + (maincoefs(125 + add2) * fl2n_i * used_dens) + 
    (maincoefs(126 + add2) * fl1_i * used_dens) + (maincoefs(127 + add2) * szb2o_i * fl2n_i) + 
    (maincoefs(128 + add2) * szc2o_i * fl2n_i) + 0 + (maincoefs(130 + add2) * szb1_i * fl1_i) + 
    (maincoefs(131 + add2) * szb2o_i * fl1_i) + (maincoefs(132 + add2) * szb1_i * fl2n_i) + 
    (maincoefs(133 + add2) * szc1_i * fl1_i) + (maincoefs(134 + add2) * szc2o_i * fl1_i) + 
    (maincoefs(135 + add2) * szc1_i * fl2n_i) + (maincoefs(136 + add2) * szb2o_i * aage2_i) + 
    (maincoefs(137 + add2) * szc2o_i * aage2_i) + (maincoefs(138 + add2) * used_dens * aage2_i) + 
    (maincoefs(139 + add2) * szb1_i * aage2_i) + (maincoefs(140 + add2) * szc1_i * aage2_i);
    
  double partiv = (maincoefs(141 + add2) * inda_2 * szb2o_i) + (maincoefs(142 + add2) * inda_2 * szc2o_i) + 
    (maincoefs(143 + add2) * inda_2 * used_dens) + (maincoefs(144 + add2) * inda_1 * szb1_i) + 
    (maincoefs(145 + add2) * inda_1 * szc1_i) + (maincoefs(146 + add2) * inda_1 * szb2o_i) + 
    (maincoefs(147 + add2) * inda_1 * szc2o_i) + (maincoefs(148 + add2) * inda_2 * szb1_i) + 
    (maincoefs(149 + add2) * inda_2 * szc1_i) + (maincoefs(150 + add2) * inda_1 * used_dens);
    
  double partv = (maincoefs(151 + add2) * indb_2 * szb2o_i) + (maincoefs(152 + add2) * indb_2 * szc2o_i) + 
    (maincoefs(153 + add2) * indb_2 * used_dens) + (maincoefs(154 + add2) * indb_1 * szb1_i) + 
    (maincoefs(155 + add2) * indb_1 * szc1_i) + (maincoefs(156 + add2) * indb_1 * szb2o_i) + 
    (maincoefs(157 + add2) * indb_1 * szc2o_i) + (maincoefs(158 + add2) * indb_2 * szb1_i) + 
    (maincoefs(159 + add2) * indb_2 * szc1_i) + (maincoefs(160 + add2) * indb_1 * used_dens);
    
  double partvi = (maincoefs(161 + add2) * indc_2 * szb2o_i) + (maincoefs(162 + add2) * indc_2 * szc2o_i) + 
    (maincoefs(163 + add2) * indc_2 * used_dens) + (maincoefs(164 + add2) * indc_1 * szb1_i) + 
    (maincoefs(165 + add2) * indc_1 * szc1_i) + (maincoefs(166 + add2) * indc_1 * szb2o_i) + 
    (maincoefs(167 + add2) * indc_1 * szc2o_i) + (maincoefs(168 + add2) * indc_2 * szb1_i) + 
    (maincoefs(169 + add2) * indc_2 * szc1_i) + (maincoefs(170 + add2) * indc_1 * used_dens);
    
  double partvii = (maincoefs(171 + add2) * inda_2 * sz1_i) + (maincoefs(172 + add2) * indb_2 * sz1_i) + 
    (maincoefs(173 + add2) * indc_2 * sz1_i) + (maincoefs(174 + add2) * inda_1 * sz2o_i) + 
    (maincoefs(175 + add2) * indb_1 * sz2o_i) + (maincoefs(176 + add2) * indc_1 * sz2o_i) + 
    (maincoefs(177 + add2) * inda_2 * fl1_i) + (maincoefs(178 + add2) * indb_2 * fl1_i) + 
    (maincoefs(179 + add2) * indc_2 * fl1_i) + (maincoefs(180 + add2) * inda_1 * fl2n_i) + 
    (maincoefs(181 + add2) * indb_1 * fl2n_i) + (maincoefs(182 + add2) * indc_1 * fl2n_i);
  
  double albatross = parti + partii + partiii + partiv + partv + partvi + partvii;
  
  return albatross;
}

//' Counts Numbers of Elements in Each Random Individual Covariate Portion of
//' Model
//' 
//' Function \code{foi_counter()} counts the number of elements in each random
//' individual covariate and returns that as a vector.
//' 
//' @name foi_counter
//' 
//' @param modelproxy A list holding the contents of a model processed with
//' function \code{\link{.modelextract}()}
//' @param zi A logical value indicating whether to focus on the zero-inflation
//' parameters.
//' 
//' @return A 6 element vector holding the numbers of elements in each random
//' individual covariate in a model (either the cont portion or the zi portion).
//' 
//' @keywords internal
//' @noRd
arma::ivec foi_counter(List modelproxy, bool zi) {
  
  arma::ivec return_vec(6, fill::zeros);
  
  if (!zi) {
    arma::vec modelinda2r = as<arma::vec>(modelproxy["indcova2s"]);
    arma::vec modelinda1r = as<arma::vec>(modelproxy["indcova1s"]);
    arma::vec modelindb2r = as<arma::vec>(modelproxy["indcovb2s"]);
    arma::vec modelindb1r = as<arma::vec>(modelproxy["indcovb1s"]);
    arma::vec modelindc2r = as<arma::vec>(modelproxy["indcovc2s"]);
    arma::vec modelindc1r = as<arma::vec>(modelproxy["indcovc1s"]);
    
    int v1_l = modelinda2r.n_elem;
    int v2_l = modelinda1r.n_elem;
    int v3_l = modelindb2r.n_elem;
    int v4_l = modelindb1r.n_elem;
    int v5_l = modelindc2r.n_elem;
    int v6_l = modelindc1r.n_elem;
    
    return_vec = {v1_l, v2_l, v3_l, v4_l, v5_l, v6_l};
  } else {
    arma::vec modelinda2r = as<arma::vec>(modelproxy["zeroindcova2s"]);
    arma::vec modelinda1r = as<arma::vec>(modelproxy["zeroindcova1s"]);
    arma::vec modelindb2r = as<arma::vec>(modelproxy["zeroindcovb2s"]);
    arma::vec modelindb1r = as<arma::vec>(modelproxy["zeroindcovb1s"]);
    arma::vec modelindc2r = as<arma::vec>(modelproxy["zeroindcovc2s"]);
    arma::vec modelindc1r = as<arma::vec>(modelproxy["zeroindcovc1s"]);
    
    int v1_l = modelinda2r.n_elem;
    int v2_l = modelinda1r.n_elem;
    int v3_l = modelindb2r.n_elem;
    int v4_l = modelindb1r.n_elem;
    int v5_l = modelindc2r.n_elem;
    int v6_l = modelindc1r.n_elem;
    
    return_vec = {v1_l, v2_l, v3_l, v4_l, v5_l, v6_l};
  }
  
  return return_vec;
}

//' Create Vector of Random Individual Covariate Terms
//' 
//' Function \code{flightoficarus()} creates vectors of random covariate
//' terms.
//' 
//' @name flightoficarus
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector of numeric values for random categorical terms. The order
//' is: 1) cov a time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1,
//' 5) cov c time 2, and 6) cov c time 1. Rows may vary, but must be the same
//' length for each model.
//' 
//' @keywords internal
//' @noRd
NumericVector flightoficarus(List modelproxy) {
  NumericVector modelinda2r = as<NumericVector>(modelproxy["indcova2s"]);
  NumericVector modelinda1r = as<NumericVector>(modelproxy["indcova1s"]);
  NumericVector modelindb2r = as<NumericVector>(modelproxy["indcovb2s"]);
  NumericVector modelindb1r = as<NumericVector>(modelproxy["indcovb1s"]);
  NumericVector modelindc2r = as<NumericVector>(modelproxy["indcovc2s"]);
  NumericVector modelindc1r = as<NumericVector>(modelproxy["indcovc1s"]);
  
  int v1_l = modelinda2r.length();
  int v2_l = modelinda1r.length();
  int v3_l = modelindb2r.length();
  int v4_l = modelindb1r.length();
  int v5_l = modelindc2r.length();
  int v6_l = modelindc1r.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
  
  NumericVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Concatenated Vector of Random Individual Covariate Term Names
//' 
//' Function \code{bootson()} creates a concatenated string vector holding all
//' covariate term names.
//' 
//' @name bootson
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector holding all covariate name terms. The order is: 1) cov a
//' time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1, 5) cov c time 2,
//' and 6) cov c time 1. Note that the element order is the same as in function
//' \code{\link{.flightoficarus}()}.
//' 
//' @keywords internal
//' @noRd
StringVector bootson(List modelproxy) {
  NumericVector modelinda2r_df = as<NumericVector>(modelproxy["indcova2s"]);
  NumericVector modelinda1r_df = as<NumericVector>(modelproxy["indcova1s"]);
  NumericVector modelindb2r_df = as<NumericVector>(modelproxy["indcovb2s"]);
  NumericVector modelindb1r_df = as<NumericVector>(modelproxy["indcovb1s"]);
  NumericVector modelindc2r_df = as<NumericVector>(modelproxy["indcovc2s"]);
  NumericVector modelindc1r_df = as<NumericVector>(modelproxy["indcovc1s"]);
  
  StringVector modelinda2r_rownames = modelinda2r_df.attr("names");
  StringVector modelinda1r_rownames = modelinda1r_df.attr("names");
  StringVector modelindb2r_rownames = modelindb2r_df.attr("names");
  StringVector modelindb1r_rownames = modelindb1r_df.attr("names");
  StringVector modelindc2r_rownames = modelindc2r_df.attr("names");
  StringVector modelindc1r_rownames = modelindc1r_df.attr("names");
  
  int v1_l = modelinda2r_rownames.length();
  int v2_l = modelinda1r_rownames.length();
  int v3_l = modelindb2r_rownames.length();
  int v4_l = modelindb1r_rownames.length();
  int v5_l = modelindc2r_rownames.length();
  int v6_l = modelindc1r_rownames.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;

  StringVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r_rownames(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Vector of Random Individual Covariate Terms for Zero-Inflated Models
//' 
//' Function \code{zero_flightoficarus()} creates vectors of random covariate
//' terms from the binomial portion of a zero-inflated model.
//' 
//' @name zero_flightoficarus
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector of numeric values for random categorical terms. The order
//' is: 1) cov a time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1,
//' 5) cov c time 2, and 6) cov c time 1. Rows may vary, but must be the same
//' length for each model.
//' 
//' @keywords internal
//' @noRd
NumericVector zero_flightoficarus(List modelproxy) {
  NumericVector modelinda2r = as<NumericVector>(modelproxy["zeroindcova2s"]);
  NumericVector modelinda1r = as<NumericVector>(modelproxy["zeroindcova1s"]);
  NumericVector modelindb2r = as<NumericVector>(modelproxy["zeroindcovb2s"]);
  NumericVector modelindb1r = as<NumericVector>(modelproxy["zeroindcovb1s"]);
  NumericVector modelindc2r = as<NumericVector>(modelproxy["zeroindcovc2s"]);
  NumericVector modelindc1r = as<NumericVector>(modelproxy["zeroindcovc1s"]);
  
  int v1_l = modelinda2r.length();
  int v2_l = modelinda1r.length();
  int v3_l = modelindb2r.length();
  int v4_l = modelindb1r.length();
  int v5_l = modelindc2r.length();
  int v6_l = modelindc1r.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
  
  NumericVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Concatenated Vector of Random Individual Covariate Term Names from
//' a Zero-Inflated Model
//' 
//' Function \code{zero_bootson()} creates a concatenated string vector holding
//' all covariate term names from the binomial portion of a zero-inflated model.
//' 
//' @name zero_bootson
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector holding all covariate name terms. The order is: 1) cov a
//' time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1, 5) cov c time 2,
//' and 6) cov c time 1. Note that the element order is the same as in function
//' \code{\link{.zero_flightoficarus}()}.
//' 
//' @keywords internal
//' @noRd
StringVector zero_bootson(List modelproxy) {
  NumericVector modelinda2r_df = as<NumericVector>(modelproxy["zeroindcova2s"]);
  NumericVector modelinda1r_df = as<NumericVector>(modelproxy["zeroindcova1s"]);
  NumericVector modelindb2r_df = as<NumericVector>(modelproxy["zeroindcovb2s"]);
  NumericVector modelindb1r_df = as<NumericVector>(modelproxy["zeroindcovb1s"]);
  NumericVector modelindc2r_df = as<NumericVector>(modelproxy["zeroindcovc2s"]);
  NumericVector modelindc1r_df = as<NumericVector>(modelproxy["zeroindcovc1s"]);
  
  StringVector modelinda2r_rownames = modelinda2r_df.attr("names");
  StringVector modelinda1r_rownames = modelinda1r_df.attr("names");
  StringVector modelindb2r_rownames = modelindb2r_df.attr("names");
  StringVector modelindb1r_rownames = modelindb1r_df.attr("names");
  StringVector modelindc2r_rownames = modelindc2r_df.attr("names");
  StringVector modelindc1r_rownames = modelindc1r_df.attr("names");
  
  int v1_l = modelinda2r_rownames.length();
  int v2_l = modelinda1r_rownames.length();
  int v3_l = modelindb2r_rownames.length();
  int v4_l = modelindb1r_rownames.length();
  int v5_l = modelindc2r_rownames.length();
  int v6_l = modelindc1r_rownames.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;

  StringVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r_rownames(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Index of Element Numbers for Random Individual Covariate Terms
//' 
//' Function \code{foi_index()} creates a matrix indexing the end points of
//' each random individual covariate in the utilized vectors.
//' 
//' @name foi_index
//' 
//' @param surv_proxy Adult survival model proxy.
//' @param obs_proxy Adult observation status model proxy.
//' @param size_proxy Adult primary size model proxy.
//' @param sizeb_proxy Adult secondary size model proxy.
//' @param sizec_proxy Adult tertiary size model proxy.
//' @param repst_proxy Adult reproductive status model proxy.
//' @param fec_proxy Adult fecundity model proxy.
//' @param jsurv_proxy Juvenile survival model proxy.
//' @param jobs_proxy Juvenile observation status model proxy.
//' @param jsize_proxy Juvenile primary size model proxy.
//' @param jsizeb_proxy Juvenile secondary size model proxy.
//' @param jsizec_proxy Juvenile tertiary size model proxy.
//' @param jrepst_proxy Juvenile reproductive status model proxy.
//' @param jmatst_proxy Juvenile maturity status model proxy.
//' 
//' @return An integer matrix with 6 rows and 20 columns. The columns contain
//' the number of elements in each random individual covariate term, with the
//' row order being: 1) cov a t2, 2) cov a t1, 3) cov b t2, 4) cov b t1,
//' 5) cov c t2, and 6) cov c t1.
//' 
//' @keywords internal
//' @noRd
arma::imat foi_index(List surv_proxy, List obs_proxy, List size_proxy, 
  List sizeb_proxy, List sizec_proxy, List repst_proxy, List fec_proxy,
  List jsurv_proxy, List jobs_proxy, List jsize_proxy, List jsizeb_proxy,
  List jsizec_proxy, List jrepst_proxy, List jmatst_proxy) {
  
  arma::ivec surv_fc = foi_counter(surv_proxy, false);
  arma::ivec obs_fc = foi_counter(obs_proxy, false);
  arma::ivec size_fc = foi_counter(size_proxy, false);
  arma::ivec sizeb_fc = foi_counter(sizeb_proxy, false);
  arma::ivec sizec_fc = foi_counter(sizec_proxy, false);
  arma::ivec repst_fc = foi_counter(repst_proxy, false);
  arma::ivec fec_fc = foi_counter(fec_proxy, false);
  arma::ivec jsurv_fc = foi_counter(jsurv_proxy, false);
  arma::ivec jobs_fc = foi_counter(jobs_proxy, false);
  arma::ivec jsize_fc = foi_counter(jsize_proxy, false);
  arma::ivec jsizeb_fc = foi_counter(jsizeb_proxy, false);
  arma::ivec jsizec_fc = foi_counter(jsizec_proxy, false);
  arma::ivec jrepst_fc = foi_counter(jrepst_proxy, false);
  arma::ivec jmatst_fc = foi_counter(jmatst_proxy, false);
  arma::ivec size_fc_zi = foi_counter(size_proxy, true);
  arma::ivec sizeb_fc_zi = foi_counter(sizeb_proxy, true);
  arma::ivec sizec_fc_zi = foi_counter(sizec_proxy, true);
  arma::ivec fec_fc_zi = foi_counter(fec_proxy, true);
  arma::ivec jsize_fc_zi = foi_counter(jsize_proxy, true);
  arma::ivec jsizeb_fc_zi = foi_counter(jsizeb_proxy, true);
  arma::ivec jsizec_fc_zi = foi_counter(jsizec_proxy, true);
  
  arma::imat final_mat(6, 21, fill::zeros);
  
  for (int i = 0; i < 6; i++) {
    final_mat(i, 0) = surv_fc(i);
    final_mat(i, 1) = obs_fc(i);
    final_mat(i, 2) = size_fc(i);
    final_mat(i, 3) = sizeb_fc(i);
    final_mat(i, 4) = sizec_fc(i);
    final_mat(i, 5) = repst_fc(i);
    final_mat(i, 6) = fec_fc(i);
    final_mat(i, 7) = jsurv_fc(i);
    final_mat(i, 8) = jobs_fc(i);
    final_mat(i, 9) = jsize_fc(i);
    final_mat(i, 10) = jsizeb_fc(i);
    final_mat(i, 11) = jsizec_fc(i);
    final_mat(i, 12) = jrepst_fc(i);
    final_mat(i, 13) = jmatst_fc(i);
    final_mat(i, 14) = size_fc_zi(i);
    final_mat(i, 15) = sizeb_fc_zi(i);
    final_mat(i, 16) = sizec_fc_zi(i);
    final_mat(i, 17) = fec_fc_zi(i);
    final_mat(i, 18) = jsize_fc_zi(i);
    final_mat(i, 19) = jsizeb_fc_zi(i);
    final_mat(i, 20) = jsizec_fc_zi(i);
  }
  
  return final_mat;
}

//' Estimate Value for Vital Rate Based on Inputs
//' 
//' Function \code{preouterator()} calculates the value of the vital rate called
//' for by the function \code{jerzeibalowski()}..
//' 
//' @name preouterator
//' 
//' @param modelproxy A model_proxy object derived from function
//' \code{modelextract()}.
//' @param maincoefs The coefficients portion of the vital rate model proxy.
//' @param randindex An integer matrix indexing all random covariates for all
//' vital rates.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' and juvenile reproductive status.
//' @param vitalyear A matrix with year coefficients for all vital rates.
//' @param vitalpatch A matrix with patch coefficients for all vital rates.
//' @param chosen_r2inda A string identifying random covariate a in time t.
//' @param chosen_r1inda A string identifying random covariate a in time t-1.
//' @param chosen_r2indb A string identifying random covariate b in time t.
//' @param chosen_r1indb A string identifying random covariate b in time t-1.
//' @param chosen_r2indc A string identifying random covariate c in time t.
//' @param chosen_r1indc A string identifying random covariate c in time t-1.
//' @param status_terms A NumericVector containing, in order: fl1_i, fl2n_i,
//' sz1_i, sz2o_i, szb1_i, szb2o_i, szc1_i, szc2o_i, aage2_i, inda_1, inda_2,
//' indb_1, indb_2, indc_1, indc_2, used_dens, sz3_i, szb3_i, szc3_i,
//' binwidth3_i, binbwidth3_i, and bincwidth3_i.
//' @param modelgroups2 A vector of group slope coefficients for time t.
//' @param modelgroups1 A vector of group slope coefficients for time t-1.
//' @param modelgroups2zi A vector of zero-inflation model group slope
//' coefficients for time t.
//' @param modelgroups1zi A vector of zero-inflation model group slope
//' coefficients for time t-1.
//' @param modelyearzi A vector of zero-inflation model time slope coefficients.
//' @param modelpatchzi A vector of zero-inflation model patch slope coefficients.
//' @param modelind A vector of individual covariate slope coefficients.
//' @param modelind_rownames A string vector with the names of the individual
//' covariate coefficients.
//' @param modelindzi A vector of individual covariate slope coefficients.
//' @param modelind_rownames_zi A string vector with the names of the individual
//' covariate coefficients.
//' @param zi A logical value indicating whether model coefficients refer to the
//' zero inflation portion of a model.
//' @param sigma The sigma term in the \code{modelproxy} object.
//' @param grp2o_i Stage group number in time \emph{t}.
//' @param grp1_i Stage group number in time \emph{t}-1.
//' @param patchnumber An integer index for pop-patch.
//' @param yearnumber An integer index for monitoring occasion in time \emph{t}.
//' @param vitaldist A parameter specifying the distribution of the vital rate.
//' Current options are: Poisson (0), negative binomial (1), Gaussian (2),
//' Gamma (3), and binomial (4).
//' @param vitalrate An integer specifying the vital rate. 1 = surv, 2 = obs,
//' 3 = size, 4 = sizeb, 5 = sizec, 6 = repst, 7 = fec, 8 = jsurv, 9 = jobs,
//' 10 = jsize, 11 = jsizeb, 12 = jsizec, 13 = jrepst, 14 = jmatst.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param ipm_cdf A logical value indicating whether to use the cumulative
//' density function to estimate size transitions in continuous distributions
//' (\code{true}), or the midpoint method (\code{false}).
//' @param matrixformat An integer representing the style of matrix to develop.
//' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
//' and age-by-stage MPM (4).
//' @param fecmod A scalar multiplier for fecundity.
//' @param repentry_i Rep entry value for time t+1.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param stage2n_i Numeric index of stage in time t.
//' @param nostages The total number of stages in the stageframe.
//' @param modeltrunc An integer coding for zero-truncation status.
//' 
//' @return A class double numeric value for the vital rate being estimated.
//' 
//' @keywords internal
//' @noRd
double preouterator(List modelproxy, NumericVector maincoefs, arma::imat randindex,
  NumericVector dev_terms, NumericMatrix vitalyear, NumericMatrix vitalpatch,
  String chosen_r2inda, String chosen_r1inda, String chosen_r2indb,
  String chosen_r1indb, String chosen_r2indc, String chosen_r1indc,
  NumericVector status_terms, NumericVector modelgroups2,
  NumericVector modelgroups1, NumericVector modelgroups2zi,
  NumericVector modelgroups1zi, NumericVector modelyearzi,
  NumericVector modelpatchzi, NumericVector modelind,
  StringVector modelind_rownames, NumericVector modelindzi,
  StringVector modelind_rownames_zi, bool zi, double sigma, double grp2o_i,
  double grp1_i, int patchnumber, int yearnumber, int vitaldist, int vitalrate,
  double exp_tol, double theta_tol, bool ipm_cdf, int matrixformat,
  double fecmod, double repentry_i, bool negfec, double stage2n_i, int nostages,
  int modeltrunc) {
  
  double preout {0.0};
  double all_out {0.0};
  
  int placeholder = vitalrate - 1;
  int placeholder_zi = placeholder + 12;
  int vitaltype {0}; // Binomial vital rates
  if (vitalrate == 3 || vitalrate == 4 || vitalrate == 5) {
    vitaltype = 1; // Size
  } else if (vitalrate == 10 || vitalrate == 11 || vitalrate == 12) {
    vitaltype = 1; // Juv size
    placeholder_zi = placeholder + 9;
  } else if (vitalrate == 7) {
    vitaltype = 2; // Fecundity
    placeholder_zi = placeholder + 11;
  }
  
  // This section occurs in all vital rates
  double mainsum = rimeotam(maincoefs, status_terms(0), status_terms(1),
    status_terms(2), status_terms(3), status_terms(4), status_terms(5),
    status_terms(6), status_terms(7), status_terms(8), status_terms(9),
    status_terms(10), status_terms(11), status_terms(12), status_terms(13),
    status_terms(14), status_terms(15), zi);
  
  bool zi_processing = false;
  
  if (vitaltype == 1) {
    if (vitalrate == 3 || vitalrate == 10) {
      if (zi && status_terms(16) == 0.0) zi_processing = true;
    } else if (vitalrate == 4 || vitalrate == 11) {
      if (zi && status_terms(17) == 0.0) zi_processing = true;
    } else if (vitalrate == 5 || vitalrate == 12) {
      if (zi && status_terms(18) == 0.0) zi_processing = true;
    } 
  } else if (vitaltype == 2) {
    if (zi && status_terms(16) == 0 && status_terms(17) == 0 &&
      status_terms(18) == 0 && vitaldist < 2) zi_processing = true;  
  }
  
  if (!zi_processing) {
    // Random covariate processing
    double chosen_randcova2 {0.0};
    if (chosen_r2inda != "none") {
      for (int indcount = 0; indcount < randindex(0, placeholder); indcount++) {
        if (chosen_r2inda == modelind_rownames(indcount)) {
          chosen_randcova2 = modelind(indcount);
        }
      }
    }
    double chosen_randcova1 {0.0};
    if (chosen_r1inda != "none") {
      int delectable_sum = randindex(0, placeholder);
      for (int indcount = 0; indcount < randindex(1, placeholder); indcount++) {
        if (chosen_r1inda == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcova1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb2 {0.0};
    if (chosen_r2indb != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder);
      for (int indcount = 0; indcount < randindex(2, placeholder); indcount++) {
        if (chosen_r2indb == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovb2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb1 {0.0};
    if (chosen_r1indb != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder);
      for (int indcount = 0; indcount < randindex(3, placeholder); indcount++) {
        if (chosen_r1indb == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovb1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc2 {0.0};
    if (chosen_r2indc != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder);
      for (int indcount = 0; indcount < randindex(4, placeholder); indcount++) {
        if (chosen_r2indc == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovc2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc1 {0.0};
    if (chosen_r1indc != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder) + randindex(4, placeholder);
      for (int indcount = 0; indcount < randindex(5, placeholder); indcount++) {
        if (chosen_r1indc == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovc1 = modelind(indcount + delectable_sum);
        }
      }
    }
    
    preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
      chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
      chosen_randcovc1 + modelgroups2(grp2o_i) + modelgroups1(grp1_i) + 
      vitalpatch(patchnumber, placeholder) + vitalyear(yearnumber, placeholder) +
      dev_terms(placeholder));
      
    if (preout > exp_tol && vitaldist < 2) preout = exp_tol;
  } else {
    // Only for size and fec
    double chosen_randcova2zi {0.0};
    if (chosen_r2inda != "none") {
      for (int indcount = 0; indcount < randindex(0, placeholder_zi); indcount++) {
        if (chosen_r2inda == modelind_rownames_zi(indcount)) {
          chosen_randcova2zi = modelindzi(indcount);
        }
      }
    }
    double chosen_randcova1zi {0.0};
    if (chosen_r1inda != "none") {
      int delectable_sum = randindex(0, placeholder_zi);
      for (int indcount = 0; indcount < randindex(1, placeholder_zi); indcount++) {
        if (chosen_r1inda == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcova1zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb2zi {0.0};
    if (chosen_r2indb != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi);
      for (int indcount = 0; indcount < randindex(2, placeholder_zi); indcount++) {
        if (chosen_r2indb == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovb2zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb1zi {0.0};
    if (chosen_r1indb != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
        randindex(2, placeholder_zi);
      for (int indcount = 0; indcount < randindex(3, placeholder_zi); indcount++) {
        if (chosen_r1indb == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovb1zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc2zi {0.0};
    if (chosen_r2indc != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
        randindex(2, placeholder_zi) + randindex(3, placeholder_zi);
      for (int indcount = 0; indcount < randindex(4, placeholder_zi); indcount++) {
        if (chosen_r2indc == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovc2zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc1zi {0.0};
    if (chosen_r1indc != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
        randindex(2, placeholder_zi) + randindex(3, placeholder_zi) + randindex(4, placeholder_zi);
      for (int indcount = 0; indcount < randindex(5, placeholder_zi); indcount++) {
        if (chosen_r1indc == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovc1zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    
    preout = (mainsum + chosen_randcova2zi + chosen_randcova1zi +
      chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
      chosen_randcovc1zi + modelgroups2zi(grp2o_i) + modelgroups1zi(grp1_i) + 
      modelpatchzi(patchnumber) + modelyearzi(yearnumber) +
      dev_terms(placeholder));
  }
  
  if (vitaltype == 0) {
    if (preout > exp_tol) preout = exp_tol;
      
    double pre_exp = exp(preout);
    all_out = pre_exp / (1.0 + pre_exp);
    
    // Rcout << "Binomial: preout: " << preout << " pre_exp: " << pre_exp <<
    //   " all_out: " << all_out << "\n";
  } else if (vitaltype == 1) {
    
    double Used_size3 = status_terms(16);
    double Used_binwidth3 = status_terms(19);
    
    if (vitalrate == 4) {
      Used_size3 = status_terms(17);
      Used_binwidth3 = status_terms(20);
      
    } else if (vitalrate == 5) {
      Used_size3 = status_terms(18);
      Used_binwidth3 = status_terms(21);
      
    } else if (vitalrate == 11) {
      Used_size3 = status_terms(17);
      Used_binwidth3 = status_terms(20);
      
    } else if (vitalrate == 12) {
      Used_size3 = status_terms(18);
      Used_binwidth3 = status_terms(21);
    }
    
    if (zi_processing) {
      
      if (preout > exp_tol) preout = exp_tol;
      
      double pre_exp = exp(preout);
      all_out = pre_exp / (1.0 + pre_exp);
      
      // Rcout << "ZI Binomial: preout: " << preout << " pre_exp: " << pre_exp <<
      //   " all_out: " << all_out << "\n";
      
    } else {
      if (vitaldist == 0) {
        // Poisson distribution
        
        if (preout > exp_tol) preout = exp_tol;
        double lambda = exp(preout);
        
        double upper_boundary = (Used_size3 + (Used_binwidth3 / 2));
        double upper_boundary_int = floor(upper_boundary);
        
        double lower_boundary = (Used_size3 - (Used_binwidth3 / 2));
        double lower_boundary_int = floor(lower_boundary);
        
        if (ipm_cdf) {
          if (lower_boundary_int < 0.0) lower_boundary_int = 0.0;
          
          double sizefac {1.0};
          if (upper_boundary_int > 0.0) {
            sizefac = upper_boundary_int * tgamma(upper_boundary_int);
          }
          double main_out = boost::math::tgamma((upper_boundary_int + 1), lambda) / sizefac;
          
          if (upper_boundary_int > lower_boundary_int) {
            double sizefac_low {1.0};
            if (lower_boundary_int > 0.0) {
              sizefac_low = lower_boundary_int * tgamma(lower_boundary_int);
            }
            all_out = main_out - boost::math::tgamma((lower_boundary_int + 1), lambda) / sizefac_low;
          } else {
            all_out = main_out;
          }
          
          if (modeltrunc == 1) {
            double den_corr = (1.0 - (exp(-1 * lambda)));
            all_out = all_out / den_corr;
          }
          // Rcout << "Poisson cdf: upper_boundary_int: " << upper_boundary_int << 
          //   " lower_boundary_int: " << lower_boundary_int << " lambda: " << lambda << 
          //   " all_out: " << all_out << "\n";
          
        } else {
          int y = static_cast<int>(upper_boundary_int);
          int y0 = static_cast<int>(lower_boundary_int);
          if (y0 < -1) y0 = -1;
          
          double current_prob {0.0};
          
          for (int summed_size = (y0 + 1); summed_size <= y; summed_size++) {
            double sizefac {1.0};
            if (Used_size3 > 0.0) {
              sizefac = Used_size3 * tgamma(Used_size3);
            }
            
            double den_corr {1.0};
            if (modeltrunc == 1) den_corr = (1.0 - (exp(-1 * lambda)));
            
            current_prob += ((pow(lambda, Used_size3) * exp(-1.0 * lambda)) / sizefac) / den_corr;
          }
          all_out = current_prob;
          
          // Rcout << "Poisson mid: upper_boundary_int: " << upper_boundary_int <<
          //   " lower_boundary_int: " << lower_boundary_int << " lambda: " << lambda << 
          //   " current_prob: " << current_prob << "\n";
        }
        
      } else if (vitaldist == 1) {
        // Negative binomial
        
        double mu = exp(preout);
        
        double theta = modelproxy["sigma"];
        if (NumericVector::is_na(theta)) theta = 1.0;
        if (theta > theta_tol) theta = theta_tol;
        double alpha = 1.0 / theta;
        
        double upper_boundary = (Used_size3 + (Used_binwidth3 / 2));
        double upper_boundary_int = floor(upper_boundary);
        int y = static_cast<int>(upper_boundary_int);
        
        double lower_boundary = (Used_size3 - (Used_binwidth3 / 2));
        double lower_boundary_int = floor(lower_boundary);
        int y0 = static_cast<int>(lower_boundary_int);
        if (y0 < -1) y0 = -1;
        
        double log_amu = log(alpha) + log(mu);
        double log_mid = -1.0 * theta * log(1.0 + (alpha * mu));
        double den_corr {1.0};
        if (modeltrunc == 1) den_corr = 1.0 - exp(log_mid);
        
        double current_prob {0.0};
        
        for (int summed_size = (y0 + 1); summed_size <= y; summed_size++) {
          double log_leftie = 0.0;
          for (int j = 0; j < summed_size; j++) {
            log_leftie = log(static_cast<double>(j) + theta) - log(static_cast<double>(j) + 1.0) + log_leftie;
          }
          double log_rightie = static_cast<double>(summed_size) * (log_amu - log(1.0 + (alpha * mu)));
          
          double raw_prob = log_leftie + log_mid + log_rightie;
          
          current_prob += exp(raw_prob) / den_corr;
        }
        all_out = current_prob;
        
        // Rcout << "Negbin: y: " << y << " y0: " << y0 << " alpha: " << alpha <<
        //   " mu: " << mu << " current_prob: " << current_prob << "\n";
        
      } else if (vitaldist == 2) {
        // Gaussian size distribution, assuming midpoint
        
        if (ipm_cdf) {
          double lower_size = Used_size3 - (0.5 * Used_binwidth3);
          double upper_size = Used_size3 + (0.5 * Used_binwidth3);
          
          double lower_prob = normcdf(lower_size, preout, sigma);
          double upper_prob = normcdf(upper_size, preout, sigma);
          
          all_out = upper_prob - lower_prob;
          
          // Rcout << "Gaussian cdf: upper_size: " << upper_size << " lower_size: " <<
          //   lower_size << " Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " preout: " << preout << " sigma: " <<
          //   sigma << " upper_prob: " << upper_prob << " lower_prob: " <<
          //   lower_prob << " all_out: " << all_out << "\n";
        } else {
          double sigma2 = sigma * sigma;
          
          all_out = (exp(-1 * (pow((Used_size3 - preout), 2) / (2.0 * sigma2))) / 
            ((pow((2 * M_PI), 0.5)) * sigma));
          all_out = all_out * Used_binwidth3; // This is the midpoint integration
          
          // Rcout << "Gaussian mid: Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " sigma: " << sigma << " preout: " <<
          //   preout << " all_out: " << all_out << "\n";
        }
      } else if (vitaldist == 3) {
        // Gamma size distribution, assuming midpoint
        
        double E_y = 1 / preout;
        double sigma2 = sigma * sigma;
        double alpha = 1.0 / sigma2;
        double beta = (alpha / E_y);
        
        if (ipm_cdf) {
          double lower_size = Used_size3 - (0.5 * Used_binwidth3);
          double upper_size = Used_size3 + (0.5 * Used_binwidth3);
          
          double lower_prob = boost::math::gamma_p(alpha, (beta * lower_size));
          double upper_prob = boost::math::gamma_p(alpha, (beta * upper_size));
          
          all_out = upper_prob - lower_prob;
          
          // Rcout << "Gamma cdf: upper_size: " << upper_size << " lower_size: " <<
          //   lower_size << " alpha: " << alpha << " beta: " << beta << " upper_prob: " <<
          //   upper_prob << " lower_prob: " << lower_prob << " all_out: " << all_out << "\n";
        } else {
          
          all_out = pow(beta, alpha) * (1.0 / tgamma(alpha)) * 
            pow(Used_size3, (alpha - 1.0)) * exp(-1.0 * beta * Used_size3);
          all_out = all_out * Used_binwidth3; // This is the midpoint integration
          
          // Rcout << "Gamma mid: Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " alpha: " << alpha << " beta: " << beta <<
          //   " all_out: " << all_out << "\n";
        }
      }
    }
  } else if (vitaltype == 2) {
    if (matrixformat != 2 || stage2n_i != static_cast<double>(nostages+1)) {
      if (vitaldist == 0 || vitaldist == 1) {
        // Poisson and negative binomial fecundity
        if (preout > exp_tol) preout = exp_tol;
        
        if (zi_processing) {
          
          all_out = (exp(preout) / (1.0 + exp(preout))) * fecmod * repentry_i;
          
        } else {
          
          all_out = exp(preout) * fecmod * repentry_i;
        }
      } else if (vitaldist == 2) {
        // Gaussian fecundity
        all_out = preout * fecmod * repentry_i;
        
        if (negfec && all_out < 0.0) all_out = 0.0;
        
      } else if (vitaldist == 3) {
        // Gamma fecundity
        all_out = (1.0 / preout) * fecmod * repentry_i;
      } else {
        all_out = maincoefs(0);
      }
    } else if (stage2n_i == static_cast<double>(nostages+1)) {
      // This propagates fecundity in deVries-formatted hMPMs
      if (vitaldist == 0 || vitaldist == 1) {
        // Poisson and negative binomial fecundity
        
        if (preout > exp_tol) preout = exp_tol;
            
        if (zi_processing) {
          
          all_out = (exp(preout) / (1.0 + exp(preout))) * fecmod * repentry_i;
          
        } else {
          
            all_out = exp(preout) * fecmod * repentry_i;
          
        }
      } else if (vitaldist == 2) {
        // Gaussian fecundity
        all_out = preout * fecmod * repentry_i;
        
        if (negfec && all_out < 0.0) {
          all_out = 0.0;
        }
      } else if (vitaldist == 3) {
        // Gamma fecundity
        all_out = (1.0 / preout) * fecmod * repentry_i;
      }
    } else {
      all_out = maincoefs(0);
    }
  }
  
  return(all_out);
}

//' Estimate All Elements of Function-based Population Projection Matrix
//' 
//' Function \code{jerzeibalowski()} swiftly calculates matrix elements in
//' function-based population projection matrices. Used in
//' \code{\link{flefko3}()}, \code{\link{flefko2}()}, and
//' \code{\link{aflefko2}()}.
//' 
//' @name jerzeibalowski
//' 
//' @param ppy A data frame showing the population, patch, and year of each
//' matrix to create, in order.
//' @param AllStages A large data frame giving all required inputs for vital
//' rate estimation other than the vital rate model coefficients themselves.
//' Contains a row for each ultimate matrix element.
//' @param stageframe The modified stageframe used in matrix calculations.
//' @param matrixformat An integer representing the style of matrix to develop.
//' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
//' and age-by-stage MPM (4).
//' @param survproxy List of coefficients estimated in model of survival.
//' @param obsproxy List of coefficients estimated in model of observation.
//' @param sizeproxy List of coefficients estimated in model of size.
//' @param repstproxy List of coefficients estimated in model of reproductive 
//' status.
//' @param fecproxy List of coefficients estimated in model of fecundity.
//' @param jsurvproxy List of coefficients estimated in model of juvenile
//' survival.
//' @param jobsproxy List of coefficients estimated in model of juvenile
//' observation.
//' @param jsizeproxy List of coefficients estimated in model of juvenile size.
//' @param jrepstproxy List of coefficients estimated in model of juvenile
//' reproductive status.
//' @param jmatstproxy List of coefficients estimated in model of juvenile
//' maturity probability.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' juvenile reproductive status, and juvenile maturity status.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param maxsize The maximum primary size to be used in element estimation.
//' @param maxsizeb The maximum secondary size to be used in element estimation.
//' @param maxsizec The maximum tertiary size to be used in element estimation.
//' @param firstage The first age to be included in age-by-stage MPM estimation.
//' @param finalage The final age to be included in age-by-stage MPM estimation.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param yearnumber An integer specifying which time at time \emph{t} to
//' develop matrices for. Must be in reference to the \code{listofyears} object
//' developed in the \code{R} matrix estimator function.
//' @param patchnumber An integer specifying which patch to develop matrices
//' for. Must be in reference to the \code{listofyears} object developed in the
//' \code{R} matrix estimator function.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param ipm_method A string indicating which method should be used to
//' estimate size transitions in cases with continuous distributions. Options
//' include \code{"midpoint"}, which uses the midpoint method, and \code{"cdf"},
//' which uses the cumulative density function.
//' @param err_check A logical value indicating whether to export a matrix of
//' conditional probabilities used to develop the \code{U} matrix. Defaults to
//' \code{FALSE}.
//' @param simplicity A logical value indicating whether to output all three
//' matrices (\code{FALSE}), and just matrices \code{U} and \code{F}
//' (\code{TRUE}). Defaults to \code{FALSE}.
//' 
//' @return A list with 2, 3, or 4 elements. If \code{simplicity} is set to
//' \code{FALSE}, then the first 3 elements are matrices, including the main MPM
//' (A), the survival-transition matrix (U), and a fecundity matrix (F). If
//' simplicity is set to \code{TRUE}, then only the survival-transition matrix
//' (U) and fecundity matrix (F) are output. If \code{err_check} is set to
//' \code{TRUE}, then another element is added, which is a 7 column matrix
//' showing survival probability, observation probability, reproduction
//' probability, sizea transition probability, sizeb transition probability,
//' sizec transition probability, and juvenile transition probability to
//' maturity for each element of the final MPM. It is possible that, due to
//' evolving development strategy, further columns are output, as well.
//' 
//' @section Notes:
//' The data frame AllStages introduces variables used in size and fecundity
//' calculations. This DataFrame is broken up into long vectors composed of
//' input sizes and related variables for these calculations. The "model" Lists
//' bring in the vital rate models, and include random coefficients where
//' needed. We also have a number of extra variables, that include such info as
//' whether to use the Poisson, negative binomial, Gamma, or Gaussian
//' distributions for size and fecundity calculations. If \code{sizedist},
//' \code{sizebdist}, \code{sizecdist}, or \code{fecdist} equals 0, 1, 2, or 3,
//' then the Poisson, negative binomial, Gaussian, or Gamma is used,
//' respectively.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.jerzeibalowski)]]
List jerzeibalowski(DataFrame AllStages, DataFrame stageframe, int matrixformat,
  List survproxy, List obsproxy, List sizeproxy, List sizebproxy,
  List sizecproxy, List repstproxy, List fecproxy, List jsurvproxy,
  List jobsproxy, List jsizeproxy, List jsizebproxy, List jsizecproxy,
  List jrepstproxy, List jmatstproxy, NumericVector f2_inda,
  NumericVector f1_inda, NumericVector f2_indb, NumericVector f1_indb,
  NumericVector f2_indc, NumericVector f1_indc, StringVector r2_inda,
  StringVector r1_inda, StringVector r2_indb, StringVector r1_indb,
  StringVector r2_indc, StringVector r1_indc, NumericVector dev_terms,
  double dens, double fecmod, double maxsize, double maxsizeb, double maxsizec,
  unsigned int firstage, unsigned int finalage, bool negfec, int yearnumber,
  int patchnumber, double exp_tol = 700.0, double theta_tol = 100000000.0,
  String ipm_method = "cdf", bool err_check = false, bool simplicity = false) {
  
  NumericMatrix out; // Initialization
  bool ipm_cdf = true;
  if (ipm_method == "midpoint") ipm_cdf = false;
  
  // Determines the size of the matrix
  StringVector stagenames = stageframe["stage"];
  int nostages = stagenames.length();
  unsigned long matrixdim {0};
  
  int nostages_counter = nostages;
  for (int i = 0; i < nostages_counter; i++) {
    if (stringcompare_hard(as<std::string>(stagenames(i)), "AlmostBorn")) nostages -= 1;  
    if (stringcompare_hard(as<std::string>(stagenames(i)), "Dead")) nostages -= 1;
  }
  
  if (matrixformat == 1) { // Ehrlen-format hMPM
    matrixdim = nostages * nostages;
  } else if (matrixformat == 2) { // deVries-format hMPM
    matrixdim = nostages * (nostages + 1);
  } else if (matrixformat == 3) { // ahMPM
    matrixdim = nostages;
  } else if (matrixformat == 4) { // age-by-stage MPM
    matrixdim = nostages * (finalage - firstage + 1);
  }
  
  // Proxy model imports and settings
  bool sizezero = as<bool>(sizeproxy["zero_inflated"]);
  bool sizebzero = as<bool>(sizebproxy["zero_inflated"]);
  bool sizeczero = as<bool>(sizecproxy["zero_inflated"]);
  bool feczero = as<bool>(fecproxy["zero_inflated"]);
  bool jsizezero = as<bool>(jsizeproxy["zero_inflated"]);
  bool jsizebzero = as<bool>(jsizebproxy["zero_inflated"]);
  bool jsizeczero = as<bool>(jsizecproxy["zero_inflated"]);
  
  bool sizetrunc = as<bool>(sizeproxy["zero_truncated"]);
  bool sizebtrunc = as<bool>(sizebproxy["zero_truncated"]);
  bool sizectrunc = as<bool>(sizecproxy["zero_truncated"]);
  bool fectrunc = as<bool>(fecproxy["zero_truncated"]);
  bool jsizetrunc = as<bool>(jsizeproxy["zero_truncated"]);
  bool jsizebtrunc = as<bool>(jsizebproxy["zero_truncated"]);
  bool jsizectrunc = as<bool>(jsizecproxy["zero_truncated"]);
  
  NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
  NumericVector obscoefs = as<NumericVector>(obsproxy["coefficients"]);
  NumericVector sizecoefs = as<NumericVector>(sizeproxy["coefficients"]);
  NumericVector sizebcoefs = as<NumericVector>(sizebproxy["coefficients"]);
  NumericVector sizeccoefs = as<NumericVector>(sizecproxy["coefficients"]);
  NumericVector repstcoefs = as<NumericVector>(repstproxy["coefficients"]);
  NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
  NumericVector jsurvcoefs = as<NumericVector>(jsurvproxy["coefficients"]);
  NumericVector jobscoefs = as<NumericVector>(jobsproxy["coefficients"]);
  NumericVector jsizecoefs = as<NumericVector>(jsizeproxy["coefficients"]);
  NumericVector jsizebcoefs = as<NumericVector>(jsizebproxy["coefficients"]);
  NumericVector jsizeccoefs = as<NumericVector>(jsizecproxy["coefficients"]);
  NumericVector jrepstcoefs = as<NumericVector>(jrepstproxy["coefficients"]);
  NumericVector jmatstcoefs = as<NumericVector>(jmatstproxy["coefficients"]);
  
  double survsigma = as<double>(survproxy["sigma"]);
  double obssigma = as<double>(obsproxy["sigma"]);
  double sizesigma = as<double>(sizeproxy["sigma"]);
  double sizebsigma = as<double>(sizebproxy["sigma"]);
  double sizecsigma = as<double>(sizecproxy["sigma"]);
  double repstsigma = as<double>(repstproxy["sigma"]);
  double fecsigma = as<double>(fecproxy["sigma"]);
  double jsurvsigma = as<double>(jsurvproxy["sigma"]);
  double jobssigma = as<double>(jobsproxy["sigma"]);
  double jsizesigma = as<double>(jsizeproxy["sigma"]);
  double jsizebsigma = as<double>(jsizebproxy["sigma"]);
  double jsizecsigma = as<double>(jsizecproxy["sigma"]);
  double jrepstsigma = as<double>(jrepstproxy["sigma"]);
  double jmatstsigma = as<double>(jmatstproxy["sigma"]);
  
  int survdist = as<int>(survproxy["dist"]);
  int obsdist = as<int>(obsproxy["dist"]);
  int sizedist = as<int>(sizeproxy["dist"]);
  int sizebdist = as<int>(sizebproxy["dist"]);
  int sizecdist = as<int>(sizecproxy["dist"]);
  int repstdist = as<int>(repstproxy["dist"]);
  int fecdist = as<int>(fecproxy["dist"]);
  int jsurvdist = as<int>(jsurvproxy["dist"]);
  int jobsdist = as<int>(jobsproxy["dist"]);
  int jsizedist = as<int>(jsizeproxy["dist"]);
  int jsizebdist = as<int>(jsizebproxy["dist"]);
  int jsizecdist = as<int>(jsizecproxy["dist"]);
  int jrepstdist = as<int>(jrepstproxy["dist"]);
  int jmatstdist = as<int>(jmatstproxy["dist"]);
  
  if (NumericVector::is_na(sizesigma)) {
    if (sizedist == 1) {
      sizesigma = 1.0;
    } else {
      sizesigma = 0.0;
    }
  }
  if (NumericVector::is_na(sizebsigma)) {
    if (sizebdist == 1) {
      sizebsigma = 1.0;
    } else {
      sizebsigma = 0.0;
    }
  }
  if (NumericVector::is_na(sizecsigma)) {
    if (sizecdist == 1) {
      sizecsigma = 1.0;
    } else {
      sizecsigma = 0.0;
    }
  }
  if (NumericVector::is_na(fecsigma)) {
    if (fecdist == 1) {
      fecsigma = 1.0;
    } else {
      fecsigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizesigma)) {
    if (sizedist == 1) {
      jsizesigma = 1.0;
    } else {
      jsizesigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizebsigma)) {
    if (sizebdist == 1) {
      jsizebsigma = 1.0;
    } else {
      jsizebsigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizecsigma)) {
    if (sizecdist == 1) {
      jsizecsigma = 1.0;
    } else {
      jsizecsigma = 0.0;
    }
  }
  
  NumericMatrix vital_year = revelations(survproxy, obsproxy, sizeproxy, sizebproxy,
    sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy, jsizeproxy,
    jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 1);
  
  NumericMatrix vital_patch = revelations(survproxy, obsproxy, sizeproxy,
    sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
    jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 2);
  
  arma::imat rand_index = foi_index(survproxy, obsproxy, sizeproxy, sizebproxy,
    sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy, jsizeproxy,
    jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy);
  
  // NumericVector imports from model_proxy objects
  NumericVector sizeyearzi = as<NumericVector>(sizeproxy["zeroyear"]);
  NumericVector sizebyearzi = as<NumericVector>(sizebproxy["zeroyear"]);
  NumericVector sizecyearzi = as<NumericVector>(sizecproxy["zeroyear"]);
  NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
  NumericVector jsizeyearzi = as<NumericVector>(jsizeproxy["zeroyear"]);
  NumericVector jsizebyearzi = as<NumericVector>(jsizebproxy["zeroyear"]);
  NumericVector jsizecyearzi = as<NumericVector>(jsizecproxy["zeroyear"]);
  
  NumericVector dud_yearzi(sizeyearzi.length());
  
  NumericVector unisyzi = unique(sizeyearzi);
  NumericVector unisyzbi = unique(sizebyearzi);
  NumericVector unisyzci = unique(sizecyearzi);
  NumericVector unijsyzi = unique(jsizeyearzi);
  NumericVector unijsyzbi = unique(jsizebyearzi);
  NumericVector unijsyzci = unique(jsizecyearzi);
  NumericVector unifeci = unique(fecyearzi);
  
  NumericVector sizepatchzi = as<NumericVector>(sizeproxy["zeropatch"]);
  NumericVector sizebpatchzi = as<NumericVector>(sizebproxy["zeropatch"]);
  NumericVector sizecpatchzi = as<NumericVector>(sizecproxy["zeropatch"]);
  NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
  NumericVector jsizepatchzi = as<NumericVector>(jsizeproxy["zeropatch"]);
  NumericVector jsizebpatchzi = as<NumericVector>(jsizebproxy["zeropatch"]);
  NumericVector jsizecpatchzi = as<NumericVector>(jsizecproxy["zeropatch"]);
  
  NumericVector dud_patchzi(sizepatchzi.length());
  
  NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
  NumericVector obsgroups2 = as<NumericVector>(obsproxy["groups2"]);
  NumericVector sizegroups2 = as<NumericVector>(sizeproxy["groups2"]);
  NumericVector sizebgroups2 = as<NumericVector>(sizebproxy["groups2"]);
  NumericVector sizecgroups2 = as<NumericVector>(sizecproxy["groups2"]);
  NumericVector repstgroups2 = as<NumericVector>(repstproxy["groups2"]);
  NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
  NumericVector jsurvgroups2 = as<NumericVector>(jsurvproxy["groups2"]);
  NumericVector jobsgroups2 = as<NumericVector>(jobsproxy["groups2"]);
  NumericVector jsizegroups2 = as<NumericVector>(jsizeproxy["groups2"]);
  NumericVector jsizebgroups2 = as<NumericVector>(jsizebproxy["groups2"]);
  NumericVector jsizecgroups2 = as<NumericVector>(jsizecproxy["groups2"]);
  NumericVector jrepstgroups2 = as<NumericVector>(jrepstproxy["groups2"]);
  NumericVector jmatstgroups2 = as<NumericVector>(jmatstproxy["groups2"]);
  
  NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
  NumericVector obsgroups1 = as<NumericVector>(obsproxy["groups1"]);
  NumericVector sizegroups1 = as<NumericVector>(sizeproxy["groups1"]);
  NumericVector sizebgroups1 = as<NumericVector>(sizebproxy["groups1"]);
  NumericVector sizecgroups1 = as<NumericVector>(sizecproxy["groups1"]);
  NumericVector repstgroups1 = as<NumericVector>(repstproxy["groups1"]);
  NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
  NumericVector jsurvgroups1 = as<NumericVector>(jsurvproxy["groups1"]);
  NumericVector jobsgroups1 = as<NumericVector>(jobsproxy["groups1"]);
  NumericVector jsizegroups1 = as<NumericVector>(jsizeproxy["groups1"]);
  NumericVector jsizebgroups1 = as<NumericVector>(jsizebproxy["groups1"]);
  NumericVector jsizecgroups1 = as<NumericVector>(jsizecproxy["groups1"]);
  NumericVector jrepstgroups1 = as<NumericVector>(jrepstproxy["groups1"]);
  NumericVector jmatstgroups1 = as<NumericVector>(jmatstproxy["groups1"]);
  
  NumericVector sizegroups2zi = as<NumericVector>(sizeproxy["zerogroups2"]);
  NumericVector sizebgroups2zi = as<NumericVector>(sizebproxy["zerogroups2"]);
  NumericVector sizecgroups2zi = as<NumericVector>(sizecproxy["zerogroups2"]);
  NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
  NumericVector jsizegroups2zi = as<NumericVector>(jsizeproxy["zerogroups2"]);
  NumericVector jsizebgroups2zi = as<NumericVector>(jsizebproxy["zerogroups2"]);
  NumericVector jsizecgroups2zi = as<NumericVector>(jsizecproxy["zerogroups2"]);
  
  NumericVector dud_groups2zi(jsizecyearzi.length());
  
  NumericVector sizegroups1zi = as<NumericVector>(sizeproxy["zerogroups1"]);
  NumericVector sizebgroups1zi = as<NumericVector>(sizebproxy["zerogroups1"]);
  NumericVector sizecgroups1zi = as<NumericVector>(sizecproxy["zerogroups1"]);
  NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
  NumericVector jsizegroups1zi = as<NumericVector>(jsizeproxy["zerogroups1"]);
  NumericVector jsizebgroups1zi = as<NumericVector>(jsizebproxy["zerogroups1"]);
  NumericVector jsizecgroups1zi = as<NumericVector>(jsizecproxy["zerogroups1"]);
  
  NumericVector dud_groups1zi(jsizecyearzi.length());
  
  NumericVector survind = flightoficarus(survproxy);
  NumericVector obsind = flightoficarus(obsproxy);
  NumericVector sizeind = flightoficarus(sizeproxy);
  NumericVector sizebind = flightoficarus(sizebproxy);
  NumericVector sizecind = flightoficarus(sizecproxy);
  NumericVector repstind = flightoficarus(repstproxy);
  NumericVector fecind = flightoficarus(fecproxy);
  NumericVector jsurvind = flightoficarus(jsurvproxy);
  NumericVector jobsind = flightoficarus(jobsproxy);
  NumericVector jsizeind = flightoficarus(jsizeproxy);
  NumericVector jsizebind = flightoficarus(jsizebproxy);
  NumericVector jsizecind = flightoficarus(jsizecproxy);
  NumericVector jrepstind = flightoficarus(jrepstproxy);
  NumericVector jmatstind = flightoficarus(jmatstproxy);
  
  NumericVector sizeindzi = zero_flightoficarus(sizeproxy);
  NumericVector sizebindzi = zero_flightoficarus(sizebproxy);
  NumericVector sizecindzi = zero_flightoficarus(sizecproxy);
  NumericVector fecindzi = zero_flightoficarus(fecproxy);
  NumericVector jsizeindzi = zero_flightoficarus(jsizeproxy);
  NumericVector jsizebindzi = zero_flightoficarus(jsizebproxy);
  NumericVector jsizecindzi = zero_flightoficarus(jsizecproxy);
  
  StringVector survind_rownames = bootson(survproxy);
  StringVector obsind_rownames = bootson(obsproxy);
  StringVector sizeind_rownames = bootson(sizeproxy);
  StringVector sizebind_rownames = bootson(sizebproxy);
  StringVector sizecind_rownames = bootson(sizecproxy);
  StringVector repstind_rownames = bootson(repstproxy);
  StringVector fecind_rownames = bootson(fecproxy);
  StringVector jsurvind_rownames = bootson(jsurvproxy);
  StringVector jobsind_rownames = bootson(jobsproxy);
  StringVector jsizeind_rownames = bootson(jsizeproxy);
  StringVector jsizebind_rownames = bootson(jsizebproxy);
  StringVector jsizecind_rownames = bootson(jsizecproxy);
  StringVector jrepstind_rownames = bootson(jrepstproxy);
  StringVector jmatstind_rownames = bootson(jmatstproxy);
  
  StringVector sizeind_rownames_zi = zero_bootson(sizeproxy);
  StringVector sizebind_rownames_zi = zero_bootson(sizebproxy);
  StringVector sizecind_rownames_zi = zero_bootson(sizecproxy);
  StringVector fecind_rownames_zi = zero_bootson(fecproxy);
  StringVector jsizeind_rownames_zi = zero_bootson(jsizeproxy);
  StringVector jsizebind_rownames_zi = zero_bootson(jsizebproxy);
  StringVector jsizecind_rownames_zi = zero_bootson(jsizecproxy);
  
  // AllStages import and settings
  Rcpp::NumericVector stage3_num = as<NumericVector>(AllStages["stage3"]);
  Rcpp::NumericVector stage2n_num = as<NumericVector>(AllStages["stage2n"]);
  Rcpp::NumericVector stage2o_num = as<NumericVector>(AllStages["stage2o"]);
  arma::vec stage3 = as<arma::vec>(stage3_num);
  arma::vec stage2n = as<arma::vec>(stage2n_num);
  arma::vec stage2o = as<arma::vec>(stage2o_num);
  
  Rcpp::NumericVector sz3 = as<NumericVector>(AllStages["size3"]);
  Rcpp::NumericVector sz2n = as<NumericVector>(AllStages["size2n"]);
  Rcpp::NumericVector sz2o = as<NumericVector>(AllStages["size2o"]);
  Rcpp::NumericVector sz1 = as<NumericVector>(AllStages["size1"]);
  Rcpp::NumericVector szb3 = as<NumericVector>(AllStages["sizeb3"]);
  Rcpp::NumericVector szb2n = as<NumericVector>(AllStages["sizeb2n"]);
  Rcpp::NumericVector szb2o = as<NumericVector>(AllStages["sizeb2o"]);
  Rcpp::NumericVector szb1 = as<NumericVector>(AllStages["sizeb1"]);
  Rcpp::NumericVector szc3 = as<NumericVector>(AllStages["sizec3"]);
  Rcpp::NumericVector szc2n = as<NumericVector>(AllStages["sizec2n"]);
  Rcpp::NumericVector szc2o = as<NumericVector>(AllStages["sizec2o"]);
  Rcpp::NumericVector szc1 = as<NumericVector>(AllStages["sizec1"]);
  Rcpp::NumericVector ob3 = as<NumericVector>(AllStages["obs3"]);
  Rcpp::NumericVector fl3 = as<NumericVector>(AllStages["rep3"]);
  Rcpp::NumericVector fl2n = as<NumericVector>(AllStages["rep2n"]);
  Rcpp::NumericVector fl2o = as<NumericVector>(AllStages["rep2o"]);
  Rcpp::NumericVector fl1 = as<NumericVector>(AllStages["rep1"]);
  Rcpp::NumericVector mat3 = as<NumericVector>(AllStages["mat3"]);
  Rcpp::NumericVector mat2n = as<NumericVector>(AllStages["mat2n"]);
  Rcpp::NumericVector mat2o = as<NumericVector>(AllStages["mat2o"]);
  Rcpp::NumericVector mat1 = as<NumericVector>(AllStages["mat1"]);
  Rcpp::NumericVector immat2n = as<NumericVector>(AllStages["imm2n"]);
  Rcpp::NumericVector immat2o = as<NumericVector>(AllStages["imm2o"]);
  Rcpp::NumericVector immat1 = as<NumericVector>(AllStages["imm1"]);
  
  Rcpp::NumericVector repentry = as<NumericVector>(AllStages["repentry3"]);
  Rcpp::NumericVector indata2n = as<NumericVector>(AllStages["indata2n"]);
  Rcpp::NumericVector indata2o = as<NumericVector>(AllStages["indata2o"]);
  Rcpp::NumericVector binwidth3 = as<NumericVector>(AllStages["binwidth"]);
  Rcpp::NumericVector binbwidth3 = as<NumericVector>(AllStages["binbwidth"]);
  Rcpp::NumericVector bincwidth3 = as<NumericVector>(AllStages["bincwidth"]);
  Rcpp::NumericVector actualage2 = as<NumericVector>(AllStages["actualage"]);
  
  Rcpp::NumericVector grp3 = as<NumericVector>(AllStages["group3"]);
  Rcpp::NumericVector grp2n = as<NumericVector>(AllStages["group2n"]);
  Rcpp::NumericVector grp2o = as<NumericVector>(AllStages["group2o"]);
  Rcpp::NumericVector grp1 = as<NumericVector>(AllStages["group1"]);
  
  Rcpp::NumericVector ovestt_num = as<NumericVector>(AllStages["ovest_t"]);
  arma::vec ovestt = as<arma::vec>(ovestt_num);
  
  Rcpp::NumericVector ovestf_num = as<NumericVector>(AllStages["ovest_f"]);
  arma::vec ovestf = as<arma::vec>(ovestf_num);
  
  Rcpp::NumericVector indata = as<NumericVector>(AllStages["indata"]);
  Rcpp::NumericVector ovgivent = as<NumericVector>(AllStages["ovgiven_t"]);
  Rcpp::NumericVector ovgivenf = as<NumericVector>(AllStages["ovgiven_f"]);
  
  Rcpp::NumericVector ovsurvmult = as<NumericVector>(AllStages["ovsurvmult"]);
  Rcpp::NumericVector ovfecmult = as<NumericVector>(AllStages["ovfecmult"]);
  
  Rcpp::IntegerVector index321_int = as<IntegerVector>(AllStages["index321"]);
  arma::uvec index321 = as<arma::uvec>(index321_int);
  
  Rcpp::IntegerVector aliveandequal = as<IntegerVector>(AllStages["aliveandequal"]); // Used to be NumericVector
  
  int n = stage3.n_elem;
  
  arma::uvec replacetvec = find(ovestt != -1.0);
  arma::uvec replacefvec = find(ovestf != -1.0);
  int replacementst = replacetvec.n_elem;
  int replacementsf = replacefvec.n_elem;
  int repindex {0};
  int properindex {0};
  int proxyindex {0};
  
  // Determination of choices of fixed and random individual covariates
  double inda1 = f1_inda(yearnumber);
  double indb1 = f1_indb(yearnumber);
  double indc1 = f1_indc(yearnumber);
  double inda2 = f2_inda(yearnumber);
  double indb2 = f2_indb(yearnumber);
  double indc2 = f2_indc(yearnumber);
  
  String chosen_r2inda = r2_inda(yearnumber);
  String chosen_r1inda = r1_inda(yearnumber);
  String chosen_r2indb = r2_indb(yearnumber);
  String chosen_r1indb = r1_indb(yearnumber);
  String chosen_r2indc = r2_indc(yearnumber);
  String chosen_r1indc = r1_indc(yearnumber);
  
  // The output matrix to collect conditional probabilities
  // Matrix out is 0 matrix with n rows & 6 columns: 0 surv, 1 obs, 2 repst,
  // 3 size, 4 size_b, 5 size_c, 6 matst, >6 are test variables
  if (err_check) {
    NumericMatrix zeroform(n, 7);
    out = zeroform;
    CharacterVector out_names = {"surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst"};
    colnames(out) = out_names;
  }
  NumericVector out_vec(7);
  
  arma::mat survtransmat(matrixdim, matrixdim, fill::zeros);
  arma::mat fectransmat(matrixdim, matrixdim, fill::zeros);
  
  double fec_addedcoefs = sum(feccoefs);
  double jsurv_coefsadded = sum(jsurvcoefs);
  double mat_predicted {0.0};
  unsigned int k {0};
  // The following loop runs through each line of AllStages, and so runs through
  // each estimable element in the matrix
  for(int i = 0; i < n; i++) {
    out_vec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    k = aliveandequal(i);
    mat_predicted = 0.0;
    
    if (err_check) out(i, 6) = 1.0; // Initialization of maturity status probability for typical case
    
    Rcpp::NumericVector statusterms = {fl1(i), fl2n(i), sz1(i), sz2o(i),
      szb1(i), szb2o(i), szc1(i), szc2o(i), actualage2(i), inda1, inda2, indb1,
      indb2, indc1, indc2, dens, sz3(i), szb3(i), szc3(i), binwidth3(i),
      binbwidth3(i), bincwidth3(i)};
    
    if (ovgivent(i) == -1 && indata(i) == 1 && stage2n(i) == stage2o(i)) {
      if ((mat2n(i) == 1 && mat3(i) == 1) || (mat2o(i) == 1 && mat3(i) == 1)) {
        
        // Adult survival transitions
        if (survdist < 5) {
          out_vec(0) = preouterator(survproxy, survcoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, survgroups2,
            survgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            survind, survind_rownames, sizeindzi, sizeind_rownames_zi, false,
            survsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 1, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
        } else {
          out_vec(0) = survcoefs(0);
        }
        if (err_check) out(i, 0) = out_vec(0);
        
        if (obsdist < 5) {
          out_vec(1) = preouterator(obsproxy, obscoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, obsgroups2,
            obsgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            obsind, obsind_rownames, sizeindzi, sizeind_rownames_zi, false,
            obssigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 2, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
          
        } else {
          out_vec(1) = obscoefs(0);
        }
        if (err_check) out(i, 1) = out_vec(1);
        
        if (ob3(i) == 1 || obsdist == 5) {
          
          if (sizedist < 5) {
            bool used_sizezero = false;
            if (sizezero && sz3(i) == 0) used_sizezero = sizezero;
            
            out_vec(3) = preouterator(sizeproxy, sizecoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, sizegroups2, sizegroups1, sizegroups2zi, sizegroups1zi,
              sizeyearzi, sizepatchzi, sizeind, sizeind_rownames, sizeindzi,
              sizeind_rownames_zi, used_sizezero, sizesigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, sizedist, 3, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
              sizetrunc);
              
          } else {
            out_vec(3) = 1.0;
          }
          if (err_check) out(i, 3) = out_vec(3);
          
          if (sizebdist < 5) {
            bool used_sizebzero = false;
            if (sizebzero && szb3(i) == 0) used_sizebzero = sizebzero;
            
            out_vec(4) = preouterator(sizebproxy, sizebcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, sizebgroups2, sizebgroups1, sizebgroups2zi,
              sizebgroups1zi, sizebyearzi, sizebpatchzi, sizebind,
              sizebind_rownames, sizebindzi, sizebind_rownames_zi, used_sizebzero,
              sizebsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist,
              4, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, sizebtrunc);
          } else {
            out_vec(4) = 1.0;
          }
          if (err_check) out(i, 4) = out_vec(4);
          
          if (sizecdist < 5) {
            bool used_sizeczero = false;
            if (sizeczero && szc3(i) == 0) used_sizeczero = sizeczero;
            
            out_vec(5) = preouterator(sizecproxy, sizeccoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, sizecgroups2, sizecgroups1, sizecgroups2zi,
              sizecgroups1zi, sizecyearzi, sizecpatchzi, sizecind,
              sizecind_rownames, sizecindzi, sizecind_rownames_zi, used_sizeczero,
              sizecsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist,
              5, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, sizectrunc);
          } else {
            out_vec(5) = 1.0;
          }
          if (err_check) out(i, 5) = out_vec(5);
          
          if (repstdist < 5) {
            out_vec(2) = preouterator(repstproxy, repstcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, repstgroups2, repstgroups1, dud_groups2zi, dud_groups1zi,
              dud_yearzi, dud_patchzi, repstind, repstind_rownames, sizeindzi,
              sizeind_rownames_zi, false, repstsigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, 4, 6, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages, 0);
              
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - out_vec(2);
            }
          } else {
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - repstcoefs(0);
            } else if (fl3(i) == 1) {
              out_vec(2) = repstcoefs(0);
            } else {
              out_vec(2) = 0.0;
            }
          }
          if (err_check) out(i, 2) = out_vec(2);
          
        } else {
          out_vec(1) = 1.0 - out_vec(1);
          out_vec(2) = 1.0;
          out_vec(3) = 1.0;
          out_vec(4) = 1.0;
          out_vec(5) = 1.0;
          out_vec(6) = 1.0;
          
          if (err_check) {
            out(i, 1) = out_vec(1);
            out(i, 2) = out_vec(2);
            out(i, 3) = out_vec(3);
            out(i, 4) = out_vec(4);
            out(i, 5) = out_vec(5);
            out(i, 6) = out_vec(6);
          }
        }
        survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
          out_vec(4) * out_vec(5) * out_vec(6);
        
      } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurv_coefsadded != 0.0) {
        // Juvenile to adult transitions
        
        if (jmatstdist < 5) {
          mat_predicted = preouterator(jmatstproxy, jmatstcoefs, rand_index,
            dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
            chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms,
            jmatstgroups2, jmatstgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi,
            dud_patchzi, jmatstind, jmatstind_rownames, jsizeindzi,
            jsizeind_rownames_zi, false, jmatstsigma, grp2o(i), grp1(i),
            patchnumber, yearnumber, 4, 21, exp_tol, theta_tol, ipm_cdf, matrixformat,
            fecmod, repentry(i), negfec, stage2n(i), nostages, 0);
          
          if (mat3(i) > 0.5) {
            out_vec(6) = mat_predicted;
          } else {
            out_vec(6) = 1 - mat_predicted;
          }
        } else {
          if (mat3(i) > 0.5) {
            out_vec(6) = 1;
          } else {
            out_vec(6) = 0;
          }
        }
        if (err_check) out(i, 6) = out_vec(6);
        
        if (jsurvdist < 5) {
          out_vec(0) = preouterator(jsurvproxy, jsurvcoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsurvgroups2,
            jsurvgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            jsurvind, jsurvind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
            jsurvsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 8, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
        } else {
          out_vec(0) = jsurvcoefs(0);
        }
        if (err_check) out(i, 0) = out_vec(0);
        
        if (jobsdist < 5) {
          out_vec(1) = preouterator(jobsproxy, jobscoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jobsgroups2,
            jobsgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            jobsind, jobsind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
            jobssigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 9, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
        } else {
          out_vec(1) = jobscoefs(0);
        }
        if (err_check) out(i, 1) = out_vec(1);
        
        if (ob3(i) == 1 || jobsdist == 5) {
          if (jsizedist < 5) {
            out_vec(3) = preouterator(jsizeproxy, jsizecoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsizegroups2,
              jsizegroups1, jsizegroups2zi, jsizegroups1zi, jsizeyearzi, jsizepatchzi,
              jsizeind, jsizeind_rownames, jsizeindzi, jsizeind_rownames_zi, jsizezero,
              jsizesigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizedist, 10,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, jsizetrunc);
          } else {
            out_vec(3) = 1.0;
          }
          if (err_check) out(i, 3) = out_vec(3);
          
          if (jsizebdist < 5) {
            out_vec(4) = preouterator(jsizebproxy, jsizebcoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsizebgroups2,
              jsizebgroups1, jsizebgroups2zi, jsizebgroups1zi, jsizebyearzi, jsizebpatchzi,
              jsizebind, jsizebind_rownames, jsizebindzi, jsizebind_rownames_zi, jsizebzero,
              jsizebsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist, 11,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, jsizebtrunc);
          } else {
            out_vec(4) = 1.0;
          }
          if (err_check) out(i, 4) = out_vec(4);
          
          if (jsizecdist < 5) {
            out_vec(5) = preouterator(jsizecproxy, jsizeccoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda,chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsizecgroups2,
              jsizecgroups1, jsizecgroups2zi, jsizecgroups1zi, jsizecyearzi, jsizecpatchzi,
              jsizecind, jsizecind_rownames, jsizecindzi, jsizecind_rownames_zi, jsizeczero,
              jsizecsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist, 12,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, jsizectrunc);
          } else {
            out_vec(5) = 1.0;
          }
          if (err_check) out(i, 5) = out_vec(5);
          
          if (jrepstdist < 5) {
            out_vec(2) = preouterator(jrepstproxy, jrepstcoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jrepstgroups2,
              jrepstgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
              jrepstind, jrepstind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              jrepstsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 13, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0);
              
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - out_vec(2);
            }
          } else {
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - jrepstcoefs(0);
            } else if (fl3(i) == 1) {
              out_vec(2) = jrepstcoefs(0);
            } else {
              out_vec(2) = 0.0;
            }
          }
          if (err_check) out(i, 2) = out_vec(2);
          
        } else {
          out_vec(1) = 1.0 - out_vec(1);
          out_vec(2) = 1.0;
          out_vec(3) = 1.0;
          out_vec(4) = 1.0;
          out_vec(5) = 1.0;
          out_vec(6) = 1.0;
          
          if (err_check) {
            out(i, 1) = out_vec(1);
            out(i, 2) = out_vec(2);
            out(i, 3) = out_vec(3);
            out(i, 4) = out_vec(4);
            out(i, 5) = out_vec(5);
            out(i, 6) = out_vec(6);
          }
        }
        
        survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
          out_vec(4) * out_vec(5) * out_vec(6);
      }
    } else if (ovgivent(i) != -1) {
      // All other transitions
      
      survtransmat(k) = ovgivent(i);
    }
    
    // This next block calculates fecundity
    if (indata2n(i) == 1 && fec_addedcoefs != 0.0) {
      if (fl2o(i) > 0.0 && ovgivenf(i) == -1.0) {
        
        fectransmat(k) = preouterator(fecproxy, feccoefs, rand_index, dev_terms,
          vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
          chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, fecgroups2,
          fecgroups1, fecgroups2zi, fecgroups1zi, fecyearzi, fecpatchzi, fecind,
          fecind_rownames, fecindzi, fecind_rownames_zi, feczero, fecsigma,
          grp2o(i), grp1(i), patchnumber, yearnumber, fecdist, 7, exp_tol,
          theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
          stage2n(i), nostages, fectrunc);
        
      } else if (ovgivenf(i) != -1 ) {
        fectransmat(k) = ovgivenf(i);
      }
    } else if (ovgivenf(i) != -1 ) {
      fectransmat(k) = ovgivenf(i);
    }
  }
  
  double ov_mult {0};
  if (replacementst > 0) {
    for (int i = 0; i < replacementst; i++) {
      
      repindex = replacetvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestt(repindex));
      
      if (rightindex.n_elem > 0) {
        proxyindex = aliveandequal(rightindex(0));
        
        ov_mult = ovsurvmult(i);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        survtransmat(properindex) = survtransmat(proxyindex) * ov_mult;
      }
    }
  }
  
  if (replacementsf > 0) {
    for (int i = 0; i < replacementsf; i++) {
      
      repindex = replacefvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestf(repindex));
      
      if (rightindex.n_elem > 0) {
        proxyindex = aliveandequal(rightindex(0));
        
        ov_mult = ovfecmult(i);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        fectransmat(properindex) = fectransmat(proxyindex) * ov_mult;
      }
    }
  }
  
  // Final output
  List output(4);
  
  if (!simplicity) {
    arma::mat amatrix = survtransmat + fectransmat;
    output(0) = amatrix;
  } else {
    output(0) = R_NilValue;
  }
  
  output(1) = survtransmat;
  output(2) = fectransmat;
  
  if (err_check) {
    output(3) = out;
  } else {
    output(3) = R_NilValue;
  }
  CharacterVector output_names = {"A", "U", "F", "out"};
  output.attr("names") = output_names;
  
  return output;
}

//' Estimate All Elements of Function-based Population Projection Sparse Matrix
//' 
//' Function \code{jerzeibalowsk_sp()} swiftly calculates matrix elements in
//' function-based population projection matrices. Used in
//' \code{f_projection3()}.
//' 
//' @name jerzeibalowsk_sp
//' 
//' @param ppy A data frame showing the population, patch, and year of each
//' matrix to create, in order.
//' @param AllStages A large data frame giving all required inputs for vital
//' rate estimation other than the vital rate model coefficients themselves.
//' Contains a row for each ultimate matrix element.
//' @param stageframe The modified stageframe used in matrix calculations.
//' @param matrixformat An integer representing the style of matrix to develop.
//' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
//' and age-by-stage MPM (4).
//' @param survproxy List of coefficients estimated in model of survival.
//' @param obsproxy List of coefficients estimated in model of observation.
//' @param sizeproxy List of coefficients estimated in model of size.
//' @param repstproxy List of coefficients estimated in model of reproductive 
//' status.
//' @param fecproxy List of coefficients estimated in model of fecundity.
//' @param jsurvproxy List of coefficients estimated in model of juvenile
//' survival.
//' @param jobsproxy List of coefficients estimated in model of juvenile
//' observation.
//' @param jsizeproxy List of coefficients estimated in model of juvenile size.
//' @param jrepstproxy List of coefficients estimated in model of juvenile
//' reproductive status.
//' @param jmatstproxy List of coefficients estimated in model of juvenile
//' maturity probability.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' juvenile reproductive status, and juvenile maturity status.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param maxsize The maximum primary size to be used in element estimation.
//' @param maxsizeb The maximum secondary size to be used in element estimation.
//' @param maxsizec The maximum tertiary size to be used in element estimation.
//' @param firstage The first age to be included in age-by-stage MPM estimation.
//' @param finalage The final age to be included in age-by-stage MPM estimation.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param yearnumber An integer specifying which time at time \emph{t} to
//' develop matrices for. Must be in reference to the \code{listofyears} object
//' developed in the \code{R} matrix estimator function.
//' @param patchnumber An integer specifying which patch to develop matrices
//' for. Must be in reference to the \code{listofyears} object developed in the
//' \code{R} matrix estimator function.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param ipm_method A string indicating which method should be used to
//' estimate size transitions in cases with continuous distributions. Options
//' include \code{"midpoint"}, which uses the midpoint method, and \code{"cdf"},
//' which uses the cumulative density function.
//' @param err_check A logical value indicating whether to export a matrix of
//' conditional probabilities used to develop the \code{U} matrix. Defaults to
//' \code{FALSE}.
//' @param simplicity A logical value indicating whether to output all three
//' matrices (\code{FALSE}), and just matrices \code{U} and \code{F}
//' (\code{TRUE}). Defaults to \code{FALSE}.
//' 
//' @return A list with 2, 3, or 4 elements. If \code{simplicity} is set to
//' \code{FALSE}, then the first 3 elements are matrices, including the main MPM
//' (A), the survival-transition matrix (U), and a fecundity matrix (F). If
//' simplicity is set to \code{TRUE}, then only the survival-transition matrix
//' (U) and fecundity matrix (F) are output. If \code{err_check} is set to
//' \code{TRUE}, then another element is added, which is a 7 column matrix
//' showing survival probability, observation probability, reproduction
//' probability, sizea transition probability, sizeb transition probability,
//' sizec transition probability, and juvenile transition probability to
//' maturity for each element of the final MPM. It is possible that, due to
//' evolving development strategy, further columns are output, as well. All
//' output matrices are sparse matrices.
//' 
//' @section Notes:
//' The data frame AllStages introduces variables used in size and fecundity
//' calculations. This DataFrame is broken up into long vectors composed of
//' input sizes and related variables for these calculations. The "model" Lists
//' bring in the vital rate models, and include random coefficients where
//' needed. We also have a number of extra variables, that include such info as
//' whether to use the Poisson, negative binomial, Gamma, or Gaussian
//' distributions for size and fecundity calculations. If \code{sizedist},
//' \code{sizebdist}, \code{sizecdist}, or \code{fecdist} equals 0, 1, 2, or 3,
//' then the Poisson, negative binomial, Gaussian, or Gamma is used,
//' respectively.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.jerzeibalowski_sp)]]
List jerzeibalowski_sp(DataFrame AllStages, DataFrame stageframe, int matrixformat,
  List survproxy, List obsproxy, List sizeproxy, List sizebproxy,
  List sizecproxy, List repstproxy, List fecproxy, List jsurvproxy,
  List jobsproxy, List jsizeproxy, List jsizebproxy, List jsizecproxy,
  List jrepstproxy, List jmatstproxy, NumericVector f2_inda,
  NumericVector f1_inda, NumericVector f2_indb, NumericVector f1_indb,
  NumericVector f2_indc, NumericVector f1_indc, StringVector r2_inda,
  StringVector r1_inda, StringVector r2_indb, StringVector r1_indb,
  StringVector r2_indc, StringVector r1_indc, NumericVector dev_terms,
  double dens, double fecmod, double maxsize, double maxsizeb, double maxsizec,
  unsigned int firstage, unsigned int finalage, bool negfec, int yearnumber,
  int patchnumber, double exp_tol = 700.0, double theta_tol = 100000000.0,
  String ipm_method = "cdf", bool err_check = false, bool simplicity = false) {
  
  NumericMatrix out; // Initialization
  bool ipm_cdf = true;
  if (ipm_method == "midpoint") ipm_cdf = false;
  
  // Determines the size of the matrix
  StringVector stagenames = stageframe["stage"];
  int nostages = stagenames.length();
  unsigned long matrixdim {0};
  
  int nostages_counter = nostages;
  for (int i = 0; i < nostages_counter; i++) {
    if (stringcompare_hard(as<std::string>(stagenames(i)), "AlmostBorn")) nostages -= 1;  
    if (stringcompare_hard(as<std::string>(stagenames(i)), "Dead")) nostages -= 1;
  }
  
  if (matrixformat == 1) { // Ehrlen-format hMPM
    matrixdim = nostages * nostages;
  } else if (matrixformat == 2) { // deVries-format hMPM
    matrixdim = nostages * (nostages + 1);
  } else if (matrixformat == 3) { // ahMPM
    matrixdim = nostages;
  } else if (matrixformat == 4) { // age-by-stage MPM
    matrixdim = nostages * (finalage - firstage + 1);
  }
  
  // Proxy model imports and settings
  bool sizezero = as<bool>(sizeproxy["zero_inflated"]);
  bool sizebzero = as<bool>(sizebproxy["zero_inflated"]);
  bool sizeczero = as<bool>(sizecproxy["zero_inflated"]);
  bool feczero = as<bool>(fecproxy["zero_inflated"]);
  bool jsizezero = as<bool>(jsizeproxy["zero_inflated"]);
  bool jsizebzero = as<bool>(jsizebproxy["zero_inflated"]);
  bool jsizeczero = as<bool>(jsizecproxy["zero_inflated"]);
  
  bool sizetrunc = as<bool>(sizeproxy["zero_truncated"]);
  bool sizebtrunc = as<bool>(sizebproxy["zero_truncated"]);
  bool sizectrunc = as<bool>(sizecproxy["zero_truncated"]);
  bool fectrunc = as<bool>(fecproxy["zero_truncated"]);
  bool jsizetrunc = as<bool>(jsizeproxy["zero_truncated"]);
  bool jsizebtrunc = as<bool>(jsizebproxy["zero_truncated"]);
  bool jsizectrunc = as<bool>(jsizecproxy["zero_truncated"]);
  
  NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
  NumericVector obscoefs = as<NumericVector>(obsproxy["coefficients"]);
  NumericVector sizecoefs = as<NumericVector>(sizeproxy["coefficients"]);
  NumericVector sizebcoefs = as<NumericVector>(sizebproxy["coefficients"]);
  NumericVector sizeccoefs = as<NumericVector>(sizecproxy["coefficients"]);
  NumericVector repstcoefs = as<NumericVector>(repstproxy["coefficients"]);
  NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
  NumericVector jsurvcoefs = as<NumericVector>(jsurvproxy["coefficients"]);
  NumericVector jobscoefs = as<NumericVector>(jobsproxy["coefficients"]);
  NumericVector jsizecoefs = as<NumericVector>(jsizeproxy["coefficients"]);
  NumericVector jsizebcoefs = as<NumericVector>(jsizebproxy["coefficients"]);
  NumericVector jsizeccoefs = as<NumericVector>(jsizecproxy["coefficients"]);
  NumericVector jrepstcoefs = as<NumericVector>(jrepstproxy["coefficients"]);
  NumericVector jmatstcoefs = as<NumericVector>(jmatstproxy["coefficients"]);
  
  double survsigma = as<double>(survproxy["sigma"]);
  double obssigma = as<double>(obsproxy["sigma"]);
  double sizesigma = as<double>(sizeproxy["sigma"]);
  double sizebsigma = as<double>(sizebproxy["sigma"]);
  double sizecsigma = as<double>(sizecproxy["sigma"]);
  double repstsigma = as<double>(repstproxy["sigma"]);
  double fecsigma = as<double>(fecproxy["sigma"]);
  double jsurvsigma = as<double>(jsurvproxy["sigma"]);
  double jobssigma = as<double>(jobsproxy["sigma"]);
  double jsizesigma = as<double>(jsizeproxy["sigma"]);
  double jsizebsigma = as<double>(jsizebproxy["sigma"]);
  double jsizecsigma = as<double>(jsizecproxy["sigma"]);
  double jrepstsigma = as<double>(jrepstproxy["sigma"]);
  double jmatstsigma = as<double>(jmatstproxy["sigma"]);
  
  int survdist = as<int>(survproxy["dist"]);
  int obsdist = as<int>(obsproxy["dist"]);
  int sizedist = as<int>(sizeproxy["dist"]);
  int sizebdist = as<int>(sizebproxy["dist"]);
  int sizecdist = as<int>(sizecproxy["dist"]);
  int repstdist = as<int>(repstproxy["dist"]);
  int fecdist = as<int>(fecproxy["dist"]);
  int jsurvdist = as<int>(jsurvproxy["dist"]);
  int jobsdist = as<int>(jobsproxy["dist"]);
  int jsizedist = as<int>(jsizeproxy["dist"]);
  int jsizebdist = as<int>(jsizebproxy["dist"]);
  int jsizecdist = as<int>(jsizecproxy["dist"]);
  int jrepstdist = as<int>(jrepstproxy["dist"]);
  int jmatstdist = as<int>(jmatstproxy["dist"]);
  
  if (NumericVector::is_na(sizesigma)) {
    if (sizedist == 1) {
      sizesigma = 1.0;
    } else {
      sizesigma = 0.0;
    }
  }
  if (NumericVector::is_na(sizebsigma)) {
    if (sizebdist == 1) {
      sizebsigma = 1.0;
    } else {
      sizebsigma = 0.0;
    }
  }
  if (NumericVector::is_na(sizecsigma)) {
    if (sizecdist == 1) {
      sizecsigma = 1.0;
    } else {
      sizecsigma = 0.0;
    }
  }
  if (NumericVector::is_na(fecsigma)) {
    if (fecdist == 1) {
      fecsigma = 1.0;
    } else {
      fecsigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizesigma)) {
    if (sizedist == 1) {
      jsizesigma = 1.0;
    } else {
      jsizesigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizebsigma)) {
    if (sizebdist == 1) {
      jsizebsigma = 1.0;
    } else {
      jsizebsigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizecsigma)) {
    if (sizecdist == 1) {
      jsizecsigma = 1.0;
    } else {
      jsizecsigma = 0.0;
    }
  }
  
  NumericMatrix vital_year = revelations(survproxy, obsproxy, sizeproxy, sizebproxy,
    sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy, jsizeproxy,
    jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 1);
  
  NumericMatrix vital_patch = revelations(survproxy, obsproxy, sizeproxy,
    sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
    jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 2);
  
  arma::imat rand_index = foi_index(survproxy, obsproxy, sizeproxy, sizebproxy,
    sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy, jsizeproxy,
    jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy);
  
  // NumericVector imports from model_proxy objects
  NumericVector sizeyearzi = as<NumericVector>(sizeproxy["zeroyear"]);
  NumericVector sizebyearzi = as<NumericVector>(sizebproxy["zeroyear"]);
  NumericVector sizecyearzi = as<NumericVector>(sizecproxy["zeroyear"]);
  NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
  NumericVector jsizeyearzi = as<NumericVector>(jsizeproxy["zeroyear"]);
  NumericVector jsizebyearzi = as<NumericVector>(jsizebproxy["zeroyear"]);
  NumericVector jsizecyearzi = as<NumericVector>(jsizecproxy["zeroyear"]);
  
  NumericVector dud_yearzi(sizeyearzi.length());
  
  NumericVector unisyzi = unique(sizeyearzi);
  NumericVector unisyzbi = unique(sizebyearzi);
  NumericVector unisyzci = unique(sizecyearzi);
  NumericVector unijsyzi = unique(jsizeyearzi);
  NumericVector unijsyzbi = unique(jsizebyearzi);
  NumericVector unijsyzci = unique(jsizecyearzi);
  NumericVector unifeci = unique(fecyearzi);
  
  NumericVector sizepatchzi = as<NumericVector>(sizeproxy["zeropatch"]);
  NumericVector sizebpatchzi = as<NumericVector>(sizebproxy["zeropatch"]);
  NumericVector sizecpatchzi = as<NumericVector>(sizecproxy["zeropatch"]);
  NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
  NumericVector jsizepatchzi = as<NumericVector>(jsizeproxy["zeropatch"]);
  NumericVector jsizebpatchzi = as<NumericVector>(jsizebproxy["zeropatch"]);
  NumericVector jsizecpatchzi = as<NumericVector>(jsizecproxy["zeropatch"]);
  
  NumericVector dud_patchzi(sizepatchzi.length());
  
  NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
  NumericVector obsgroups2 = as<NumericVector>(obsproxy["groups2"]);
  NumericVector sizegroups2 = as<NumericVector>(sizeproxy["groups2"]);
  NumericVector sizebgroups2 = as<NumericVector>(sizebproxy["groups2"]);
  NumericVector sizecgroups2 = as<NumericVector>(sizecproxy["groups2"]);
  NumericVector repstgroups2 = as<NumericVector>(repstproxy["groups2"]);
  NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
  NumericVector jsurvgroups2 = as<NumericVector>(jsurvproxy["groups2"]);
  NumericVector jobsgroups2 = as<NumericVector>(jobsproxy["groups2"]);
  NumericVector jsizegroups2 = as<NumericVector>(jsizeproxy["groups2"]);
  NumericVector jsizebgroups2 = as<NumericVector>(jsizebproxy["groups2"]);
  NumericVector jsizecgroups2 = as<NumericVector>(jsizecproxy["groups2"]);
  NumericVector jrepstgroups2 = as<NumericVector>(jrepstproxy["groups2"]);
  NumericVector jmatstgroups2 = as<NumericVector>(jmatstproxy["groups2"]);
  
  NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
  NumericVector obsgroups1 = as<NumericVector>(obsproxy["groups1"]);
  NumericVector sizegroups1 = as<NumericVector>(sizeproxy["groups1"]);
  NumericVector sizebgroups1 = as<NumericVector>(sizebproxy["groups1"]);
  NumericVector sizecgroups1 = as<NumericVector>(sizecproxy["groups1"]);
  NumericVector repstgroups1 = as<NumericVector>(repstproxy["groups1"]);
  NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
  NumericVector jsurvgroups1 = as<NumericVector>(jsurvproxy["groups1"]);
  NumericVector jobsgroups1 = as<NumericVector>(jobsproxy["groups1"]);
  NumericVector jsizegroups1 = as<NumericVector>(jsizeproxy["groups1"]);
  NumericVector jsizebgroups1 = as<NumericVector>(jsizebproxy["groups1"]);
  NumericVector jsizecgroups1 = as<NumericVector>(jsizecproxy["groups1"]);
  NumericVector jrepstgroups1 = as<NumericVector>(jrepstproxy["groups1"]);
  NumericVector jmatstgroups1 = as<NumericVector>(jmatstproxy["groups1"]);
  
  NumericVector sizegroups2zi = as<NumericVector>(sizeproxy["zerogroups2"]);
  NumericVector sizebgroups2zi = as<NumericVector>(sizebproxy["zerogroups2"]);
  NumericVector sizecgroups2zi = as<NumericVector>(sizecproxy["zerogroups2"]);
  NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
  NumericVector jsizegroups2zi = as<NumericVector>(jsizeproxy["zerogroups2"]);
  NumericVector jsizebgroups2zi = as<NumericVector>(jsizebproxy["zerogroups2"]);
  NumericVector jsizecgroups2zi = as<NumericVector>(jsizecproxy["zerogroups2"]);
  
  NumericVector dud_groups2zi(jsizecyearzi.length());
  
  NumericVector sizegroups1zi = as<NumericVector>(sizeproxy["zerogroups1"]);
  NumericVector sizebgroups1zi = as<NumericVector>(sizebproxy["zerogroups1"]);
  NumericVector sizecgroups1zi = as<NumericVector>(sizecproxy["zerogroups1"]);
  NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
  NumericVector jsizegroups1zi = as<NumericVector>(jsizeproxy["zerogroups1"]);
  NumericVector jsizebgroups1zi = as<NumericVector>(jsizebproxy["zerogroups1"]);
  NumericVector jsizecgroups1zi = as<NumericVector>(jsizecproxy["zerogroups1"]);
  
  NumericVector dud_groups1zi(jsizecyearzi.length());
  
  NumericVector survind = flightoficarus(survproxy);
  NumericVector obsind = flightoficarus(obsproxy);
  NumericVector sizeind = flightoficarus(sizeproxy);
  NumericVector sizebind = flightoficarus(sizebproxy);
  NumericVector sizecind = flightoficarus(sizecproxy);
  NumericVector repstind = flightoficarus(repstproxy);
  NumericVector fecind = flightoficarus(fecproxy);
  NumericVector jsurvind = flightoficarus(jsurvproxy);
  NumericVector jobsind = flightoficarus(jobsproxy);
  NumericVector jsizeind = flightoficarus(jsizeproxy);
  NumericVector jsizebind = flightoficarus(jsizebproxy);
  NumericVector jsizecind = flightoficarus(jsizecproxy);
  NumericVector jrepstind = flightoficarus(jrepstproxy);
  NumericVector jmatstind = flightoficarus(jmatstproxy);
  
  NumericVector sizeindzi = zero_flightoficarus(sizeproxy);
  NumericVector sizebindzi = zero_flightoficarus(sizebproxy);
  NumericVector sizecindzi = zero_flightoficarus(sizecproxy);
  NumericVector fecindzi = zero_flightoficarus(fecproxy);
  NumericVector jsizeindzi = zero_flightoficarus(jsizeproxy);
  NumericVector jsizebindzi = zero_flightoficarus(jsizebproxy);
  NumericVector jsizecindzi = zero_flightoficarus(jsizecproxy);
  
  StringVector survind_rownames = bootson(survproxy);
  StringVector obsind_rownames = bootson(obsproxy);
  StringVector sizeind_rownames = bootson(sizeproxy);
  StringVector sizebind_rownames = bootson(sizebproxy);
  StringVector sizecind_rownames = bootson(sizecproxy);
  StringVector repstind_rownames = bootson(repstproxy);
  StringVector fecind_rownames = bootson(fecproxy);
  StringVector jsurvind_rownames = bootson(jsurvproxy);
  StringVector jobsind_rownames = bootson(jobsproxy);
  StringVector jsizeind_rownames = bootson(jsizeproxy);
  StringVector jsizebind_rownames = bootson(jsizebproxy);
  StringVector jsizecind_rownames = bootson(jsizecproxy);
  StringVector jrepstind_rownames = bootson(jrepstproxy);
  StringVector jmatstind_rownames = bootson(jmatstproxy);
  
  StringVector sizeind_rownames_zi = zero_bootson(sizeproxy);
  StringVector sizebind_rownames_zi = zero_bootson(sizebproxy);
  StringVector sizecind_rownames_zi = zero_bootson(sizecproxy);
  StringVector fecind_rownames_zi = zero_bootson(fecproxy);
  StringVector jsizeind_rownames_zi = zero_bootson(jsizeproxy);
  StringVector jsizebind_rownames_zi = zero_bootson(jsizebproxy);
  StringVector jsizecind_rownames_zi = zero_bootson(jsizecproxy);
  
  // AllStages import and settings
  Rcpp::NumericVector stage3_num = as<NumericVector>(AllStages["stage3"]);
  Rcpp::NumericVector stage2n_num = as<NumericVector>(AllStages["stage2n"]);
  Rcpp::NumericVector stage2o_num = as<NumericVector>(AllStages["stage2o"]);
  arma::vec stage3 = as<arma::vec>(stage3_num);
  arma::vec stage2n = as<arma::vec>(stage2n_num);
  arma::vec stage2o = as<arma::vec>(stage2o_num);
  
  Rcpp::NumericVector sz3 = as<NumericVector>(AllStages["size3"]);
  Rcpp::NumericVector sz2n = as<NumericVector>(AllStages["size2n"]);
  Rcpp::NumericVector sz2o = as<NumericVector>(AllStages["size2o"]);
  Rcpp::NumericVector sz1 = as<NumericVector>(AllStages["size1"]);
  Rcpp::NumericVector szb3 = as<NumericVector>(AllStages["sizeb3"]);
  Rcpp::NumericVector szb2n = as<NumericVector>(AllStages["sizeb2n"]);
  Rcpp::NumericVector szb2o = as<NumericVector>(AllStages["sizeb2o"]);
  Rcpp::NumericVector szb1 = as<NumericVector>(AllStages["sizeb1"]);
  Rcpp::NumericVector szc3 = as<NumericVector>(AllStages["sizec3"]);
  Rcpp::NumericVector szc2n = as<NumericVector>(AllStages["sizec2n"]);
  Rcpp::NumericVector szc2o = as<NumericVector>(AllStages["sizec2o"]);
  Rcpp::NumericVector szc1 = as<NumericVector>(AllStages["sizec1"]);
  Rcpp::NumericVector ob3 = as<NumericVector>(AllStages["obs3"]);
  Rcpp::NumericVector fl3 = as<NumericVector>(AllStages["rep3"]);
  Rcpp::NumericVector fl2n = as<NumericVector>(AllStages["rep2n"]);
  Rcpp::NumericVector fl2o = as<NumericVector>(AllStages["rep2o"]);
  Rcpp::NumericVector fl1 = as<NumericVector>(AllStages["rep1"]);
  Rcpp::NumericVector mat3 = as<NumericVector>(AllStages["mat3"]);
  Rcpp::NumericVector mat2n = as<NumericVector>(AllStages["mat2n"]);
  Rcpp::NumericVector mat2o = as<NumericVector>(AllStages["mat2o"]);
  Rcpp::NumericVector mat1 = as<NumericVector>(AllStages["mat1"]);
  Rcpp::NumericVector immat2n = as<NumericVector>(AllStages["imm2n"]);
  Rcpp::NumericVector immat2o = as<NumericVector>(AllStages["imm2o"]);
  Rcpp::NumericVector immat1 = as<NumericVector>(AllStages["imm1"]);
  
  Rcpp::NumericVector repentry = as<NumericVector>(AllStages["repentry3"]);
  Rcpp::NumericVector indata2n = as<NumericVector>(AllStages["indata2n"]);
  Rcpp::NumericVector indata2o = as<NumericVector>(AllStages["indata2o"]);
  Rcpp::NumericVector binwidth3 = as<NumericVector>(AllStages["binwidth"]);
  Rcpp::NumericVector binbwidth3 = as<NumericVector>(AllStages["binbwidth"]);
  Rcpp::NumericVector bincwidth3 = as<NumericVector>(AllStages["bincwidth"]);
  Rcpp::NumericVector actualage2 = as<NumericVector>(AllStages["actualage"]);
  
  Rcpp::NumericVector grp3 = as<NumericVector>(AllStages["group3"]);
  Rcpp::NumericVector grp2n = as<NumericVector>(AllStages["group2n"]);
  Rcpp::NumericVector grp2o = as<NumericVector>(AllStages["group2o"]);
  Rcpp::NumericVector grp1 = as<NumericVector>(AllStages["group1"]);
  
  Rcpp::NumericVector ovestt_num = as<NumericVector>(AllStages["ovest_t"]);
  arma::vec ovestt = as<arma::vec>(ovestt_num);
  
  Rcpp::NumericVector ovestf_num = as<NumericVector>(AllStages["ovest_f"]);
  arma::vec ovestf = as<arma::vec>(ovestf_num);
  
  Rcpp::NumericVector indata = as<NumericVector>(AllStages["indata"]);
  Rcpp::NumericVector ovgivent = as<NumericVector>(AllStages["ovgiven_t"]);
  Rcpp::NumericVector ovgivenf = as<NumericVector>(AllStages["ovgiven_f"]);
  
  Rcpp::NumericVector ovsurvmult = as<NumericVector>(AllStages["ovsurvmult"]);
  Rcpp::NumericVector ovfecmult = as<NumericVector>(AllStages["ovfecmult"]);
  
  Rcpp::IntegerVector index321_int = as<IntegerVector>(AllStages["index321"]);
  arma::uvec index321 = as<arma::uvec>(index321_int);
  
  Rcpp::IntegerVector aliveandequal = as<IntegerVector>(AllStages["aliveandequal"]); // Used to be NumericVector
  
  int n = stage3.n_elem;
  
  arma::uvec replacetvec = find(ovestt != -1.0);
  arma::uvec replacefvec = find(ovestf != -1.0);
  int replacementst = replacetvec.n_elem;
  int replacementsf = replacefvec.n_elem;
  int repindex {0};
  int properindex {0};
  int proxyindex {0};
  
  // Determination of choices of fixed and random individual covariates
  double inda1 = f1_inda(yearnumber);
  double indb1 = f1_indb(yearnumber);
  double indc1 = f1_indc(yearnumber);
  double inda2 = f2_inda(yearnumber);
  double indb2 = f2_indb(yearnumber);
  double indc2 = f2_indc(yearnumber);
  
  String chosen_r2inda = r2_inda(yearnumber);
  String chosen_r1inda = r1_inda(yearnumber);
  String chosen_r2indb = r2_indb(yearnumber);
  String chosen_r1indb = r1_indb(yearnumber);
  String chosen_r2indc = r2_indc(yearnumber);
  String chosen_r1indc = r1_indc(yearnumber);
  
  // The output matrix to collect conditional probabilities
  // Matrix out is 0 matrix with n rows & 6 columns: 0 surv, 1 obs, 2 repst,
  // 3 size, 4 size_b, 5 size_c, 6 matst, >6 are test variables
  if (err_check) {
    NumericMatrix zeroform(n, 7);
    out = zeroform;
    CharacterVector out_names = {"surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst"};
    colnames(out) = out_names;
  }
  NumericVector out_vec(7);
  
  arma::sp_mat survtransmat(matrixdim, matrixdim);
  arma::sp_mat fectransmat(matrixdim, matrixdim);
  
  double fec_addedcoefs = sum(feccoefs);
  double jsurv_coefsadded = sum(jsurvcoefs);
  double mat_predicted {0.0};
  unsigned int k {0};
  // The following loop runs through each line of AllStages, and so runs through
  // each estimable element in the matrix
  for(int i = 0; i < n; i++) {
    out_vec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    k = aliveandequal(i);
    mat_predicted = 0.0;
    
    if (err_check) out(i, 6) = 1.0; // Initialization of maturity status probability for typical case
    
    Rcpp::NumericVector statusterms = {fl1(i), fl2n(i), sz1(i), sz2o(i),
      szb1(i), szb2o(i), szc1(i), szc2o(i), actualage2(i), inda1, inda2, indb1,
      indb2, indc1, indc2, dens, sz3(i), szb3(i), szc3(i), binwidth3(i),
      binbwidth3(i), bincwidth3(i)};
    
    if (ovgivent(i) == -1 && indata(i) == 1 && stage2n(i) == stage2o(i)) {
      if ((mat2n(i) == 1 && mat3(i) == 1) || (mat2o(i) == 1 && mat3(i) == 1)) {
        
        // Adult survival transitions
        if (survdist < 5) {
          out_vec(0) = preouterator(survproxy, survcoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, survgroups2,
            survgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            survind, survind_rownames, sizeindzi, sizeind_rownames_zi, false,
            survsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 1, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
        } else {
          out_vec(0) = survcoefs(0);
        }
        if (err_check) out(i, 0) = out_vec(0);
        
        if (obsdist < 5) {
          out_vec(1) = preouterator(obsproxy, obscoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, obsgroups2,
            obsgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            obsind, obsind_rownames, sizeindzi, sizeind_rownames_zi, false,
            obssigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 2, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
          
        } else {
          out_vec(1) = obscoefs(0);
        }
        if (err_check) out(i, 1) = out_vec(1);
        
        if (ob3(i) == 1 || obsdist == 5) {
          
          if (sizedist < 5) {
            bool used_sizezero = false;
            if (sizezero && sz3(i) == 0) used_sizezero = sizezero;
            
            out_vec(3) = preouterator(sizeproxy, sizecoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, sizegroups2, sizegroups1, sizegroups2zi, sizegroups1zi,
              sizeyearzi, sizepatchzi, sizeind, sizeind_rownames, sizeindzi,
              sizeind_rownames_zi, used_sizezero, sizesigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, sizedist, 3, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages,
              sizetrunc);
              
          } else {
            out_vec(3) = 1.0;
          }
          if (err_check) out(i, 3) = out_vec(3);
          
          if (sizebdist < 5) {
            bool used_sizebzero = false;
            if (sizebzero && szb3(i) == 0) used_sizebzero = sizebzero;
            
            out_vec(4) = preouterator(sizebproxy, sizebcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, sizebgroups2, sizebgroups1, sizebgroups2zi,
              sizebgroups1zi, sizebyearzi, sizebpatchzi, sizebind,
              sizebind_rownames, sizebindzi, sizebind_rownames_zi, used_sizebzero,
              sizebsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist,
              4, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, sizebtrunc);
          } else {
            out_vec(4) = 1.0;
          }
          if (err_check) out(i, 4) = out_vec(4);
          
          if (sizecdist < 5) {
            bool used_sizeczero = false;
            if (sizeczero && szc3(i) == 0) used_sizeczero = sizeczero;
            
            out_vec(5) = preouterator(sizecproxy, sizeccoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, sizecgroups2, sizecgroups1, sizecgroups2zi,
              sizecgroups1zi, sizecyearzi, sizecpatchzi, sizecind,
              sizecind_rownames, sizecindzi, sizecind_rownames_zi, used_sizeczero,
              sizecsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist,
              5, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, sizectrunc);
          } else {
            out_vec(5) = 1.0;
          }
          if (err_check) out(i, 5) = out_vec(5);
          
          if (repstdist < 5) {
            out_vec(2) = preouterator(repstproxy, repstcoefs, rand_index,
              dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
              chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
              statusterms, repstgroups2, repstgroups1, dud_groups2zi, dud_groups1zi,
              dud_yearzi, dud_patchzi, repstind, repstind_rownames, sizeindzi,
              sizeind_rownames_zi, false, repstsigma, grp2o(i), grp1(i),
              patchnumber, yearnumber, 4, 6, exp_tol, theta_tol, ipm_cdf,
              matrixformat, fecmod, repentry(i), negfec, stage2n(i), nostages, 0);
              
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - out_vec(2);
            }
          } else {
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - repstcoefs(0);
            } else if (fl3(i) == 1) {
              out_vec(2) = repstcoefs(0);
            } else {
              out_vec(2) = 0.0;
            }
          }
          if (err_check) out(i, 2) = out_vec(2);
          
        } else {
          out_vec(1) = 1.0 - out_vec(1);
          out_vec(2) = 1.0;
          out_vec(3) = 1.0;
          out_vec(4) = 1.0;
          out_vec(5) = 1.0;
          out_vec(6) = 1.0;
          
          if (err_check) {
            out(i, 1) = out_vec(1);
            out(i, 2) = out_vec(2);
            out(i, 3) = out_vec(3);
            out(i, 4) = out_vec(4);
            out(i, 5) = out_vec(5);
            out(i, 6) = out_vec(6);
          }
        }
        survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
          out_vec(4) * out_vec(5) * out_vec(6);
        
      } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurv_coefsadded != 0.0) {
        // Juvenile to adult transitions
        
        if (jmatstdist < 5) {
          mat_predicted = preouterator(jmatstproxy, jmatstcoefs, rand_index,
            dev_terms, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
            chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms,
            jmatstgroups2, jmatstgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi,
            dud_patchzi, jmatstind, jmatstind_rownames, jsizeindzi,
            jsizeind_rownames_zi, false, jmatstsigma, grp2o(i), grp1(i),
            patchnumber, yearnumber, 4, 21, exp_tol, theta_tol, ipm_cdf, matrixformat,
            fecmod, repentry(i), negfec, stage2n(i), nostages, 0);
          
          if (mat3(i) > 0.5) {
            out_vec(6) = mat_predicted;
          } else {
            out_vec(6) = 1 - mat_predicted;
          }
        } else {
          if (mat3(i) > 0.5) {
            out_vec(6) = 1;
          } else {
            out_vec(6) = 0;
          }
        }
        if (err_check) out(i, 6) = out_vec(6);
        
        if (jsurvdist < 5) {
          out_vec(0) = preouterator(jsurvproxy, jsurvcoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsurvgroups2,
            jsurvgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            jsurvind, jsurvind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
            jsurvsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 8, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
        } else {
          out_vec(0) = jsurvcoefs(0);
        }
        if (err_check) out(i, 0) = out_vec(0);
        
        if (jobsdist < 5) {
          out_vec(1) = preouterator(jobsproxy, jobscoefs, rand_index, dev_terms,
            vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
            chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jobsgroups2,
            jobsgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
            jobsind, jobsind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
            jobssigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 9, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0);
        } else {
          out_vec(1) = jobscoefs(0);
        }
        if (err_check) out(i, 1) = out_vec(1);
        
        if (ob3(i) == 1 || jobsdist == 5) {
          if (jsizedist < 5) {
            out_vec(3) = preouterator(jsizeproxy, jsizecoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsizegroups2,
              jsizegroups1, jsizegroups2zi, jsizegroups1zi, jsizeyearzi, jsizepatchzi,
              jsizeind, jsizeind_rownames, jsizeindzi, jsizeind_rownames_zi, jsizezero,
              jsizesigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizedist, 10,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, jsizetrunc);
          } else {
            out_vec(3) = 1.0;
          }
          if (err_check) out(i, 3) = out_vec(3);
          
          if (jsizebdist < 5) {
            out_vec(4) = preouterator(jsizebproxy, jsizebcoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsizebgroups2,
              jsizebgroups1, jsizebgroups2zi, jsizebgroups1zi, jsizebyearzi, jsizebpatchzi,
              jsizebind, jsizebind_rownames, jsizebindzi, jsizebind_rownames_zi, jsizebzero,
              jsizebsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist, 11,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, jsizebtrunc);
          } else {
            out_vec(4) = 1.0;
          }
          if (err_check) out(i, 4) = out_vec(4);
          
          if (jsizecdist < 5) {
            out_vec(5) = preouterator(jsizecproxy, jsizeccoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda,chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jsizecgroups2,
              jsizecgroups1, jsizecgroups2zi, jsizecgroups1zi, jsizecyearzi, jsizecpatchzi,
              jsizecind, jsizecind_rownames, jsizecindzi, jsizecind_rownames_zi, jsizeczero,
              jsizecsigma, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist, 12,
              exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, jsizectrunc);
          } else {
            out_vec(5) = 1.0;
          }
          if (err_check) out(i, 5) = out_vec(5);
          
          if (jrepstdist < 5) {
            out_vec(2) = preouterator(jrepstproxy, jrepstcoefs, rand_index, dev_terms,
              vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
              chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, jrepstgroups2,
              jrepstgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi,
              jrepstind, jrepstind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              jrepstsigma, grp2o(i), grp1(i), patchnumber, yearnumber, 4, 13, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0);
              
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - out_vec(2);
            }
          } else {
            if (fl3(i) == 0) {
              out_vec(2) = 1.0 - jrepstcoefs(0);
            } else if (fl3(i) == 1) {
              out_vec(2) = jrepstcoefs(0);
            } else {
              out_vec(2) = 0.0;
            }
          }
          if (err_check) out(i, 2) = out_vec(2);
          
        } else {
          out_vec(1) = 1.0 - out_vec(1);
          out_vec(2) = 1.0;
          out_vec(3) = 1.0;
          out_vec(4) = 1.0;
          out_vec(5) = 1.0;
          out_vec(6) = 1.0;
          
          if (err_check) {
            out(i, 1) = out_vec(1);
            out(i, 2) = out_vec(2);
            out(i, 3) = out_vec(3);
            out(i, 4) = out_vec(4);
            out(i, 5) = out_vec(5);
            out(i, 6) = out_vec(6);
          }
        }
        
        survtransmat(k) = out_vec(0) * out_vec(1) * out_vec(2) * out_vec(3) *
          out_vec(4) * out_vec(5) * out_vec(6);
      }
    } else if (ovgivent(i) != -1) {
      // All other transitions
      
      survtransmat(k) = ovgivent(i);
    }
    
    // This next block calculates fecundity
    if (indata2n(i) == 1 && fec_addedcoefs != 0.0) {
      if (fl2o(i) > 0.0 && ovgivenf(i) == -1.0) {
        
        fectransmat(k) = preouterator(fecproxy, feccoefs, rand_index, dev_terms,
          vital_year, vital_patch, chosen_r2inda, chosen_r1inda, chosen_r2indb,
          chosen_r1indb, chosen_r2indc, chosen_r1indc, statusterms, fecgroups2,
          fecgroups1, fecgroups2zi, fecgroups1zi, fecyearzi, fecpatchzi, fecind,
          fecind_rownames, fecindzi, fecind_rownames_zi, feczero, fecsigma,
          grp2o(i), grp1(i), patchnumber, yearnumber, fecdist, 7, exp_tol,
          theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
          stage2n(i), nostages, fectrunc);
        
      } else if (ovgivenf(i) != -1 ) {
        fectransmat(k) = ovgivenf(i);
      }
    } else if (ovgivenf(i) != -1 ) {
      fectransmat(k) = ovgivenf(i);
    }
  }
  
  double ov_mult {0};
  if (replacementst > 0) {
    for (int i = 0; i < replacementst; i++) {
      
      repindex = replacetvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestt(repindex));
      
      if (rightindex.n_elem > 0) {
        proxyindex = aliveandequal(rightindex(0));
        
        ov_mult = ovsurvmult(i);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        survtransmat(properindex) = survtransmat(proxyindex) * ov_mult;
      }
    }
  }
  
  if (replacementsf > 0) {
    for (int i = 0; i < replacementsf; i++) {
      
      repindex = replacefvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestf(repindex));
      
      if (rightindex.n_elem > 0) {
        proxyindex = aliveandequal(rightindex(0));
        
        ov_mult = ovfecmult(i);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        fectransmat(properindex) = fectransmat(proxyindex) * ov_mult;
      }
    }
  }
  
  // Final output
  List output(4);
  
  if (!simplicity) {
    arma::sp_mat amatrix = survtransmat + fectransmat;
    output(0) = amatrix;
  } else {
    output(0) = R_NilValue;
  }
  
  output(1) = survtransmat;
  output(2) = fectransmat;
  
  if (err_check) {
    output(3) = out;
  } else {
    output(3) = R_NilValue;
  }
  CharacterVector output_names = {"A", "U", "F", "out"};
  output.attr("names") = output_names;
  
  return output;
}

//' Create Historically Structured Version of ahMPM
//' 
//' Function \code{.thefifthhousemate()} takes an ahistorical MPM as input, and
//' uses the \code{allstages} index to create a historically structured version
//' of it.
//' 
//' @name thefifthhousemate
//' 
//' @param mpm The original ahMPM, supplied as a \code{lefkoMat} object.
//' @param allstages The index dataframe developed by
//' \code{\link{.simplepizzle}()}.
//' @param stageframe The ahistorical stageframe supplied by
//' \code{\link{.simplepizzle}()}.
//' @param format Integer indicating whether historical matrices should be in
//' (1) Ehrlen or (2) deVries format.
//' 
//' @return This will return a list of lists. The first list is composed of all
//' new \code{A} matrices. The second list is composed of all new \code{U}
//' matrices. The third list is composed of all new \code{F} matrices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.thefifthhousemate)]]
Rcpp::List thefifthhousemate (List mpm, DataFrame allstages,
  DataFrame stageframe, int format) {
  Rcpp::List old_Umats = mpm["U"];
  Rcpp::List old_Fmats = mpm["F"];
  
  Rcpp::IntegerVector stageid = stageframe["stage_id"];
  int nostages = stageid.length();
  int nocols = nostages * nostages;
  
  if (format == 2) nocols = (nostages -1) * nostages;
  
  Rcpp::IntegerVector old_index = allstages["index21"];
  Rcpp::IntegerVector new_indexu = allstages["index321u"];
  Rcpp::IntegerVector new_indexf = allstages["index321f"];
  
  int num_mats = old_Umats.length();
  int index_elems = new_indexu.length();
  
  Rcpp::List new_Umats(num_mats);
  Rcpp::List new_Fmats(num_mats);
  Rcpp::List new_Amats(num_mats);
  
  arma::mat new_U(nocols, nocols, fill::zeros);
  arma::mat new_F(nocols, nocols, fill::zeros);
  arma::mat new_A(nocols, nocols, fill::zeros);
  arma::mat old_U(nostages, nostages, fill::zeros);
  arma::mat old_F(nostages, nostages, fill::zeros);
  
  for (int i = 0; i < num_mats; i++) {
    new_U.zeros();
    new_F.zeros();
    new_A.zeros();
    old_U.zeros();
    old_F.zeros();
    
    old_U = as<arma::mat>(old_Umats(i));
    old_F = as<arma::mat>(old_Fmats(i));
    
    for (int j = 0; j < index_elems; j++) {
      if (new_indexu(j) > -1.0) {
        new_U(new_indexu(j)) = old_U(old_index(j));
      }
      if (new_indexf(j) > -1.0) {
        new_F(new_indexf(j)) = old_F(old_index(j));
      }
    }
    
    new_A = new_U + new_F;
    new_Umats(i) = new_U;
    new_Fmats(i) = new_F;
    new_Amats(i) = new_A;
  }
  
  Rcpp::List output = List::create(Named("A") = new_Amats, _["U"] = new_Umats,
    _["F"] = new_Fmats);
  return output;
}

//' Creates Matrices of Year and Patch Terms in Leslie Models
//' 
//' Function \code{revelations_leslie()} creates a matrix holding either the
//' year or patch coefficients from Leslie vital rate models. This reduces
//' memory load in function \code{\link{motherbalowski}()}.
//' 
//' @name revelations_leslie
//' 
//' @param survproxy The proxy vital rate model covering survival from the main
//' matrix estimator function.
//' @param fecproxy The proxy vital rate model covering fecundity from the main
//' matrix estimator function.
//' 
//' @return A matrix with 2 columns corresponding to the number of vital rates
//' and number of columns equal to the number of year or patches.
//' 
//' @keywords internal
//' @noRd
NumericMatrix revelations_leslie(List survproxy, List fecproxy, int mat_switch) {
  
  NumericMatrix final_mat;
  
  if (mat_switch == 1) {
    NumericVector survyear = as<NumericVector>(survproxy["years"]);
    NumericVector fecyear = as<NumericVector>(fecproxy["years"]);
    
    int matrows = survyear.length();
    
    NumericMatrix year_mat(matrows, 2);
    year_mat(_, 0) = survyear;
    year_mat(_, 1) = fecyear;
    
    final_mat = year_mat;
    
  } else if (mat_switch == 2) {
    
    NumericVector survpatch = as<NumericVector>(survproxy["patches"]);
    NumericVector fecpatch = as<NumericVector>(fecproxy["patches"]);
    
    int matrows = survpatch.length();
    
    NumericMatrix patch_mat(matrows, 2);
    patch_mat(_, 0) = survpatch;
    patch_mat(_, 1) = fecpatch;

    final_mat = patch_mat;
  }
  
  return final_mat;
}

//' Create Index of Element Numbers for Random Individual Covariate Terms in
//' Leslie Models
//' 
//' Function \code{foi_index_leslie()} creates a matrix indexing the end points
//' of each random individual covariate in the utilized vectors. Used in
//' function \code{\link{motherbalowski}()}.
//' 
//' @name foi_index_leslie
//' 
//' @param surv_proxy Adult survival model proxy.
//' @param fec_proxy Adult fecundity model proxy.
//' 
//' @return An integer matrix with 6 rows and 3 columns. The columns contain the
//' number of elements in each random individual covariate term, with the row
//' order being: 1) cov a t2, 2) cov a t1, 3) cov b t2, 4) cov b t1,
//' 5) cov c t2, and 6) cov c t1.
//' 
//' @keywords internal
//' @noRd
arma::imat foi_index_leslie(List surv_proxy, List fec_proxy) {
  
  arma::ivec surv_fc = foi_counter(surv_proxy, false);
  arma::ivec fec_fc = foi_counter(fec_proxy, false);
  arma::ivec fec_fc_zi = foi_counter(fec_proxy, true);
  
  arma::imat final_mat(6, 3, fill::zeros);
  
  for (int i = 0; i < 6; i++) {
    final_mat(i, 0) = surv_fc(i);
    final_mat(i, 1) = fec_fc(i);
    final_mat(i, 2) = fec_fc_zi(i);
  }
  
  return final_mat;
}

//' Estimate All Elements of Function-based Population Projection Matrix
//' 
//' Function \code{motherbalowski()} swiftly calculates matrix elements in
//' function-based Leslie population projection matrices. Used in
//' \code{\link{fleslie}()}.
//' 
//' @param actualages An integer vector of all actual ages to be included in the
//' matrices, in order.
//' @param ageframe The modified stageframe used in matrix calculations.
//' @param survproxy List of coefficients estimated in model of survival.
//' @param fecproxy List of coefficients estimated in model of fecundity.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param surv_dev A numeric value indicating the deviation to the linear
//' model of survival input by the user.
//' @param fec_dev A numeric value indicating the deviation to the linear
//' model of fecundity input by the user.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param finalage The final age to be included in Leslie MPM estimation.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to \code{0}.
//' @param yearnumber An integer specifying which time at time \emph{t} to
//' develop matrices for. Must be in reference to the \code{listofyears} object
//' developed in the \code{R} matrix estimator function.
//' @param patchnumber An integer specifying which patch to develop matrices
//' for. Must be in reference to the \code{listofyears} object developed in the
//' \code{R} matrix estimator function.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' 
//' @return A list of 3 matrices, including the main MPM (A), the survival-
//' transition matrix (U), and a fecundity matrix (F).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.motherbalowski)]]
List motherbalowski(IntegerVector actualages, DataFrame ageframe, List survproxy,
  List fecproxy, NumericVector f2_inda, NumericVector f1_inda,
  NumericVector f2_indb, NumericVector f1_indb, NumericVector f2_indc,
  NumericVector f1_indc, StringVector r2_inda, StringVector r1_inda,
  StringVector r2_indb, StringVector r1_indb, StringVector r2_indc,
  StringVector r1_indc, double surv_dev, double fec_dev, double dens,
  double fecmod, unsigned int finalage, bool negfec,
  int yearnumber, int patchnumber, double exp_tol = 700.0,
  double theta_tol = 100000000.0, bool simplicity = false) {
  
  // Determines the size of the matrix
  StringVector sf_agenames = as<StringVector>(ageframe["stage"]);
  IntegerVector sf_minage = as<IntegerVector>(ageframe["min_age"]);
  IntegerVector sf_maxage = as<IntegerVector>(ageframe["max_age"]);
  IntegerVector sf_repstatus = as<IntegerVector>(ageframe["repstatus"]);
  int noages = actualages.length();
  
  bool cont = false;
  if (sf_maxage(noages - 1) == NA_INTEGER) {
    cont = true;
  }
  
  // Proxy model imports and settings
  NumericVector survcoefs = as<NumericVector>(survproxy["coefficients"]);
  NumericVector feccoefs = as<NumericVector>(fecproxy["coefficients"]);
  
  bool feczero = as<bool>(fecproxy["zero_inflated"]);
  int survdist = as<int>(survproxy["dist"]);
  int fecdist = as<int>(fecproxy["dist"]);
  double fecsigma = as<double>(fecproxy["sigma"]);
  
  if (NumericVector::is_na(fecsigma)) {
    if (fecdist == 1) {
      fecsigma = 1.0;
    } else {
      fecsigma = 0.0;
    }
  }

  NumericMatrix vital_year = revelations_leslie(survproxy, fecproxy, 1);
  NumericMatrix vital_patch = revelations_leslie(survproxy, fecproxy, 2);
  
  NumericVector fecyearzi = as<NumericVector>(fecproxy["zeroyear"]);
  NumericVector fecpatchzi = as<NumericVector>(fecproxy["zeropatch"]);
  NumericVector survgroups2 = as<NumericVector>(survproxy["groups2"]);
  NumericVector fecgroups2 = as<NumericVector>(fecproxy["groups2"]);
  NumericVector survgroups1 = as<NumericVector>(survproxy["groups1"]);
  NumericVector fecgroups1 = as<NumericVector>(fecproxy["groups1"]);
  NumericVector fecgroups2zi = as<NumericVector>(fecproxy["zerogroups2"]);
  NumericVector fecgroups1zi = as<NumericVector>(fecproxy["zerogroups1"]);
  
  NumericVector survind = flightoficarus(survproxy);
  NumericVector fecind = flightoficarus(fecproxy);
  NumericVector fecindzi = zero_flightoficarus(fecproxy);
  
  arma::imat rand_index = foi_index_leslie(survproxy, fecproxy);
  
  StringVector survind_rownames = bootson(survproxy);
  StringVector fecind_rownames = bootson(fecproxy);
  StringVector fecind_rownames_zi = zero_bootson(fecproxy);
  
  // Determination of choices of fixed and random individual covariates
  double inda1 = f1_inda(yearnumber);
  double indb1 = f1_indb(yearnumber);
  double indc1 = f1_indc(yearnumber);
  double inda2 = f2_inda(yearnumber);
  double indb2 = f2_indb(yearnumber);
  double indc2 = f2_indc(yearnumber);
  
  String chosen_r2inda = r2_inda(yearnumber);
  String chosen_r1inda = r1_inda(yearnumber);
  String chosen_r2indb = r2_indb(yearnumber);
  String chosen_r1indb = r1_indb(yearnumber);
  String chosen_r2indc = r2_indc(yearnumber);
  String chosen_r1indc = r1_indc(yearnumber);
  
  // The output matrices
  arma::mat survtransmat(noages, noages, fill::zeros);
  arma::mat fectransmat(noages, noages, fill::zeros);
  
  // The following loop runs through each age, and so runs through
  // each estimable element in the matrix
  double fec_addedcoefs = sum(feccoefs);
  for(int i = 0; i < noages; i++) {
    // Adult survival transitions
    
    double preout {0.0};
    
    if (survdist < 5) {
      
      double chosen_randcova2 {0.0};
      if (chosen_r2inda != "none") {
        for (int indcount = 0; indcount < rand_index(0, 0); indcount++) {
          if (chosen_r2inda == survind_rownames(indcount)) {
            chosen_randcova2 = survind(indcount);
          }
        }
      }
      double chosen_randcova1 {0.0};
      if (chosen_r1inda != "none") {
        int delectable_sum = rand_index(0, 0);
        for (int indcount = 0; indcount < rand_index(1, 0); indcount++) {
          if (chosen_r1inda == survind_rownames(indcount + delectable_sum)) {
            chosen_randcova1 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovb2 {0.0};
      if (chosen_r2indb != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0);
        for (int indcount = 0; indcount < rand_index(2, 0); indcount++) {
          if (chosen_r2indb == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovb2 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovb1 {0.0};
      if (chosen_r1indb != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0);
        for (int indcount = 0; indcount < rand_index(3, 0); indcount++) {
          if (chosen_r1indb == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovb1 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovc2 {0.0};
      if (chosen_r2indc != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
          rand_index(3, 0);
        for (int indcount = 0; indcount < rand_index(4, 0); indcount++) {
          if (chosen_r2indc == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovc2 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovc1 {0.0};
      if (chosen_r1indc != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
          rand_index(3, 0) + rand_index(4, 0);
        for (int indcount = 0; indcount < rand_index(5, 0); indcount++) {
          if (chosen_r1indc == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovc1 = survind(indcount + delectable_sum);
          }
        }
      }
      
      double mainsum = rimeotam(survcoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, static_cast<double>(actualages(i)), inda1, inda2, indb1, indb2,
        indc1, indc2, dens, false);
      
      preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
        chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
        chosen_randcovc1 + survgroups2(0) + survgroups1(0) + 
        vital_patch(patchnumber, 0) + vital_year(yearnumber, 0) + surv_dev);

      if (preout > exp_tol) preout = exp_tol; // This catches numbers too high to be dealt with properly
      if (i < (noages - 1)) {
        survtransmat(i+1, i) = exp(preout) / (1.0 + exp(preout));
      } else {
        if (cont) {
          survtransmat(i, i) = exp(preout) / (1.0 + exp(preout));
        }
      }
    } else {
      if (i < (noages - 1)) {
        survtransmat(i+1, i) = survcoefs(0);
      } else {
        survtransmat(i, i) = survcoefs(0);
      }
    }
    
    // This next block calculates fecundity
    if (fec_addedcoefs != 0.0) {
      if (sf_repstatus(i) == 1) {
        
        double chosen_randcova2 {0.0};
        if (chosen_r2inda != "none") {
          for (int indcount = 0; indcount < rand_index(0, 6); indcount++) {
            if (chosen_r2inda == fecind_rownames(indcount)) {
              chosen_randcova2 = fecind(indcount);
            }
          }
        }
        double chosen_randcova1 {0.0};
        if (chosen_r1inda != "none") {
          int delectable_sum = rand_index(0, 6);
          for (int indcount = 0; indcount < rand_index(1, 6); indcount++) {
            if (chosen_r1inda == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcova1 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb2 {0.0};
        if (chosen_r2indb != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6);
          for (int indcount = 0; indcount < rand_index(2, 6); indcount++) {
            if (chosen_r2indb == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovb2 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb1 {0.0};
        if (chosen_r1indb != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6);
          for (int indcount = 0; indcount < rand_index(3, 6); indcount++) {
            if (chosen_r1indb == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovb1 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc2 {0.0};
        if (chosen_r2indc != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
            rand_index(3, 6);
          for (int indcount = 0; indcount < rand_index(4, 6); indcount++) {
            if (chosen_r2indc == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovc2 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc1 {0.0};
        if (chosen_r1indc != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
            rand_index(3, 6) + rand_index(4, 6);
          for (int indcount = 0; indcount < rand_index(5, 6); indcount++) {
            if (chosen_r1indc == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovc1 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcova2zi {0.0};
        if (chosen_r2inda != "none") {
          for (int indcount = 0; indcount < rand_index(0, 16); indcount++) {
            if (chosen_r2inda == fecind_rownames_zi(indcount)) {
              chosen_randcova2zi = fecindzi(indcount);
            }
          }
        }
        double chosen_randcova1zi {0.0};
        if (chosen_r1inda != "none") {
          int delectable_sum = rand_index(0, 16);
          for (int indcount = 0; indcount < rand_index(1, 16); indcount++) {
            if (chosen_r1inda == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcova1zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb2zi {0.0};
        if (chosen_r2indb != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16);
          for (int indcount = 0; indcount < rand_index(2, 16); indcount++) {
            if (chosen_r2indb == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovb2zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb1zi {0.0};
        if (chosen_r1indb != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16);
          for (int indcount = 0; indcount < rand_index(3, 16); indcount++) {
            if (chosen_r1indb == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovb1zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc2zi {0.0};
        if (chosen_r2indc != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
            rand_index(3, 16);
          for (int indcount = 0; indcount < rand_index(4, 16); indcount++) {
            if (chosen_r2indc == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovc2zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc1zi {0.0};
        if (chosen_r1indc != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
            rand_index(3, 16) + rand_index(4, 16);
          for (int indcount = 0; indcount < rand_index(5, 16); indcount++) {
            if (chosen_r1indc == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovc1zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
            
        double preoutx {0.0};
        
        if (fecdist < 4) {
          if (feczero) {
            
            double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, static_cast<double>(i), inda1, inda2, indb1, indb2, indc1,
              indc2, dens, true);
            
            preoutx = (mainsum + chosen_randcova2zi + chosen_randcova1zi +
              chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
              chosen_randcovc1zi + fecgroups2zi(0) + fecgroups1zi(0) + 
              fecpatchzi(patchnumber) + fecyearzi(yearnumber) + fec_dev);
            
          } else {
            
            double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, static_cast<double>(i), inda1, inda2, indb1, indb2, indc1,
              indc2, dens, false);
            
            preoutx = (mainsum + chosen_randcova2 + chosen_randcova1 +
              chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
              chosen_randcovc1 + fecgroups2(0) + fecgroups1(0) + 
              vital_patch(patchnumber, 1) + vital_year(yearnumber, 1) + fec_dev);
          }
          
          if (fecdist == 0 || fecdist == 1) {
            // Poisson and negative binomial fecundity
            
            if (feczero) {
              if (preoutx > exp_tol) preoutx = exp_tol;
              
              fectransmat(0, i) = (exp(preoutx) / (1.0 + exp(preoutx))) * fecmod;
            } else {
              if (preoutx > exp_tol) preoutx = exp_tol;
              
              fectransmat(0, i) = exp(preoutx) * fecmod;
            }
          } else if (fecdist == 2) {
            // Gaussian fecundity
            fectransmat(0, i) = preoutx * fecmod;
            
            if (negfec && fectransmat(0, i) < 0.0) {
              fectransmat(0, i) = 0.0;
            }
          } else if (fecdist == 3) {
            // Gamma fecundity
            fectransmat(0, i) = (1.0 / preoutx) * fecmod;
          }
          
        } else {
          fectransmat(0, i) = feccoefs(0);
        }
      }
    }
  }
  
  List output;
  
  if (simplicity) {
    output = List::create(_["U"] = survtransmat, _["F"] = fectransmat);
  } else {
    arma::mat amatrix = survtransmat + fectransmat;
    output = List::create(Named("A") = amatrix, _["U"] = survtransmat,
      _["F"] = fectransmat);
  }
  
  return output;
}

//' Extract Key Components from Simple Numerical Model
//' 
//' This function creates a skeleton list needed for functions
//' \code{jerzeibalowski()} and \code{motherbalowski()}, when a vital rate model
//' is simply a scalar.
//' 
//' @name numeric_extractor
//' 
//' @param object A numerical value, typical \code{1} or \code{0}.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}. Will generally be
//' \code{numeric}.}
//' \item{family}{The response distribution. Here, given as \code{constant}.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}. Here, given as \code{NULL}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.
//' Here, given as \code{NULL}.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}. Here, given as
//' \code{NULL}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables. Not used in \code{lm}/\code{glm}/\code{negbin} objects.
//' Here, given as \code{NULL}.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}. Here, given as \code{NULL}.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model. Here, given as \code{NULL}.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}. Here, given as \code{NULL}.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to
//' \code{1.0}.}
//' \item{theta}{The estimated theta, if the response is negative binomial.
//' Otherwise, will equal \code{1.0}.}
//' 
//' @keywords internal
//' @noRd
List numeric_extractor(NumericVector object) {
  String object_class = "numeric";
  String resp_family = "constant";
  int dist = 5;
  
  NumericVector coefs(1);
  coefs(0) = object(0);
  CharacterVector all_vars = {"Intercept"};
  
  double sigma {1.0};
  double theta {1.0};
  
  Rcpp::List output = List::create(_["class"] = object_class, _["family"] = resp_family,
    _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
    _["all_vars"] = all_vars, _["fixed_vars"] = all_vars, _["fixed_slopes"] = coefs,
    _["fixed_zi_vars"] = R_NilValue, _["fixed_zi_slopes"] = R_NilValue,
    _["random_vars"] = R_NilValue, _["random_slopes"] = R_NilValue,
    _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
    _["sigma"] = sigma, _["theta"] = theta);

  return output;
}

//' Extract Key Components of lm/glm/negbin Objects
//' 
//' This function extracts the components of an \code{lm}, \code{glm}, or
//' \code{negbin} (function \code{glm.nb()}) object needed for functions
//' \code{jerzeibalowski()} and \code{motherbalowski()}.
//' 
//' @name glm_extractor
//' 
//' @param object An \code{lm}, \code{glm}, or \code{negbin} object.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}. Will generally be either
//' \code{lm} or \code{glm}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated. Not used in \code{lm}/\code{glm}/\code{negbin} objects.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated. Not used in \code{lm}/\code{glm}/\code{negbin} objects.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables. Not used in \code{lm}/\code{glm}/\code{negbin} objects.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
//' objects.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model. Not used in \code{lm}/\code{glm}/
//' \code{negbin} objects.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
//' objects.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to
//' \code{1.0}.}
//' \item{theta}{The estimated theta, if the response is negative binomial.
//' Otherwise, will equal \code{1.0}.}
//' 
//' @section Notes:
//' Output from function \code{glm.nb()} is technically of class \code{negbin},
//' but is treated as class \code{glm} here.
//' 
//' @keywords internal
//' @noRd
List glm_extractor(List object) {
  StringVector input_class = object.attr("class");
  std::string object_class = as<std::string>(input_class[0]);
  String resp_family = "gaussian";
  int dist = 2;
  double theta = 1.0;
  
  if (stringcompare_hard(object_class, "glm")) {
    List big_resp = object["family"];
    resp_family = as<String>(big_resp["family"]);
    
    List str_check_pois = stringcompare_soft(resp_family, "poisson");
    List str_check_negb = stringcompare_soft(resp_family, "negbin");
    List str_check_gamm = stringcompare_soft(resp_family, "gamma");
    List str_check_bin = stringcompare_soft(resp_family, "binomial");
    if (str_check_pois["contains"]) {
      dist = 0;
    } else if (str_check_negb["contains"]) {
      dist = 1;
      theta = object["theta"];
    } else if (str_check_gamm["contains"]) {
      dist = 3;
    } else if (str_check_bin["contains"]) {
      dist = 4;
    }
  } else if (stringcompare_hard(object_class, "negbin")) {
    resp_family = "negbin";
    dist = 1;
    theta = object["theta"];
  }
  
  NumericVector coefs = object["coefficients"];
  CharacterVector all_vars = coefs.attr("names");
  
  NumericVector residuals = object["residuals"];
  int no_residuals = residuals.length();
  double sum_squared_residuals {0.0};
  double sigma {1.0};
  
  for (int i = 0; i < no_residuals; i++) {
    sum_squared_residuals += (residuals(i) * residuals(i));
  }
  if (no_residuals > 2) {
    sigma = sum_squared_residuals / (static_cast<double>(no_residuals) - 2.0);
    sigma = sqrt(sigma);
  }
  
  Rcpp::List output = List::create(_["class"] = object_class, _["family"] = resp_family,
    _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
    _["all_vars"] = all_vars, _["fixed_vars"] = all_vars, _["fixed_slopes"] = coefs,
    _["fixed_zi_vars"] = R_NilValue, _["fixed_zi_slopes"] = R_NilValue,
    _["random_vars"] = R_NilValue, _["random_slopes"] = R_NilValue,
    _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
    _["sigma"] = sigma, _["theta"] = theta);
  
  return output;
}

//' Extract Key Components of vglm Objects
//' 
//' This function extracts the components of a \code{vglm} object needed for
//' functions \code{jerzeibalowski()} and \code{motherbalowski()}.
//' 
//' @name vglm_extractor
//' 
//' @param object A \code{vglm} object.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated. Not used in \code{vglm} objects.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated. Always \code{TRUE} for \code{vglm} objects.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables. Not used in \code{vglm} objects.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}. Not used in \code{vglm} objects.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model. Not used in \code{vglm} objects.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}. Not used in \code{vglm} objects.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.}
//' \item{theta}{The estimated theta, if the response is negative binomial.
//' Otherwise, will equal \code{1.0}.}
//' 
//' @keywords internal
//' @noRd
List vglm_extractor(S4 object) {
  int dist {0};
  double theta = object.slot("dispersion");
  
  S4 m_family = object.slot("family");
  String resp_family = m_family.slot("vfamily");
  
  List str_check_pois = stringcompare_soft(resp_family, "poisson");
  List str_check_negb = stringcompare_soft(resp_family, "negbin");
  if (str_check_pois["contains"]) {
    dist = 0;
  } else if (str_check_negb["contains"]) {
    dist = 1;
  } else {
    throw Rcpp::exception("Response distribution not recognized.", false);
  }
  
  NumericVector fixed_slopes = object.slot("coefficients");
  CharacterVector fixed_vars = fixed_slopes.attr("names");
  
  List all_terms = object.slot("terms");
  List all_terms_stuff = all_terms["terms"];
  NumericMatrix ats = all_terms_stuff.attr("factors");
  CharacterVector all_vars = rownames(ats);
  
  Rcpp::List output = List::create(_["class"] = "vglm", _["family"] = resp_family,
    _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = true,
    _["all_vars"] = all_vars, _["fixed_vars"] = fixed_vars, _["fixed_slopes"] = fixed_slopes,
    _["fixed_zi_vars"] = R_NilValue, _["fixed_zi_slopes"] = R_NilValue,
    _["random_vars"] = R_NilValue, _["random_slopes"] = R_NilValue,
    _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
    _["sigma"] = 1.0, _["theta"] = theta);
  
  return output;
}

//' Extract Key Components of zeroinfl Objects
//' 
//' This function extracts the components of a \code{zeroinfl} object needed for
//' functions \code{jerzeibalowski()} and \code{motherbalowski()} to work.
//' 
//' @name zeroinfl_extractor
//' 
//' @param object A \code{zeroinfl} object.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}. Will always be
//' \code{glmmTMB}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated. Defaults to \code{FALSE} for \code{zeroinfl} objects.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}.}
//' \item{random_vars}{A string vector holding the names of the random
//' variables. Not used in \code{zeroinfl} objects.}
//' \item{random_slopes}{A numeric vector holding the slope coefficients of the
//' random variables, in the same order as \code{random_vars}. Not used in
//' \code{zeroinfl} objects.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model. Not used in \code{zeroinfl} objects.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}. Not used in \code{zeroinfl} objects.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
//' Equivalent output to lme4's \code{sigma()} function.}
//' \item{theta}{The scale parameter theta used in the negative binomial
//' distribution. Defaults to \code{1.0}.}
//' 
//' @section Notes:
//' This function will only work in the case where random terms are given as
//' \code{(1 | ranterm)}, where \code{ranterm} is the name of the random
//' variable.
//' 
//' @keywords internal
//' @noRd
List zeroinfl_extractor(List object) {
  String model_family = object["dist"];
  int dist {0};
  double theta {1.0};
  bool zi = true;
  
  List str_check_pois = stringcompare_soft(model_family, "poisson");
  List str_check_negb = stringcompare_soft(model_family, "negbin");
  
  if (str_check_pois["contains"]) {
    dist = 0;
    
  } else if (str_check_negb["contains"]) {
    dist = 1;
    theta = object["theta"];
    
  } else {
    throw Rcpp::exception("Unrecognized response distribution.", false);
  }
  
  List model_part = object["model"];
  CharacterVector all_vars = model_part.attr("names");
  
  List all_coefs = object["coefficients"];
  NumericVector fixed_slopes = all_coefs["count"];
  NumericVector zi_slopes = all_coefs["zero"];
  CharacterVector fixed_terms = fixed_slopes.attr("names");
  CharacterVector zi_terms = zi_slopes.attr("names");
  
  List output = List::create(_["class"] = "zeroinfl", _["family"] = model_family,
    _["dist"] = dist, _["zero_inflated"] = zi, _["zero_truncated"] = false,
    _["all_vars"] = all_vars, _["fixed_vars"] = fixed_terms, _["fixed_slopes"] = fixed_slopes,
    _["fixed_zi_vars"] = zi_terms, _["fixed_zi_slopes"] = zi_slopes,
    _["random_vars"] = R_NilValue, _["random_slopes"] = R_NilValue,
    _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
    _["sigma"] = 1.0, _["theta"] = theta);
  
  return output;
}

//' Extract Key Components of merMod Objects
//' 
//' This function extracts the components of a \code{merMod} object needed for
//' functions \code{jerzeibalowski()} and \code{motherbalowski()} to work.
//' 
//' @name lme4_extractor
//' 
//' @param object A \code{merMod} object.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}. Will generally be either
//' \code{lmerMod} or \code{glmerMod}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated. Not used in lme4.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated. Not used in lme4.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables. Not used in lme4 objects.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}. Not used in lme4 objects.}
//' \item{random_vars}{A string vector holding the names of the random
//' variables.}
//' \item{random_slopes}{A numeric vector holding the slope coefficients of the
//' random variables, in the same order as \code{random_vars}.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model. Not used in lme4.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}. Not used in lme4.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
//' Equivalent output to lme4's \code{sigma()} function.}
//' \item{theta}{Not used in lme4 output. Defaults to \code{1.0}.}
//' 
//' @section Notes:
//' This function will only work in the case where random terms are given as
//' \code{(1 | ranterm)}, where \code{ranterm} is the name of the random
//' variable.
//' 
//' @keywords internal
//' @noRd
List lme4_extractor(S4 object) {
  std::string object_class = object.attr("class");
  NumericVector coefs = object.slot("beta");
  String resp_family = "gaussian";
  int dist = 2;
  List pp = object.slot("pp");
  List cnms = object.slot("cnms");
  StringVector ran_names = cnms.attr("names");
  List flist = object.slot("flist");
  int no_ran_terms = flist.length();
  double sigma {1.0};
  double theta {1.0};
  
  if (stringcompare_hard(object_class, "glmerMod")) {
    List big_resp = object.slot("resp");
    List resp_family_S3 = big_resp["family"];
    resp_family = as<String>(resp_family_S3["family"]);
    
    List str_check_pois = stringcompare_soft(resp_family, "poisson");
    List str_check_negb = stringcompare_soft(resp_family, "negbin");
    List str_check_gamm = stringcompare_soft(resp_family, "gamma");
    List str_check_bin = stringcompare_soft(resp_family, "binomial");
    if (str_check_pois["contains"]) {
      dist = 0;
    } else if (str_check_negb["contains"]) {
      dist = 1;
    } else if (str_check_gamm["contains"]) {
      dist = 3;
    } else if (str_check_bin["contains"]) {
      dist = 4;
    }
  }
  
  // Names of all variables
  DataFrame oframe = as<DataFrame>(object.slot("frame"));
  StringVector all_var_names = oframe.attr("names");
  
  // Fixed variables
  NumericMatrix ppX = pp["X"];
  StringVector coef_names = colnames(ppX);
  
  // Random variables
  // "b" = crossprod(PR$Lambdat, object@u), # == Lambda %*% u
  arma::sp_mat lambdat_sp = as<sp_mat>(pp["Lambdat"]);
  arma::mat lambdat = arma::mat(lambdat_sp);
  NumericVector u = object.slot("u");
  arma::vec ucolvec = arma::vec(u);
  
  // b gives the random terms in the case where all random terms are (1 | ranterm)
  // We probably need to fill a matrix by row, and then get the col names as the random variables
  arma::mat b = lambdat * ucolvec;
  NumericVector ran_slopes(b.begin(), b.end());
  
  CharacterVector ran_term_index(ran_slopes.length());
  int ran_term_index_counter {0};
  List ran_term_list(no_ran_terms);
  List ran_index_list(no_ran_terms);
  
  // Names of random variables
  for (int i = 0; i < no_ran_terms; i++) {
    IntegerVector new_term = flist[i];
    CharacterVector new_term_names = new_term.attr("levels");
    
    int new_term_names_length = new_term_names.length();
    NumericVector new_value(new_term_names_length);
    
    for (int j = 0; j < new_term_names_length; j++) {
      //ran_term_index(j) = new_term_names(j);
      new_value(j) = ran_slopes(ran_term_index_counter);
      ran_term_index_counter++;
    }
    new_value.attr("names") = new_term_names;
    ran_term_list(i) = new_value;
    ran_index_list(i) = new_value.attr("names");
  }
  ran_index_list.attr("names") = ran_names;
  ran_term_list.attr("names") = ran_names;
  
  List devcomp = object.slot("devcomp");
  NumericVector cmp = devcomp["cmp"];
  IntegerVector dd = devcomp["dims"];
  CharacterVector cmp_names = cmp.attr("names");
  CharacterVector dd_names = dd.attr("names");
  int cmp_length = cmp.length();
  int dd_length = dd.length();
  int useSc_place = 0;
  int REML_place = 0;
  bool useSc = false;
  bool REML = false;
  
  for (int i = 0; i < dd_length; i++) {
    if (stringcompare_hard(as<std::string>(dd_names(i)), "useSc")) {
      useSc_place = i;
    } else if (stringcompare_hard(as<std::string>(dd_names(i)), "REML")) {
      REML_place = i;
    }
  }
  if (dd(useSc_place) > 0) {
    useSc = true;
    
    if (dd(REML_place) > 0) {
      REML = true;
    }
  }
  
  if (useSc) {
    for (int i = 0; i < cmp_length; i++) {
      if (REML) {
        if (stringcompare_hard(as<std::string>(cmp_names(i)), "sigmaREML")) {
          sigma = cmp(i);
        }
      } else {
        if (stringcompare_hard(as<std::string>(cmp_names(i)), "sigmaML")) {
          sigma = cmp(i);
        }
      }
    }
  }
  
  Rcpp::List output = List::create(_["class"] = object_class, _["family"] = resp_family,
    _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
    _["all_vars"] = all_var_names, _["fixed_vars"] = coef_names, _["fixed_slopes"] = coefs,
    _["fixed_zi_vars"] = R_NilValue, _["fixed_zi_slopes"] = R_NilValue,
    _["random_vars"] = ran_index_list, _["random_slopes"] = ran_term_list,
    _["random_zi_vars"] = R_NilValue, _["random_zi_slopes"] = R_NilValue,
    _["sigma"] = sigma, _["theta"] = theta);
  
  return output;
}

//' Extract Key Components of glmmTMB Objects
//' 
//' This function extracts the components of a \code{glmmTMB} object needed for
//' functions \code{jerzeibalowski()} and \code{motherbalowski()} to work.
//' 
//' @name glmmTMB_extractor
//' 
//' @param object A \code{glmmTMB} object.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}. Will always be
//' \code{glmmTMB}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}.}
//' \item{random_vars}{A string vector holding the names of the random
//' variables.}
//' \item{random_slopes}{A numeric vector holding the slope coefficients of the
//' random variables, in the same order as \code{random_vars}.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
//' Equivalent output to lme4's \code{sigma()} function.}
//' \item{theta}{The scale parameter theta used in the negative binomial
//' distribution. Defaults to \code{1.0}.}
//' 
//' @section Notes:
//' This function will only work in the case where random terms are given as
//' \code{(1 | ranterm)}, where \code{ranterm} is the name of the random
//' variable.
//' 
//' @keywords internal
//' @noRd
List glmmTMB_extractor(List object) {
  
  List the_green_flist;
  List the_green_ziflist;
  CharacterVector ran_zi_vars;
  
  List dissident_aggressor = object["fit"];
  NumericVector parameter_values = dissident_aggressor["par"];
  CharacterVector term_names = parameter_values.attr("names");
  
  NumericVector all_the_stuff = dissident_aggressor["parfull"];
  CharacterVector names_of_all_the_stuff = all_the_stuff.attr("names");
  int length_of_all_the_stuff = all_the_stuff.length();
  
  String class_thistime = {"glmmTMB"};
  int dist = 2;
  double sigma {1.0};
  double theta {1.0};
  
  List variance_crap = object["sdr"];
  NumericVector random_values = variance_crap["par.random"];
  CharacterVector random_values_tags = random_values.attr("names");
  int no_random_values = random_values.length(); // This has ALL random coefficients
  LogicalVector random_b(no_random_values);
  LogicalVector random_bzi(no_random_values);
  int no_ranb {0};
  int no_ranbzi {0};
  int no_ran_vars {0};
  int no_rz_vars {0};  
  
  for (int i = 0; i < no_random_values; i++) {
    if (stringcompare_hard(as<std::string>(random_values_tags(i)), "b")) {
      no_ranb++;
      random_b(i) = 1;
    } else if (stringcompare_hard(as<std::string>(random_values_tags(i)), "bzi")) {
      no_ranbzi++;
      random_bzi(i) = 1;
    }
  }
  
  DataFrame exploding_reloading = as<DataFrame>(object["frame"]);
  CharacterVector all_vars = exploding_reloading.attr("names");
  
  List the_green_manalishi = object["modelInfo"];
  CharacterVector ran_vars = the_green_manalishi["grpVar"];
  List the_green_reTrms = the_green_manalishi["reTrms"];
  List the_green_cond = the_green_reTrms["cond"];
  List the_green_zi = the_green_reTrms["zi"];
  
  if (no_ranb > 0) {
    the_green_flist = the_green_cond["flist"];
    no_ran_vars = the_green_flist.length();
  }
  if (no_ranbzi > 0) {
    the_green_ziflist = the_green_zi["flist"];
    no_rz_vars = the_green_ziflist.length();
    ran_zi_vars = the_green_ziflist.attr("names");
  }
  
  List obj = object["obj"];
  Environment env = obj["env"];
  List data = env["data"];
  NumericMatrix Xmat = data["X"];
  CharacterVector all_fixed_terms = colnames(Xmat);
  NumericMatrix Xmatzi = data["Xzi"];
  CharacterVector all_zi_terms = colnames(Xmatzi);
  
  int total_fixed_slopes {0};
  int total_zi_slopes {0};
  
  // Extract random slopes
  List ran_term_list(no_ran_vars);
  List ran_slope_list(no_ran_vars);
  List ran_zi_term_list;
  List ran_zi_slope_list;
  int ran_term_counter {0};
  
  // Extract cond random coefficients
  if (no_ranb > 0) {
    for (int i = 0; i < no_ran_vars; i++) {
      
      IntegerVector current_factor = the_green_flist[i];
      CharacterVector current_names = current_factor.attr("levels");
      
      ran_term_list(i) = current_names;
      int current_names_length = current_names.length();
      NumericVector current_slopes(current_names_length);
      
      for (int j = 0; j < current_names_length; j++) {
        if (random_b(ran_term_counter) > 0) {
          current_slopes(j) = random_values(ran_term_counter);
        }
        ran_term_counter++;
      }
      ran_slope_list(i) = current_slopes;
    }
    ran_term_list.attr("names") = ran_vars;
    ran_slope_list.attr("names") = ran_vars;
  }
  
  // Extract zi random coefficients
  if (no_ranbzi > 0) {
    List rztl(no_rz_vars);
    List rstl(no_rz_vars);
    
    for (int i = 0; i < no_rz_vars; i++) {
      IntegerVector current_factor = the_green_ziflist[i];
      CharacterVector current_names = current_factor.attr("levels");
      
      rztl(i) = current_names;
      int current_names_length = current_names.length();
      NumericVector current_slopes(current_names_length);
      
      for (int j = 0; j < current_names_length; j++) {
        if (random_bzi(ran_term_counter) > 0) {
          current_slopes(j) = random_values(ran_term_counter);
        }
        ran_term_counter++;
      }
      rstl(i) = current_slopes;
    }
    
    ran_zi_term_list = rztl;
    ran_zi_slope_list = rstl;
    ran_zi_term_list.attr("names") = ran_zi_vars;
    ran_zi_slope_list.attr("names") = ran_zi_vars;
  }
  
  // Extract dispersion parameter
  double dispersion {1.0};
  for (int i = 0; i < length_of_all_the_stuff; i++) {
    if (stringcompare_hard(as<std::string>(names_of_all_the_stuff(i)), "betad")) {
      dispersion = all_the_stuff(i);
    }
  }
  
  // Response distribution tests
  bool trunc = false;
  bool zi = false;
  
  List with_the_two_pronged_crown = the_green_manalishi["family"];
  String model_family = as<String>(with_the_two_pronged_crown["family"]);
  
  List str_check_pois = stringcompare_soft(model_family, "poisson");
  List str_check_negb = stringcompare_soft(model_family, "negbin");
  List str_check_nbin = stringcompare_soft(model_family, "nbinom");
  List str_check_gamm = stringcompare_soft(model_family, "Gamma");
  List str_check_bin = stringcompare_soft(model_family, "binomial");
  List str_check_trun = stringcompare_soft(model_family, "trunc");
  
  if (str_check_pois["contains"]) {
    dist = 0;
    
    if (str_check_trun["contains"]) {
      trunc = true;
    }
  } else if (str_check_negb["contains"] || str_check_nbin["contains"]) {
    dist = 1;
    theta = exp(dispersion);
    
    if (str_check_trun["contains"]) {
      trunc = true;
    }
  } else if (str_check_gamm["contains"]) {
    dist = 3;
    sigma = exp(-0.5 * dispersion);
    
    List str_check_zi = stringcompare_soft(model_family, "zi");
    if (str_check_zi["contains"]) {
      zi = true;
    }
  } else if (str_check_bin["contains"]) {
    dist = 4;
  } else if (dist == 2) {
    sigma = exp(0.5 * dispersion);
  } else {
    throw Rcpp::exception("Unrecognized response distribution.", false);
  }
  
  // Most zero-inflation tests (apart from ziGamma)
  int no_params = parameter_values.length();
  
  for (int i = 0; i < no_params; i++) {
    if (term_names(i) == "betazi" && zi == false) {
      //model_family += zi_addition;
      zi = true;
    }
    
    if (term_names(i) == "beta") {
      
      total_fixed_slopes++;
    } else if (term_names(i) == "betazi") {
      total_zi_slopes++;
    }
  }
  
  // Extract fixed and zero-inflated slopes
  NumericVector fixed_slopes(total_fixed_slopes);
  NumericVector zi_slopes(total_zi_slopes);
  CharacterVector fixed_terms(total_fixed_slopes);
  CharacterVector zi_terms(total_zi_slopes);
  int fixed_slope_counter {0};
  int zi_slope_counter {0};
  
  for (int i = 0; i < no_params; i++) {
    if (term_names(i) == "beta") {
      fixed_slopes(fixed_slope_counter) = parameter_values(i);
      fixed_terms(fixed_slope_counter) = all_fixed_terms(fixed_slope_counter);
      fixed_slope_counter++;
    } else if (term_names(i) == "betazi") {
      zi_slopes(zi_slope_counter) = parameter_values(i);
      zi_terms(zi_slope_counter) = all_zi_terms(zi_slope_counter);
      zi_slope_counter++;
    }
  }
  
  List output = List::create(_["class"] = class_thistime, _["family"] = model_family,
    _["dist"] = dist, _["zero_inflated"] = zi, _["zero_truncated"] = trunc,
    _["all_vars"] = all_vars, _["fixed_vars"] = fixed_terms, _["fixed_slopes"] = fixed_slopes,
    _["fixed_zi_vars"] = zi_terms, _["fixed_zi_slopes"] = zi_slopes,
    _["random_vars"] = ran_term_list, _["random_slopes"] = ran_slope_list,
    _["random_zi_vars"] = ran_zi_term_list, _["random_zi_slopes"] = ran_zi_slope_list,
    _["sigma"] = sigma, _["theta"] = theta);
  
  return output;
}

//' Function Extracting Core Components From S3 Vital Rate Models
//' 
//' Function \code{S3_extractor()} extracts all needed terms from S3 objects
//' used as vital rate models.
//' 
//' @name S3_extractor
//' 
//' @param object An S3 vital rate model. Currently, this should be output from
//' functions \code{lm()}, \code{glm()}, \code{glm.nb()}, \code{zeroinfl()},
//' and \code{glmmTMB()}.
//' 
//' @return A list describing the vital rate model in standard output required
//' from function \code{modelextract()} to operate. Elements currently include:
//' \item{class}{The exact class of \code{object}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}.}
//' \item{random_vars}{A string vector holding the names of the random
//' variables.}
//' \item{random_slopes}{A numeric vector holding the slope coefficients of the
//' random variables, in the same order as \code{random_vars}.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
//' Equivalent output to lme4's \code{sigma()} function.}
//' \item{theta}{The scale parameter theta used in the negative binomial
//' distribution. Defaults to \code{1.0}.}
//' 
//' @section Notes:
//' This function currently handles models developed with functions \code{lm()}
//' and \code{glm()} from package \code{stats}, function \code{glm.nb()} from
//' package \code{MASS}, function \code{zeroinfl()} from package \code{pscl},
//' and function \code{glmmTMB()} from package \code{glmmTMB}.
//' 
//' @keywords internal
//' @noRd
List S3_extractor(List object) {
  StringVector model_class = object.attr("class");
  int model_type {0}; // 0 = unrecognized, 1 = lm/glm/negbin, 2 = zeroinfl, 3 = glmmTMB
  
  List output;
  
  for (int i = 0; i < model_class.length(); i++) {
    if (stringcompare_hard(as<std::string>(model_class(i)), "lm")) {
      model_type = 1;
    } else if (stringcompare_hard(as<std::string>(model_class(i)), "zeroinfl")) {
      model_type = 2;
    } else if (stringcompare_hard(as<std::string>(model_class(i)), "glmmTMB")) {
      model_type = 3;
    }
  }
  
  if (model_type == 1) {
    output = glm_extractor(object);
  } else if (model_type == 2) {
    output = zeroinfl_extractor(object);
  } else if (model_type == 3) {
    output = glmmTMB_extractor(object);
  } else {
    throw Rcpp::exception("Model type unrecognized.", false);
  }
  
  return output;
}

//' Function Extracting Core Components From S4 Vital Rate Models
//' 
//' Function \code{S4_extractor()} extracts all needed terms from S4 objects
//' used as vital rate models.
//' 
//' @name S4_extractor
//' 
//' @param object An S4 vital rate model. Currently, this should be output from
//' functions \code{vglm()}, \code{lmer()}, and \code{glmer()}.
//' 
//' @return A list describing the vital rate model in standard output required
//' from function \code{modelextract()} to operate. Elements currently include:
//' \item{class}{The exact class of \code{object}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}.}
//' \item{random_vars}{A string vector holding the names of the random
//' variables.}
//' \item{random_slopes}{A numeric vector holding the slope coefficients of the
//' random variables, in the same order as \code{random_vars}.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the sero-inflation model.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to 1.0.
//' Equivalent output to lme4's \code{sigma()} function.}
//' \item{theta}{The scale parameter theta used in the negative binomial
//' distribution. Defaults to \code{1.0}.}
//' 
//' @section Notes:
//' This function currently handles models developed with function \code{vglm()}
//' from package \code{VGAM}, and functions \code{lmer()} and \code{glmer()}
//' from package \code{lme4}.
//' 
//' @keywords internal
//' @noRd
List S4_extractor(S4 object) {
  String model_class = object.attr("class");
  int model_type {0}; // 0 = unrecognized, 1 = vglm, 2 = merMod
  
  List output;
  
  if (stringcompare_hard(model_class, "vglm")) {
    model_type = 1;
    output = vglm_extractor(object);
  } else if (stringcompare_hard(model_class, "lmerMod") || 
      stringcompare_hard(model_class, "glmerMod")) {
    model_type = 2;
    output = lme4_extractor(object);
  } else {
    throw Rcpp::exception("Model type unrecognized.", false);
  }
  
  return output;
}

//' Extract Coefficients From Linear Vital Rate Models
//' 
//' Function \code{modelextract()} extracts coefficient values from linear
//' models estimated through various linear modeling functions in R, to estimate
//' vital rates in \code{lefko3}. Used to supply coefficients to
//' \code{\link{flefko3}()}, \code{\link{flefko2}()}, \code{\link{fleslie()}},
//' and \code{\link{aflefko2}()}.
//' 
//' @param object A linear model estimated through one of the methods used in
//' function \code{\link{modelsearch}()}.
//' @param paramnames Data frame giving the names of standard coefficients
//' required by matrix creation functions.
//' @param mainyears A vector of the names of the monitoring occasions.
//' @param mainpatches A vector of the names of the patches. Should be \code{NA}
//' if no patches specified.
//' @param maingroups A vector of the names of all stage groups.
//' @param mainindcova A vector denoting values of individual covariate
//' \code{a}, when that individual covariate is categorical.
//' @param mainindcovb A vector denoting values of individual covariate
//' \code{b}, when that individual covariate is categorical.
//' @param mainindcovc A vector denoting values of individual covariate
//' \code{c}, when that individual covariate is categorical.
//' 
//' @return This function returns a list with the following elements:
//' \item{coefficients}{Vector of fixed effect coefficients.}
//' \item{years}{Vector of occasion coefficients, typically random.}
//' \item{zeroyear}{Vector of zero-inflated occasion coefficients, typically
//' random.}
//' \item{patches}{Vector of patch coefficients, typically random.}
//' \item{zeropatch}{Vector of zero-inflated patch coefficients, typically
//' random.}
//' \item{groups2}{Vector of group coefficients for time t.}
//' \item{groups1}{Vector of group coefficients for time t-1.}
//' \item{zerogroups2}{Vector of zero-inflated group coefficients for time t.}
//' \item{zerogroups1}{Vector of zero-inflated group coefficients for time t-1.}
//' \item{indcova2s}{Vector of individual covariate \code{a} values for time t.}
//' \item{indcova1s}{Vector of individual covariate \code{a} values for time t-1.}
//' \item{indcovb2s}{Vector of individual covariate \code{b} values for time t.}
//' \item{indcovb1s}{Vector of individual covariate \code{b} values for time t-1.}
//' \item{indcovc2s}{Vector of individual covariate \code{c} values for time t.}
//' \item{indcovc1s}{Vector of individual covariate \code{c} values for time t-1.}
//' \item{zeroindcova2s}{Vector of zero-inflated individual covariate \code{a}
//' values for time t.}
//' \item{zeroindcova1s}{Vector of zero-inflated individual covariate \code{a}
//' values for time t-1.}
//' \item{zeroindcovb2s}{Vector of zero-inflated individual covariate \code{b}
//' values for time t.}
//' \item{zeroindcovb1s}{Vector of zero-inflated individual covariate \code{b}
//' values for time t-1.}
//' \item{zeroindcovc2s}{Vector of zero-inflated individual covariate \code{c}
//' values for time t.}
//' \item{zeroindcovc1s}{Vector of zero-inflated individual covariate \code{c}
//' values for time t-1.}
//' \item{class}{The R class of the vital rate model.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to
//' \code{1.0}. Equivalent output to package lme4's \code{sigma()} function.}
//' \item{theta}{The scale parameter theta used in the negative binomial
//' distribution. Defaults to \code{1.0}.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.modelextract)]]
List modelextract(RObject object, DataFrame paramnames, NumericVector mainyears,
  CharacterVector mainpatches, RObject maingroups, RObject mainindcova,
  RObject mainindcovb, RObject mainindcovc) {
  
  CharacterVector fixed_zi_vars;
  NumericVector fixed_zi_slopes;
  List random_vars;
  List random_slopes;
  List random_zi_vars;
  List random_zi_slopes;
  
  NumericVector indcova2s;
  NumericVector indcova1s;
  NumericVector indcovb2s;
  NumericVector indcovb1s;
  NumericVector indcovc2s;
  NumericVector indcovc1s;
  NumericVector zeroindcova2s;
  NumericVector zeroindcova1s;
  NumericVector zeroindcovb2s;
  NumericVector zeroindcovb1s;
  NumericVector zeroindcovc2s;
  NumericVector zeroindcovc1s;
  
  List core_components;
  if (object.isS4()) {
    core_components = S4_extractor(as<S4>(object));
  } else if (is<List>(object)) {
    core_components = S3_extractor(as<List>(object));
  } else {
    core_components = numeric_extractor(as<NumericVector>(object));
  }
  
  CharacterVector fixed_vars = core_components["fixed_vars"];
  NumericVector fixed_slopes = core_components["fixed_slopes"];
  int no_fixed_vars = fixed_vars.length();
  
  Nullable<CharacterVector> fixed_zi_vars_ = core_components["fixed_zi_vars"];
  Nullable<NumericVector> fixed_zi_slopes_ = core_components["fixed_zi_slopes"];
  
  if (fixed_zi_vars_.isNotNull()) {
    fixed_zi_vars = fixed_zi_vars_;
  } else {
    fixed_zi_vars = {"Intercept"};
  }
  if (fixed_zi_slopes_.isNotNull()) {
    fixed_zi_slopes = fixed_zi_slopes_;
  } else {
    fixed_zi_slopes = {0.0};
  }
  
  CharacterVector modelparam_names = paramnames["modelparams"];
  std::string year2var = as<std::string>(modelparam_names(0));
  std::string individvar = as<std::string>(modelparam_names(1));
  std::string patchvar = as<std::string>(modelparam_names(2));
  std::string surv3var = as<std::string>(modelparam_names(3));
  std::string obs3var = as<std::string>(modelparam_names(4));
  std::string size3var = as<std::string>(modelparam_names(5));
  std::string sizeb3var = as<std::string>(modelparam_names(6));
  std::string sizec3var = as<std::string>(modelparam_names(7));
  std::string repst3var = as<std::string>(modelparam_names(8));
  std::string fec3var = as<std::string>(modelparam_names(9));
  std::string fec2var = as<std::string>(modelparam_names(10));
  std::string size2var = as<std::string>(modelparam_names(11));
  std::string size1var = as<std::string>(modelparam_names(12));
  std::string sizeb2var = as<std::string>(modelparam_names(13));
  std::string sizeb1var = as<std::string>(modelparam_names(14));
  std::string sizec2var = as<std::string>(modelparam_names(15));
  std::string sizec1var = as<std::string>(modelparam_names(16));
  std::string repst2var = as<std::string>(modelparam_names(17));
  std::string repst1var = as<std::string>(modelparam_names(18));
  std::string matst3var = as<std::string>(modelparam_names(19));
  std::string matst2var = as<std::string>(modelparam_names(20));
  std::string agevar = as<std::string>(modelparam_names(21));
  std::string densityvar = as<std::string>(modelparam_names(22));
  std::string indcova2var = as<std::string>(modelparam_names(23));
  std::string indcova1var = as<std::string>(modelparam_names(24));
  std::string indcovb2var = as<std::string>(modelparam_names(25));
  std::string indcovb1var = as<std::string>(modelparam_names(26));
  std::string indcovc2var = as<std::string>(modelparam_names(27));
  std::string indcovc1var = as<std::string>(modelparam_names(28));
  std::string group2var = as<std::string>(modelparam_names(29));
  std::string group1var = as<std::string>(modelparam_names(30));
  
  int no_fixed_zi_vars = fixed_zi_vars.length();
  
  int no_years = mainyears.length();
  CharacterVector mainyears_text(mainyears);
  NumericVector year_coefs(no_years);
  NumericVector zeroyear_coefs(no_years);
  year_coefs.attr("names") = mainyears_text;
  zeroyear_coefs.attr("names") = mainyears_text;
  
  int no_patches = mainpatches.length();
  NumericVector patch_coefs(no_patches);
  NumericVector zeropatch_coefs(no_patches);
  if (no_patches < 2 && Rcpp::traits::is_na<STRSXP>(mainpatches(0))) {
    CharacterVector new_patch_names = {"pop"};
    patch_coefs.attr("names") = new_patch_names;
    zeropatch_coefs.attr("names") = new_patch_names;
  } else {
    patch_coefs.attr("names") = mainpatches;
    zeropatch_coefs.attr("names") = mainpatches;
  }
  
  CharacterVector maingroups_text(maingroups);
  int no_groups = maingroups_text.length();
  NumericVector group2_coefs(no_groups);
  NumericVector group1_coefs(no_groups);
  NumericVector zerogroup2_coefs(no_groups);
  NumericVector zerogroup1_coefs(no_groups);
  group2_coefs.attr("names") = maingroups_text;
  group1_coefs.attr("names") = maingroups_text;
  zerogroup2_coefs.attr("names") = maingroups_text;
  zerogroup1_coefs.attr("names") = maingroups_text;
  
  CharacterVector indcova_names;
  CharacterVector indcovb_names;
  CharacterVector indcovc_names;
  int no_indcova_names {0};
  int no_indcovb_names {0};
  int no_indcovc_names {0};
  
  if (is<CharacterVector>(mainindcova) || is<NumericVector>(mainindcova)) {
    indcova_names = as<CharacterVector>(mainindcova);
    no_indcova_names = indcova_names.length();
    
    bool realnames = false;
    if (no_indcova_names > 0) {
      for (int i = 0; i < no_indcova_names; i++) {
        if (!stringcompare_hard(as<std::string>(indcova_names(i)), "none")) realnames = true;
      }
    }
    
    if (!realnames) {
      no_indcova_names = 0;
    } else {
      NumericVector indcova2s_inc(no_indcova_names);
      NumericVector indcova1s_inc(no_indcova_names);
      NumericVector zeroindcova2s_inc(no_indcova_names);
      NumericVector zeroindcova1s_inc(no_indcova_names);
      
      indcova2s = indcova2s_inc;
      indcova1s = indcova1s_inc;
      zeroindcova2s = zeroindcova2s_inc;
      zeroindcova1s = zeroindcova1s_inc;
      
      indcova2s.attr("names") = indcova_names;
      indcova1s.attr("names") = indcova_names;
      zeroindcova2s.attr("names") = indcova_names;
      zeroindcova1s.attr("names") = indcova_names;
    }
  }
  
  if (is<CharacterVector>(mainindcovb) || is<NumericVector>(mainindcovb)) {
    indcovb_names = as<CharacterVector>(mainindcovb);
    no_indcovb_names = indcovb_names.length();
    
    bool realnames = false;
    if (no_indcovb_names > 0) {
      for (int i = 0; i < no_indcovb_names; i++) {
        if (!stringcompare_hard(as<std::string>(indcovb_names(i)), "none")) realnames = true;
      }
    }
    
    if (!realnames) {
      no_indcovb_names = 0;
    } else {
      NumericVector indcovb2s_inc(no_indcovb_names);
      NumericVector indcovb1s_inc(no_indcovb_names);
      NumericVector zeroindcovb2s_inc(no_indcovb_names);
      NumericVector zeroindcovb1s_inc(no_indcovb_names);
      
      indcovb2s = indcovb2s_inc;
      indcovb1s = indcovb1s_inc;
      zeroindcovb2s = zeroindcovb2s_inc;
      zeroindcovb1s = zeroindcovb1s_inc;
      
      indcovb2s.attr("names") = indcovb_names;
      indcovb1s.attr("names") = indcovb_names;
      zeroindcovb2s.attr("names") = indcovb_names;
      zeroindcovb1s.attr("names") = indcovb_names;
    }
  }
  
  if (is<CharacterVector>(mainindcovc) || is<NumericVector>(mainindcovc)) {
    indcovc_names = as<CharacterVector>(mainindcovc);
    no_indcovc_names = indcovc_names.length();
    
    bool realnames = false;
    if (no_indcovc_names > 0) {
      for (int i = 0; i < no_indcovc_names; i++) {
        if (!stringcompare_hard(as<std::string>(indcovc_names(i)), "none")) realnames = true;
      }
    }
    
    if (!realnames) {
      no_indcovc_names = 0;
    } else {
      NumericVector indcovc2s_inc(no_indcovc_names);
      NumericVector indcovc1s_inc(no_indcovc_names);
      NumericVector zeroindcovc2s_inc(no_indcovc_names);
      NumericVector zeroindcovc1s_inc(no_indcovc_names);
      
      indcovc2s = indcovc2s_inc;
      indcovc1s = indcovc1s_inc;
      zeroindcovc2s = zeroindcovc2s_inc;
      zeroindcovc1s = zeroindcovc1s_inc;
      
      indcovc2s.attr("names") = indcovc_names;
      indcovc1s.attr("names") = indcovc_names;
      zeroindcovc2s.attr("names") = indcovc_names;
      zeroindcovc1s.attr("names") = indcovc_names;
    }
  }
  
  NumericVector coef_vec(283);
  
  for (int i = 0; i < no_fixed_vars; i++) {
    for (int j = 0; j < no_years; j++) {
      if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(mainyears_text(j)))) {
        year_coefs(j) = fixed_slopes(i);
      }
    }
    for (int j = 0; j < no_patches; j++) {
      if (stringcompare_simple(as<std::string>(fixed_vars(i)), patchvar)) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(mainpatches(j)))) {
          patch_coefs(j) = fixed_slopes(i);
        }
      }
    }
    
    for (int j = 0; j < no_groups; j++) {
      if (stringcompare_simple(as<std::string>(fixed_vars(i)), group2var)) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(maingroups_text(j)))) {
          group2_coefs(j) = fixed_slopes(i);
        }
      }
    }
    for (int j = 0; j < no_groups; j++) {
      if (stringcompare_simple(as<std::string>(fixed_vars(i)), group1var)) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(maingroups_text(j)))) {
          group1_coefs(j) = fixed_slopes(i);
        }
      }
    }
    
    if (no_indcova_names > 0) {
      for (int j = 0; j < no_indcova_names; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcova2var)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(indcova_names(j)))) {
            indcova2s(j) = fixed_slopes(i);
          }
        }
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcova1var)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(indcova_names(j)))) {
            indcova1s(j) = fixed_slopes(i);
          }
        }
      }
    }
    
    if (no_indcovb_names > 0) {
      for (int j = 0; j < no_indcovb_names; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovb2var)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(indcovb_names(j)))) {
            indcovb2s(j) = fixed_slopes(i);
          }
        }
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovb1var)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(indcovb_names(j)))) {
            indcovb1s(j) = fixed_slopes(i);
          }
        }
      }
    }
    
    if (no_indcovc_names > 0) {
      for (int j = 0; j < no_indcovc_names; j++) {
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovc2var)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(indcovc_names(j)))) {
            indcovc2s(j) = fixed_slopes(i);
          }
        }
        if (stringcompare_simple(as<std::string>(fixed_vars(i)), indcovc1var)) {
          if (stringcompare_simple(as<std::string>(fixed_vars(i)), as<std::string>(indcovc_names(j)))) {
            indcovc1s(j) = fixed_slopes(i);
          }
        }
      }
    }
    
    if (stringcompare_simple(as<std::string>(fixed_vars(i)), "ntercep")) {
      coef_vec(0) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), repst1var)) {
      coef_vec(1) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), repst2var)) {
      coef_vec(2) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), size1var)) {
      coef_vec(3) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), size2var)) {
      coef_vec(4) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), repst1var, repst2var)) {
      coef_vec(5) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, size2var)) {
      coef_vec(6) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, repst1var)) {
      coef_vec(7) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, repst2var)) {
      coef_vec(8) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, repst1var)) {
      coef_vec(9) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, repst2var)) {
      coef_vec(10) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), agevar)) {
      coef_vec(11) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, agevar)) {
      coef_vec(12) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, agevar)) {
      coef_vec(13) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), repst1var, agevar)) {
      coef_vec(14) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), repst2var, agevar)) {
      coef_vec(15) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcova2var)) {
      coef_vec(16) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovb2var)) {
      coef_vec(17) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovc2var)) {
      coef_vec(18) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcova1var)) {
      coef_vec(19) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovb1var)) {
      coef_vec(20) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), indcovc1var)) {
      coef_vec(21) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, size2var)) {
      coef_vec(22) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, size2var)) {
      coef_vec(23) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, size2var)) {
      coef_vec(24) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, repst2var)) {
      coef_vec(25) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, repst2var)) {
      coef_vec(26) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, repst2var)) {
      coef_vec(27) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, size1var)) {
      coef_vec(28) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, size1var)) {
      coef_vec(29) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, size1var)) {
      coef_vec(30) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, repst1var)) {
      coef_vec(31) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, repst1var)) {
      coef_vec(32) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, repst1var)) {
      coef_vec(33) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovb2var)) {
      coef_vec(34) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovc2var)) {
      coef_vec(35) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, indcovc2var)) {
      coef_vec(36) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovb1var)) {
      coef_vec(37) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovc1var)) {
      coef_vec(38) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, indcovc1var)) {
      coef_vec(39) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovb1var)) {
      coef_vec(40) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovb2var)) {
      coef_vec(41) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, indcovc1var)) {
      coef_vec(42) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, indcovc2var)) {
      coef_vec(43) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, indcovc1var)) {
      coef_vec(44) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, indcovc2var)) {
      coef_vec(45) = fixed_slopes(i);
    }
    
    // New coefficients
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizeb2var)) {
      coef_vec(100) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizeb1var)) {
      coef_vec(101) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizec2var)) {
      coef_vec(102) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), sizec1var)) {
      coef_vec(103) = fixed_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_vars(i)), densityvar)) {
      coef_vec(104) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, sizeb2var)) {
      coef_vec(105) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, sizec2var)) {
      coef_vec(106) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizeb1var)) {
      coef_vec(107) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizec1var)) {
      coef_vec(108) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, sizec1var)) {
      coef_vec(109) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizeb2var)) {
      coef_vec(110) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizec2var)) {
      coef_vec(111) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, sizec2var)) {
      coef_vec(112) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizeb2var)) {
      coef_vec(113) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, sizec2var)) {
      coef_vec(114) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, sizec2var)) {
      coef_vec(115) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizeb1var)) {
      coef_vec(116) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, sizec1var)) {
      coef_vec(117) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, sizec1var)) {
      coef_vec(118) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size2var, densityvar)) {
      coef_vec(119) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, densityvar)) {
      coef_vec(120) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, densityvar)) {
      coef_vec(121) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), size1var, densityvar)) {
      coef_vec(122) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, densityvar)) {
      coef_vec(123) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, densityvar)) {
      coef_vec(124) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), repst2var, densityvar)) {
      coef_vec(125) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), repst1var, densityvar)) {
      coef_vec(126) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, repst2var)) {
      coef_vec(127) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, repst2var)) {
      coef_vec(128) = fixed_slopes(i);
    }
    // 129 = 0
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, repst1var)) {
      coef_vec(130) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, repst1var)) {
      coef_vec(131) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, repst2var)) {
      coef_vec(132) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, repst1var)) {
      coef_vec(133) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, repst1var)) {
      coef_vec(134) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, repst2var)) {
      coef_vec(135) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb2var, agevar)) {
      coef_vec(136) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec2var, agevar)) {
      coef_vec(137) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), densityvar, agevar)) {
      coef_vec(138) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizeb1var, agevar)) {
      coef_vec(139) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), sizec1var, agevar)) {
      coef_vec(140) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizeb2var)) {
      coef_vec(141) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizec2var)) {
      coef_vec(142) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, densityvar)) {
      coef_vec(143) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizeb1var)) {
      coef_vec(144) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizec1var)) {
      coef_vec(145) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizeb2var)) {
      coef_vec(146) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, sizec2var)) {
      coef_vec(147) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizeb1var)) {
      coef_vec(148) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, sizec1var)) {
      coef_vec(149) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, densityvar)) {
      coef_vec(150) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizeb2var)) {
      coef_vec(151) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizec2var)) {
      coef_vec(152) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, densityvar)) {
      coef_vec(153) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizeb1var)) {
      coef_vec(154) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizec1var)) {
      coef_vec(155) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizeb2var)) {
      coef_vec(156) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, sizec2var)) {
      coef_vec(157) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizeb1var)) {
      coef_vec(158) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, sizec1var)) {
      coef_vec(159) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, densityvar)) {
      coef_vec(160) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizeb2var)) {
      coef_vec(161) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizec2var)) {
      coef_vec(162) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, densityvar)) {
      coef_vec(163) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizeb1var)) {
      coef_vec(164) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizec1var)) {
      coef_vec(165) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizeb2var)) {
      coef_vec(166) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, sizec2var)) {
      coef_vec(167) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizeb1var)) {
      coef_vec(168) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, sizec1var)) {
      coef_vec(169) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, densityvar)) {
      coef_vec(170) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, size1var)) {
      coef_vec(171) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, size1var)) {
      coef_vec(172) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, size1var)) {
      coef_vec(173) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, size2var)) {
      coef_vec(174) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, size2var)) {
      coef_vec(175) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, size2var)) {
      coef_vec(176) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova2var, repst1var)) {
      coef_vec(177) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb2var, repst1var)) {
      coef_vec(178) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc2var, repst1var)) {
      coef_vec(179) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcova1var, repst2var)) {
      coef_vec(180) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovb1var, repst2var)) {
      coef_vec(181) = fixed_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_vars(i)), indcovc1var, repst2var)) {
      coef_vec(182) = fixed_slopes(i);
    }
  }
  
  for (int i = 0; i < no_fixed_zi_vars; i++) {
    for (int j = 0; j < no_years; j++) {
      if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(mainyears_text(j)))) {
        zeroyear_coefs(j) = fixed_zi_slopes(i);
      }
    }
    for (int j = 0; j < no_patches; j++) {
      if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), patchvar)) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(mainpatches(j)))) {
          zeropatch_coefs(j) = fixed_zi_slopes(i);
        }
      }
    }
    
    for (int j = 0; j < no_groups; j++) {
      if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), group2var)) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(maingroups_text(j)))) {
          zerogroup2_coefs(j) = fixed_zi_slopes(i);
        }
      }
    }
    for (int j = 0; j < no_groups; j++) {
      if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), group1var)) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(maingroups_text(j)))) {
          zerogroup1_coefs(j) = fixed_zi_slopes(i);
        }
      }
    }
    
    if (no_indcova_names > 0) {
      for (int j = 0; j < no_indcova_names; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcova2var)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(indcova_names(j)))) {
            zeroindcova2s(j) = fixed_zi_slopes(i);
          }
        }
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcova1var)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(indcova_names(j)))) {
            zeroindcova1s(j) = fixed_zi_slopes(i);
          }
        }
      }
    }
    
    if (no_indcovb_names > 0) {
      for (int j = 0; j < no_indcovb_names; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovb2var)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(indcovb_names(j)))) {
            zeroindcovb2s(j) = fixed_zi_slopes(i);
          }
        }
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovb1var)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(indcovb_names(j)))) {
            zeroindcovb1s(j) = fixed_zi_slopes(i);
          }
        }
      }
    }
    
    if (no_indcovc_names > 0) {
      for (int j = 0; j < no_indcovc_names; j++) {
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovc2var)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(indcovc_names(j)))) {
            zeroindcovc2s(j) = fixed_zi_slopes(i);
          }
        }
        if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), indcovc1var)) {
          if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), as<std::string>(indcovc_names(j)))) {
            zeroindcovc1s(j) = fixed_zi_slopes(i);
          }
        }
      }
    }
    
    if (stringcompare_simple(as<std::string>(fixed_zi_vars(i)), "ntercep")) {
      coef_vec(46) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), repst1var)) {
      coef_vec(47) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), repst2var)) {
      coef_vec(48) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), size1var)) {
      coef_vec(49) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), size2var)) {
      coef_vec(50) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst1var, repst2var)) {
      coef_vec(51) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, size2var)) {
      coef_vec(52) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, repst1var)) {
      coef_vec(53) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, repst2var)) {
      coef_vec(54) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, repst1var)) {
      coef_vec(55) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, repst2var)) {
      coef_vec(56) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), agevar)) {
      coef_vec(57) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, agevar)) {
      coef_vec(58) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, agevar)) {
      coef_vec(59) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst1var, agevar)) {
      coef_vec(60) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst2var, agevar)) {
      coef_vec(61) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcova2var)) {
      coef_vec(62) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovb2var)) {
      coef_vec(63) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovc2var)) {
      coef_vec(64) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcova1var)) {
      coef_vec(65) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovb1var)) {
      coef_vec(66) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), indcovc1var)) {
      coef_vec(67) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, size2var)) {
      coef_vec(68) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, size2var)) {
      coef_vec(69) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, size2var)) {
      coef_vec(70) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, repst2var)) {
      coef_vec(71) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, repst2var)) {
      coef_vec(72) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, repst2var)) {
      coef_vec(73) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, size1var)) {
      coef_vec(74) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, size1var)) {
      coef_vec(75) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, size1var)) {
      coef_vec(76) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, repst1var)) {
      coef_vec(77) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, repst1var)) {
      coef_vec(78) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, repst1var)) {
      coef_vec(79) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovb2var)) {
      coef_vec(80) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovc2var)) {
      coef_vec(81) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, indcovc2var)) {
      coef_vec(82) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovb1var)) {
      coef_vec(83) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovc1var)) {
      coef_vec(84) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, indcovc1var)) {
      coef_vec(85) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovb1var)) {
      coef_vec(86) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovb2var)) {
      coef_vec(87) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, indcovc1var)) {
      coef_vec(88) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, indcovc2var)) {
      coef_vec(89) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, indcovc1var)) {
      coef_vec(90) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, indcovc2var)) {
      coef_vec(91) = fixed_zi_slopes(i);
    }
    
    
    // New coefficients
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizeb2var)) {
      coef_vec(200) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizeb1var)) {
      coef_vec(201) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizec2var)) {
      coef_vec(202) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), sizec1var)) {
      coef_vec(203) = fixed_zi_slopes(i);
    }
    if (stringcompare_hard(as<std::string>(fixed_zi_vars(i)), densityvar)) {
      coef_vec(204) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, sizeb2var)) {
      coef_vec(205) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, sizec2var)) {
      coef_vec(206) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizeb1var)) {
      coef_vec(207) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizec1var)) {
      coef_vec(208) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, sizec1var)) {
      coef_vec(209) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizeb2var)) {
      coef_vec(210) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizec2var)) {
      coef_vec(211) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, sizec2var)) {
      coef_vec(212) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizeb2var)) {
      coef_vec(213) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, sizec2var)) {
      coef_vec(214) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, sizec2var)) {
      coef_vec(215) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizeb1var)) {
      coef_vec(216) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, sizec1var)) {
      coef_vec(217) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, sizec1var)) {
      coef_vec(218) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size2var, densityvar)) {
      coef_vec(219) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, densityvar)) {
      coef_vec(220) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, densityvar)) {
      coef_vec(221) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), size1var, densityvar)) {
      coef_vec(222) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, densityvar)) {
      coef_vec(223) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, densityvar)) {
      coef_vec(224) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst2var, densityvar)) {
      coef_vec(225) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), repst1var, densityvar)) {
      coef_vec(226) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, repst2var)) {
      coef_vec(227) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, repst2var)) {
      coef_vec(228) = fixed_zi_slopes(i);
    }
    // 229 = 0
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, repst1var)) {
      coef_vec(230) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, repst1var)) {
      coef_vec(231) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, repst2var)) {
      coef_vec(232) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, repst1var)) {
      coef_vec(233) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, repst1var)) {
      coef_vec(234) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, repst2var)) {
      coef_vec(235) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb2var, agevar)) {
      coef_vec(236) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec2var, agevar)) {
      coef_vec(237) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), densityvar, agevar)) {
      coef_vec(238) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizeb1var, agevar)) {
      coef_vec(239) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), sizec1var, agevar)) {
      coef_vec(240) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizeb2var)) {
      coef_vec(241) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizec2var)) {
      coef_vec(242) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, densityvar)) {
      coef_vec(243) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizeb1var)) {
      coef_vec(244) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizec1var)) {
      coef_vec(245) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizeb2var)) {
      coef_vec(246) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, sizec2var)) {
      coef_vec(247) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizeb1var)) {
      coef_vec(248) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, sizec1var)) {
      coef_vec(249) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, densityvar)) {
      coef_vec(250) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizeb2var)) {
      coef_vec(251) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizec2var)) {
      coef_vec(252) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, densityvar)) {
      coef_vec(253) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizeb1var)) {
      coef_vec(254) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizec1var)) {
      coef_vec(255) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizeb2var)) {
      coef_vec(256) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, sizec2var)) {
      coef_vec(257) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizeb1var)) {
      coef_vec(258) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, sizec1var)) {
      coef_vec(259) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, densityvar)) {
      coef_vec(260) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizeb2var)) {
      coef_vec(261) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizec2var)) {
      coef_vec(262) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, densityvar)) {
      coef_vec(263) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizeb1var)) {
      coef_vec(264) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizec1var)) {
      coef_vec(265) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizeb2var)) {
      coef_vec(266) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, sizec2var)) {
      coef_vec(267) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizeb1var)) {
      coef_vec(268) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, sizec1var)) {
      coef_vec(269) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, densityvar)) {
      coef_vec(270) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, size1var)) {
      coef_vec(271) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, size1var)) {
      coef_vec(272) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, size1var)) {
      coef_vec(273) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, size2var)) {
      coef_vec(274) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, size2var)) {
      coef_vec(275) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, size2var)) {
      coef_vec(276) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova2var, repst1var)) {
      coef_vec(277) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb2var, repst1var)) {
      coef_vec(278) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc2var, repst1var)) {
      coef_vec(279) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcova1var, repst2var)) {
      coef_vec(280) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovb1var, repst2var)) {
      coef_vec(281) = fixed_zi_slopes(i);
    }
    if (stringcompare_x(as<std::string>(fixed_zi_vars(i)), indcovc1var, repst2var)) {
      coef_vec(282) = fixed_zi_slopes(i);
    }
  }
  
  // Random slopes
  Nullable<List> random_vars_ = core_components["random_vars"];
  Nullable<List> random_zi_vars_ = core_components["random_zi_vars"];
  Nullable<List> random_slopes_ = core_components["random_slopes"];
  Nullable<List> random_zi_slopes_ = core_components["random_zi_slopes"];
  
  if (random_vars_.isNotNull()) {
    random_vars = random_vars_;
    random_slopes = random_slopes_;
    
    CharacterVector random_names = random_vars.attr("names");
    int no_random_vars = random_names.length();
    
    for (int i = 0; i < no_random_vars; i++) {
      if (stringcompare_hard(as<std::string>(random_names(i)), year2var)) {
        CharacterVector ran_year_names = random_vars[i];
        NumericVector ran_year_slopes = random_slopes[i];
        int no_ran_year_slopes = ran_year_names.length();
        
        for (int j = 0; j < no_ran_year_slopes; j++) {
          for (int k = 0; k < no_years; k++) {
            if (stringcompare_hard(as<std::string>(ran_year_names(j)), as<std::string>(mainyears_text(k)))) {
              year_coefs(k) = ran_year_slopes(j);
            }
          }
        }
      }
      
      if (stringcompare_hard(as<std::string>(random_names(i)), patchvar)) {
        CharacterVector ran_patch_names = random_vars[i];
        NumericVector ran_patch_slopes = random_slopes[i];
        int no_ran_patch_slopes = ran_patch_names.length();
        
        for (int j = 0; j < no_ran_patch_slopes; j++) {
          for (int k = 0; k < no_patches; k++) {
            if (stringcompare_hard(as<std::string>(ran_patch_names(j)), as<std::string>(mainpatches(k)))) {
              patch_coefs(k) = ran_patch_slopes(j);
            }
          }
        }
      }
      
      if (stringcompare_hard(as<std::string>(random_names(i)), indcova2var)) {
        CharacterVector ran_inda2_names = random_vars[i];
        NumericVector ran_inda2_slopes = random_slopes[i];
        int no_ran_inda2_slopes = ran_inda2_names.length();
      
        for (int j = 0; j < no_ran_inda2_slopes; j++) {
          for (int k = 0; k < no_indcova_names; k++) {
            if (stringcompare_hard(as<std::string>(ran_inda2_names(j)), as<std::string>(indcova_names(k)))) {
              indcova2s(k) = ran_inda2_slopes(j);
            }
          }
        }
      }
      if (stringcompare_hard(as<std::string>(random_names(i)), indcova1var)) {
        CharacterVector ran_inda1_names = random_vars[i];
        NumericVector ran_inda1_slopes = random_slopes[i];
        int no_ran_inda1_slopes = ran_inda1_names.length();
      
        for (int j = 0; j < no_ran_inda1_slopes; j++) {
          for (int k = 0; k < no_indcova_names; k++) {
            if (stringcompare_hard(as<std::string>(ran_inda1_names(j)), as<std::string>(indcova_names(k)))) {
              indcova1s(k) = ran_inda1_slopes(j);
            }
          }
        }
      }
      if (stringcompare_hard(as<std::string>(random_names(i)), indcovb2var)) {
        CharacterVector ran_indb2_names = random_vars[i];
        NumericVector ran_indb2_slopes = random_slopes[i];
        int no_ran_indb2_slopes = ran_indb2_names.length();
      
        for (int j = 0; j < no_ran_indb2_slopes; j++) {
          for (int k = 0; k < no_indcovb_names; k++) {
            if (stringcompare_hard(as<std::string>(ran_indb2_names(j)), as<std::string>(indcovb_names(k)))) {
              indcovb2s(k) = ran_indb2_slopes(j);
            }
          }
        }
      }
      if (stringcompare_hard(as<std::string>(random_names(i)), indcovb1var)) {
        CharacterVector ran_indb1_names = random_vars[i];
        NumericVector ran_indb1_slopes = random_slopes[i];
        int no_ran_indb1_slopes = ran_indb1_names.length();
      
        for (int j = 0; j < no_ran_indb1_slopes; j++) {
          for (int k = 0; k < no_indcovb_names; k++) {
            if (stringcompare_hard(as<std::string>(ran_indb1_names(j)), as<std::string>(indcovb_names(k)))) {
              indcovb1s(k) = ran_indb1_slopes(j);
            }
          }
        }
      }
      if (stringcompare_hard(as<std::string>(random_names(i)), indcovc2var)) {
        CharacterVector ran_indc2_names = random_vars[i];
        NumericVector ran_indc2_slopes = random_slopes[i];
        int no_ran_indc2_slopes = ran_indc2_names.length();
      
        for (int j = 0; j < no_ran_indc2_slopes; j++) {
          for (int k = 0; k < no_indcovc_names; k++) {
            if (stringcompare_hard(as<std::string>(ran_indc2_names(j)), as<std::string>(indcovc_names(k)))) {
              indcovc2s(k) = ran_indc2_slopes(j);
            }
          }
        }
      }
      if (stringcompare_hard(as<std::string>(random_names(i)), indcovc1var)) {
        CharacterVector ran_indc1_names = random_vars[i];
        NumericVector ran_indc1_slopes = random_slopes[i];
        int no_ran_indc1_slopes = ran_indc1_names.length();
      
        for (int j = 0; j < no_ran_indc1_slopes; j++) {
          for (int k = 0; k < no_indcovc_names; k++) {
            if (stringcompare_hard(as<std::string>(ran_indc1_names(j)), as<std::string>(indcovc_names(k)))) {
              indcovc1s(k) = ran_indc1_slopes(j);
            }
          }
        }
      }
      
    }
  }
  
  if (random_zi_vars_.isNotNull() && random_zi_slopes_.isNotNull()) {
    
    random_zi_vars = random_zi_vars_;
    random_zi_slopes = random_zi_slopes_;
    
    if (random_zi_slopes.length() > 0) {
      CharacterVector random_zi_names = random_zi_vars.attr("names");
      int no_random_zi_vars = random_zi_names.length();
      
      for (int i = 0; i < no_random_zi_vars; i++) {
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), year2var)) {
          CharacterVector ran_zi_year_names = random_zi_vars[i];
          NumericVector ran_zi_year_slopes = random_zi_slopes[i];
          int no_ran_zi_year_slopes = ran_zi_year_names.length();
          
          for (int j = 0; j < no_ran_zi_year_slopes; j++) {
            for (int k = 0; k < no_years; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_year_names(j)), 
                  as<std::string>(mainyears_text(k)))) {
                zeroyear_coefs(k) = ran_zi_year_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), patchvar)) {
          CharacterVector ran_zi_patch_names = random_zi_vars[i];
          NumericVector ran_zi_patch_slopes = random_zi_slopes[i];
          int no_ran_zi_patch_slopes = ran_zi_patch_names.length();
          
          for (int j = 0; j < no_ran_zi_patch_slopes; j++) {
            for (int k = 0; k < no_patches; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_patch_names(j)), as<std::string>(mainpatches(k)))) {
                zeropatch_coefs(k) = ran_zi_patch_slopes(j);
              }
            }
          }
        }
        
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcova2var)) {
          CharacterVector ran_zi_inda2_names = random_zi_vars[i];
          NumericVector ran_zi_inda2_slopes = random_zi_slopes[i];
          int no_ran_zi_inda2_slopes = ran_zi_inda2_names.length();
        
          for (int j = 0; j < no_ran_zi_inda2_slopes; j++) {
            for (int k = 0; k < no_indcova_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_inda2_names(j)),
                as<std::string>(indcova_names(k)))) {
                zeroindcova2s(k) = ran_zi_inda2_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcova1var)) {
          CharacterVector ran_zi_inda1_names = random_zi_vars[i];
          NumericVector ran_zi_inda1_slopes = random_zi_slopes[i];
          int no_ran_zi_inda1_slopes = ran_zi_inda1_names.length();
        
          for (int j = 0; j < no_ran_zi_inda1_slopes; j++) {
            for (int k = 0; k < no_indcova_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_inda1_names(j)),
                as<std::string>(indcova_names(k)))) {
                zeroindcova1s(k) = ran_zi_inda1_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovb2var)) {
          CharacterVector ran_zi_indb2_names = random_zi_vars[i];
          NumericVector ran_zi_indb2_slopes = random_zi_slopes[i];
          int no_ran_zi_indb2_slopes = ran_zi_indb2_names.length();
        
          for (int j = 0; j < no_ran_zi_indb2_slopes; j++) {
            for (int k = 0; k < no_indcovb_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_indb2_names(j)),
                as<std::string>(indcovb_names(k)))) {
                zeroindcovb2s(k) = ran_zi_indb2_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovb1var)) {
          CharacterVector ran_zi_indb1_names = random_zi_vars[i];
          NumericVector ran_zi_indb1_slopes = random_zi_slopes[i];
          int no_ran_zi_indb1_slopes = ran_zi_indb1_names.length();
        
          for (int j = 0; j < no_ran_zi_indb1_slopes; j++) {
            for (int k = 0; k < no_indcovb_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_indb1_names(j)),
                as<std::string>(indcovb_names(k)))) {
                zeroindcovb1s(k) = ran_zi_indb1_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovc2var)) {
          CharacterVector ran_zi_indc2_names = random_zi_vars[i];
          NumericVector ran_zi_indc2_slopes = random_zi_slopes[i];
          int no_ran_zi_indc2_slopes = ran_zi_indc2_names.length();
        
          for (int j = 0; j < no_ran_zi_indc2_slopes; j++) {
            for (int k = 0; k < no_indcovc_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_indc2_names(j)),
                as<std::string>(indcovc_names(k)))) {
                zeroindcovc2s(k) = ran_zi_indc2_slopes(j);
              }
            }
          }
        }
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcovc1var)) {
          CharacterVector ran_zi_indc1_names = random_zi_vars[i];
          NumericVector ran_zi_indc1_slopes = random_zi_slopes[i];
          int no_ran_zi_indc1_slopes = ran_zi_indc1_names.length();
        
          for (int j = 0; j < no_ran_zi_indc1_slopes; j++) {
            for (int k = 0; k < no_indcovc_names; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_indc1_names(j)),
                as<std::string>(indcovc_names(k)))) {
                zeroindcovc1s(k) = ran_zi_indc1_slopes(j);
              }
            }
          }
        }
        
      }
    }
  }
  
  CharacterVector noneslot {"none"};
  if (no_indcova_names == 0) {
    indcova2s = {0};
    indcova1s = {0};
    indcova2s.attr("names") = noneslot;
    indcova1s.attr("names") = noneslot;
    
    zeroindcova2s = {0};
    zeroindcova1s = {0};
    zeroindcova2s.attr("names") = noneslot;
    zeroindcova1s.attr("names") = noneslot;
  }
  if (no_indcovb_names == 0) {
    indcovb2s = {0};
    indcovb1s = {0};
    indcovb2s.attr("names") = noneslot;
    indcovb1s.attr("names") = noneslot;
    
    zeroindcovb2s = {0};
    zeroindcovb1s = {0};
    zeroindcovb2s.attr("names") = noneslot;
    zeroindcovb1s.attr("names") = noneslot;
  }
  if (no_indcovc_names == 0) {
    indcovc2s = {0};
    indcovc1s = {0};
    indcovc2s.attr("names") = noneslot;
    indcovc1s.attr("names") = noneslot;
    
    zeroindcovc2s = {0};
    zeroindcovc1s = {0};
    zeroindcovc2s.attr("names") = noneslot;
    zeroindcovc1s.attr("names") = noneslot;
  }
  
  List output(28);
  
  output(0) = coef_vec;
  output(1) = year_coefs;
  output(2) = zeroyear_coefs;
  output(3) = patch_coefs;
  output(4) = zeropatch_coefs;
  output(5) = group2_coefs;
  output(6) = group1_coefs;
  output(7) = zerogroup2_coefs;
  output(8) = zerogroup1_coefs;
  output(9) = indcova2s;
  output(10) = indcova1s;
  output(11) = indcovb2s;
  output(12) = indcovb1s;
  output(13) = indcovc2s;
  output(14) = indcovc1s;
  output(15) = zeroindcova2s;
  output(16) = zeroindcova1s;
  output(17) = zeroindcovb2s;
  output(18) = zeroindcovb1s;
  output(19) = zeroindcovc2s;
  output(20) = zeroindcovc1s;
  output(21) = as<CharacterVector>(core_components["class"]);
  output(22) = as<CharacterVector>(core_components["family"]);
  output(23) = as<IntegerVector>(core_components["dist"]);
  output(24) = as<LogicalVector>(core_components["zero_inflated"]);
  output(25) = as<LogicalVector>(core_components["zero_truncated"]);
  output(26) = as<NumericVector>(core_components["sigma"]);
  output(27) = as<NumericVector>(core_components["theta"]);

  CharacterVector output_names = {"coefficients", "years", "zeroyear", "patches",
    "zeropatch", "groups2", "groups1", "zerogroups2", "zerogroups1", "indcova2s",
    "indcova1s", "indcovb2s", "indcovb1s", "indcovc2s", "indcovc1s",
    "zeroindcova2s", "zeroindcova1s", "zeroindcovb2s", "zeroindcovb1s",
    "zeroindcovc2s", "zeroindcovc1s", "class", "family", "dist", "zero_inflated",
    "zero_truncated", "sigma", "theta"};
  output.attr("names") = output_names;
  
  return output;
}

//' Key Function Passing Models and Other Parameters to Matrix Estimators
//' 
//' Function \code{raymccooney()} takes the various vital rate models and other
//' parameters and coordinates them as input into the function-based matrix
//' estimation functions.
//' 
//' @param listofyears A data frame where the rows designate the exact order of
//' years and patches to produce matrices for.
//' @param modelsuite An object of class \code{lefkoMod}, or a similarly
//' structured list object. All 14 vital rate models and the \code{paramnames}
//' data frame are required.
//' @param mainyears A numeric vector of all times at time \emph{t}.
//' @param mainpatches A string vector of patch names.
//' @param maingroups A string vector of stage group names.
//' @param mainindcova Typically a string vector of individual covariate
//' category names.
//' @param mainindcovb Typically a string vector of individual covariate
//' category names.
//' @param mainindcovc Typically a string vector of individual covariate
//' category names.
//' @param StageFrame The stageframe object identifying the life history model
//' being operationalized.
//' @param OverWrite The overwrite table used in analysis, as modified by 
//' \code{.overwrite_reassess}. Must be processed via \code{.overwrite_reassess}
//' rather than being a raw overwrite or supplement table.
//' @param repmatrix The reproductive matrix used in analysis.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' juvenile reproductive status, and juvenile maturity status.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param firstage The first age to be used in the analysis. Should typically
//' be \code{0} for pre-breeding and \code{1} for post-breeding life history
//' models. If not building age-by-stage MPMs, then should be set to \code{0}.
//' @param finalage The final age to be used in analysis. If not building
//' age-by-stage MPMs, then should be set to \code{0}.
//' @param format Indicates whether historical matrices should be in (1) Ehrlen
//' or (2) deVries format.
//' @param style The style of analysis, where 0 is historical, 1 is ahistorical,
//' and 2 is age-by-stage.
//' @param cont Denotes whether age-by-stage matrix continues past the final
//' age.
//' @param filter An integer denoting whether to filter the DataFrame to
//' eliminate unusable rows, and if so, how to do so. Possible values: \code{0}:
//' no filtering, \code{1}: filter out rows with \code{index321 == -1}, and
//' \code{2}: filter out rows with \code{aliveandequal == -1}.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param ipm_method A string indicating which method should be used to
//' estimate size transitions in cases with continuous distributions. Options
//' include \code{"midpoint"}, which uses the midpoint method, and \code{"cdf"},
//' which uses the cumulative density function.
//' @param err_check If \code{TRUE}, then also output objects \code{prob_out}
//' and \code{allstages} for error checking purposes.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' 
//' @return A list with with up to 5 elements. In order: \code{A}: a list of A
//' matrices, or a list of \code{NULL} values if \code{simplicity = TRUE};
//' \code{U}: a list of U matrices, in the same order as \code{A}; \code{F}:
//' a list of F matrices, in the same order as \code{A}; \code{prob_out}: a list
//' of error-checking conditional probability matrices, or a list of \code{NULL}
//' values if \code{err_check = FALSE}; and \code{allstages}: a data frame
//' showing the used values of all variables used in transition calculations.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.raymccooney)]]
List raymccooney(DataFrame listofyears, List modelsuite, NumericVector mainyears,
  CharacterVector mainpatches, RObject maingroups, RObject mainindcova,
  RObject mainindcovb, RObject mainindcovc, DataFrame StageFrame,
  DataFrame OverWrite, arma::mat repmatrix, NumericVector f2_inda,
  NumericVector f1_inda, NumericVector f2_indb, NumericVector f1_indb,
  NumericVector f2_indc, NumericVector f1_indc, StringVector r2_inda,
  StringVector r1_inda, StringVector r2_indb, StringVector r1_indb,
  StringVector r2_indc, StringVector r1_indc, NumericVector dev_terms,
  double dens, double fecmod, int firstage, int finalage, int format, int style,
  int cont, int filter, bool negfec, double exp_tol = 700.0,
  double theta_tol = 100000000.0, String ipm_method = "cdf",
  bool err_check = false, bool simplicity = false) {
  
  // listofyears import and settings
  IntegerVector years = listofyears["yearorder"];
  IntegerVector patches = listofyears["patchorder"];
  int loy_length = years.length();
  
  int matrixformat {0};
  if (style == 2) {
    matrixformat = 4;
  } else if (style == 1) {
    matrixformat = 3;
  } else if (style == 0) {
    if (format == 1) {
      matrixformat = 1;
    } else if (format == 2) {
      matrixformat = 2;
    } else {
      throw Rcpp::exception("Matrix format is not recognized.", false);
    }
  } else {
    throw Rcpp::exception("Matrix style is not recognized.", false);
  }
  
  Rcpp::DataFrame allstages = theoldpizzle(StageFrame, OverWrite, repmatrix,
    firstage, finalage, format, style, cont, filter);
  
  NumericVector size3 = allstages["size3"];
  NumericVector size2n = allstages["size2n"];
  NumericVector size2o = allstages["size2o"];
  NumericVector sizeb3 = allstages["sizeb3"];
  NumericVector sizeb2n = allstages["sizeb2n"];
  NumericVector sizeb2o = allstages["sizeb2o"];
  NumericVector sizec3 = allstages["sizec3"];
  NumericVector sizec2n = allstages["sizec2n"];
  NumericVector sizec2o = allstages["sizec2o"];
  
  NumericVector maxveca = {max(size3), max(size2n), max(size2o)}; // What about NAs?
  NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
  NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
  
  double maxsize = max(maxveca); // What about NAs?
  double maxsizeb = max(maxvecb);
  double maxsizec = max(maxvecc);
  
  RObject surv_model = modelsuite["survival_model"];
  RObject obs_model = modelsuite["observation_model"];
  RObject size_model = modelsuite["size_model"];
  RObject sizeb_model = modelsuite["sizeb_model"];
  RObject sizec_model = modelsuite["sizec_model"];
  RObject repst_model = modelsuite["repstatus_model"];
  RObject fec_model = modelsuite["fecundity_model"];
  RObject jsurv_model = modelsuite["juv_survival_model"];
  RObject jobs_model = modelsuite["juv_observation_model"];
  RObject jsize_model = modelsuite["juv_size_model"];
  RObject jsizeb_model = modelsuite["juv_sizeb_model"];
  RObject jsizec_model = modelsuite["juv_sizec_model"];
  RObject jrepst_model = modelsuite["juv_reproduction_model"];
  RObject jmatst_model = modelsuite["juv_maturity_model"];
  DataFrame paramnames = as<DataFrame>(modelsuite["paramnames"]);
  
  List surv_proxy = modelextract(surv_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List obs_proxy = modelextract(obs_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List size_proxy = modelextract(size_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List sizeb_proxy = modelextract(sizeb_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List sizec_proxy = modelextract(sizec_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List repst_proxy = modelextract(repst_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List fec_proxy = modelextract(fec_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jsurv_proxy = modelextract(jsurv_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jobs_proxy = modelextract(jobs_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jsize_proxy = modelextract(jsize_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jsizeb_proxy = modelextract(jsizeb_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jsizec_proxy = modelextract(jsizec_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jrepst_proxy = modelextract(jrepst_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List jmatst_proxy = modelextract(jmatst_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  
  // Now we create the matrices and order them within the correct lsit structure
  List A_mats(loy_length);
  List F_mats(loy_length);
  List U_mats(loy_length);
  List out_mats(loy_length);
  
  int yearnumber {0};
  int patchnumber {0};
  
  for (int i = 0; i < loy_length; i++) {
    
    yearnumber = years(i) - 1;
    patchnumber = patches(i) - 1;
    
    List madsexmadrigal_oneyear = jerzeibalowski(allstages, StageFrame,
      matrixformat, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
      repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
      jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda, f1_inda, f2_indb, f1_indb,
      f2_indc, f1_indc, r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc,
      dev_terms, dens, fecmod, maxsize, maxsizeb, maxsizec, firstage, finalage,
      negfec, yearnumber, patchnumber, exp_tol, theta_tol, ipm_method, err_check,
      simplicity);
    
    if (!simplicity) A_mats(i) = madsexmadrigal_oneyear["A"];
    F_mats(i) = madsexmadrigal_oneyear["F"];
    U_mats(i) = madsexmadrigal_oneyear["U"];
    if (err_check ) out_mats(i) = madsexmadrigal_oneyear["out"];
  }
  
  List output;
  
  if (simplicity && err_check) {
    output = List::create(_["U"] = U_mats, _["F"] = F_mats, _["prob_out"] = out_mats,
      _["allstages"] = allstages);
  } else if (simplicity) {
    output = List::create(_["U"] = U_mats, _["F"] = F_mats);
  } else if (err_check) {
    output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats,
      _["prob_out"] = out_mats, _["allstages"] = allstages);
  } else {
    output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats);
  }
  
  return output;
}

//' Function Passing Models and Other Parameters to Leslie Matrix Estimator
//' 
//' This function takes the various vital rate models and other parameters and
//' coordinates them as input into function \code{fleslie()}.
//' 
//' @param listofyears A data frame where the rows designate the exact order of
//' years and patches to produce matrices for.
//' @param modelsuite An object of class \code{lefkoMod}, or a similarly
//' structured list object. Survival model, fecundity model, and the
//' \code{paramnames} data frame are required.
//' @param actualages An integer vector of all actual ages to be included in the
//' matrices, in order.
//' @param mainyears A numeric vector of all times at time \emph{t}.
//' @param mainpatches A string vector of patch names.
//' @param maingroups A string vector of stage group names.
//' @param mainindcova Typically a string vector of individual covariate
//' category names.
//' @param mainindcovb Typically a string vector of individual covariate
//' category names.
//' @param mainindcovc Typically a string vector of individual covariate
//' category names.
//' @param ageframe The modified stageframe used in matrix calculations.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, and fecundity.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param finalage The final age to be used in analysis.
//' @param cont Denotes whether age-by-stage matrix continues past the final
//' age.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param simplicity If \code{TRUE}, then only outputs matrices \code{U} and
//' \code{F}, rather than also outputting matrix \code{A}. Defaults to
//' \code{FALSE}.
//' 
//' @return A list with with up to 5 elements. In order: \code{A}: a list of A
//' matrices, or a list of \code{NULL} values if \code{simplicity = TRUE};
//' \code{U}: a list of U matrices, in the same order as \code{A}; \code{F}:
//' a list of F matrices, in the same order as \code{A}; \code{prob_out}: a list
//' of error-checking conditional probability matrices, or a list of \code{NULL}
//' values if \code{err_check = FALSE}; and \code{allstages}: a data frame
//' showing the used values of all variables used in transition calculations.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.mothermccooney)]]
List mothermccooney(DataFrame listofyears, List modelsuite, IntegerVector actualages,
  NumericVector mainyears, CharacterVector mainpatches, RObject maingroups,
  RObject mainindcova, RObject mainindcovb, RObject mainindcovc, DataFrame ageframe,
  NumericVector f2_inda, NumericVector f1_inda, NumericVector f2_indb,
  NumericVector f1_indb, NumericVector f2_indc, NumericVector f1_indc,
  StringVector r2_inda, StringVector r1_inda, StringVector r2_indb,
  StringVector r1_indb, StringVector r2_indc, StringVector r1_indc,
  NumericVector dev_terms, double dens, double fecmod, int finalage, int cont,
  bool negfec, double exp_tol = 700.0, double theta_tol = 100000000.0,
  bool err_check = false, bool simplicity = false) {
  
  // listofyears import and settings
  IntegerVector years = listofyears["yearorder"];
  IntegerVector patches = listofyears["patchorder"];
  int loy_length = years.length();
  
  // Deviation terms
  double surv_dev = dev_terms(0);
  double fec_dev = dev_terms(1);
  
  RObject surv_model = modelsuite["survival_model"];
  RObject fec_model = modelsuite["fecundity_model"];
  DataFrame paramnames = as<DataFrame>(modelsuite["paramnames"]);
  
  List surv_proxy = modelextract(surv_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  List fec_proxy = modelextract(fec_model, paramnames,
    mainyears, mainpatches, maingroups,
    mainindcova, mainindcovb, mainindcovc);
  
  // Now we create the matrices and order them within the correct list structure
  List A_mats(loy_length);
  List F_mats(loy_length);
  List U_mats(loy_length);
  
  int yearnumber {0};
  int patchnumber {0};
  
  for (int i = 0; i < loy_length; i++) {
    
    yearnumber = years(i) - 1;
    patchnumber = patches(i) - 1;
    
    List madsexmadrigal_oneyear = motherbalowski(actualages, ageframe, surv_proxy,
      fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc, r2_inda,
      r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, surv_dev, fec_dev, dens, fecmod,
      finalage, negfec, yearnumber, patchnumber, exp_tol, theta_tol, simplicity);
    
    if (!simplicity) A_mats(i) = madsexmadrigal_oneyear["A"];
    F_mats(i) = madsexmadrigal_oneyear["F"];
    U_mats(i) = madsexmadrigal_oneyear["U"];
  }
  
  List output;
  
  if (simplicity) {
    output = List::create(_["U"] = U_mats, _["F"] = F_mats);
  } else {
    output = List::create(_["A"] = A_mats, _["U"] = U_mats, _["F"] = F_mats);
  }
  return output;
}

//' Project Function-based Matrix Projection Model
//' 
//' Function \code{f_projection3()} develops and projects function-based matrix
//' models. Unlike \code{\link{projection3}()}, which uses matrices provided as
//' input via already created \code{lefkoMat} objects, function
//' \code{f_projection3()} creates matrices at each time step from vital rate
//' models and parameter inputs provided. Projections may be deterministic or
//' stochastic, and may be density dependent in either case. Replicates may also
//' be produced.
//' 
//' @name f_projection3
//' 
//' @param data The historical vertical demographic data frame used to estimate
//' vital rates (class \code{hfvdata}), which is required to initialize times and
//' patches properly. Variable names should correspond to the naming conventions
//' in \code{\link{verticalize3}()} and \code{\link{historicalize3}()}.
//' @param format An integer indicating the kind of function-based MPM to create.
//' Possible choices include: \code{1}, Ehrlen-format historical MPM; \code{2},
//' deVries-format historical MPM; \code{3}, ahistorical MPM; \code{4},
//' age-by-stage MPM; and \code{5}, Leslie (age-based) MPM.
//' @param prebreeding A logical value indicating whether the life history model
//' is a pre-breeding model. Only used in Leslie and age-by-stage MPMs. Defaults
//' to \code{TRUE}.
//' @param start_age The age from which to start the matrix. Defaults to
//' \code{NA}, in which case age \code{1} is used if \code{prebreeding = TRUE},
//' and age \code{0} is used if \code{prebreeding = FALSE}.
//' @param last_age The final age to use in the matrix. Defaults to \code{NA}, in
//' which case the highest age in the dataset is used.
//' @param fecage_min The minimum age at which reproduction is possible. Defaults
//' to \code{NA}, which is interpreted to mean that fecundity should be assessed
//' starting in the minimum age observed in the dataset.
//' @param fecage_max The maximum age at which reproduction is possible. Defaults
//' to \code{NA}, which is interpreted to mean that fecundity should be assessed
//' until the final observed age.
//' @param cont A logical value designating whether to allow continued survival
//' of individuals past the final age noted in the stageframe, using the 
//' demographic characteristics of the final age. Defaults to \code{TRUE}.
//' @param stochastic A logical value denoting whether to conduct a stochastic
//' projection or a deterministic / cyclical projection.
//' @param standardize A logical value denoting whether to re-standardize the
//' population size to \code{1.0} at each occasion. Used in density-independent
//' simulations in which it is more important to know the general trend in
//' population growth than the explicit growth rate. Defaults to \code{FALSE}.
//' @param growthonly A logical value indicating whether to produce only the
//' projected population size at each occasion (\code{TRUE}), or also to produce
//' vectors showing the stage distribution at each occasion (\code{FALSE}).
//' Defaults to \code{TRUE}.
//' @param repvalue A logical value indicating whether to calculate reproductive
//' value vectors at each time step. Can only be set to \code{TRUE} if 
//' \code{growthonly = FALSE}. Setting to \code{TRUE} may dramatically increase
//' the duration of calculations. Defaults to \code{FALSE}.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion to the nearest
//' integer. Defaults to \code{FALSE}.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0,
//' and forces fecundity to be non-negative; and \code{2} forces all column rows
//' in the survival-transition matrices to total no more than 1.0, in addition
//' to the actions outlined for option \code{1}.
//' @param ipm_method A string indicating what method to use to estimate size
//' transition probabilities, if size is treated as continuous. Options include:
//' \code{"midpoint"}, which utilizes the midpoint method; and \code{"CDF"},
//' which uses the cumulative distribution function. Defaults to \code{"CDF"}.
//' @param nreps The number of replicate projections.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' \code{10000}.
//' @param repmod A scalar multiplier of fecundity. Defaults to \code{1}.
//' @param exp_tol A numeric value used to indicate a maximum value to set
//' exponents to in the core kernel to prevent numerical overflow. Defaults to
//' \code{700}.
//' @param theta_tol A numeric value used to indicate a maximum value to theta as
//' used in the negative binomial probability density kernel. Defaults to
//' \code{100000000}, but can be reset to other values during error checking.
//' @param random_inda A logical value denoting whether to treat individual
//' covariate \code{a} as a random, categorical variable. Otherwise is treated as
//' a fixed, numeric variable. Defaults to \code{FALSE}.
//' @param random_indb A logical value denoting whether to treat individual
//' covariate \code{b} as a random, categorical variable. Otherwise is treated as
//' a fixed, numeric variable. Defaults to \code{FALSE}.
//' @param random_indc A logical value denoting whether to treat individual
//' covariate \code{c} as a random, categorical variable. Otherwise is treated as
//' a fixed, numeric variable. Defaults to \code{FALSE}.
//' @param err_check A logical value indicating whether to append extra output
//' for debugging purposes.
//' @param stageframe An object of class \code{stageframe}. These objects are
//' generated by function \code{\link{sf_create}()}, and include information on
//' the size, observation status, propagule status, reproduction status,
//' immaturity status, maturity status, stage group, size bin widths, and other
//' key characteristics of each ahistorical stage. Required for all MPM formats
//' except Leslie MPMs.
//' @param supplement An optional data frame of class \code{lefkoSD} that
//' provides supplemental data that should be incorporated into the MPM. Three
//' kinds of data may be integrated this way: transitions to be estimated via the
//' use of proxy transitions, transition overwrites from the literature or
//' supplemental studies, and transition multipliers for survival and fecundity.
//' This data frame should be produced using the \code{\link{supplemental}()}
//' function. Can be used in place of or in addition to an overwrite table (see 
//' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
//' below).
//' @param repmatrix An optional reproduction matrix. This matrix is composed
//' mostly of \code{0}s, with non-zero entries acting as element identifiers and
//' multipliers for fecundity (with \code{1} equaling full fecundity). If left
//' blank, and no \code{supplement} is provided, then \code{flefko3()} will
//' assume that all stages marked as reproductive produce offspring at 1x that of
//' estimated fecundity, and that offspring production will yield the first stage
//' noted as propagule or immature. May be the dimensions of either a historical
//' or an ahistorical matrix. If the latter, then all stages will be used in
//' occasion \emph{t}-1 for each suggested ahistorical transition.
//' @param overwrite An optional data frame developed with the
//' \code{\link{overwrite}()} function describing transitions to be overwritten
//' either with given values or with other estimated transitions. Note that this
//' function supplements overwrite data provided in \code{supplement}.
//' @param modelsuite A \code{lefkoMod} object, at minimum with all required
//' best-fit vital rate models and a \code{paramnames} data frame, and following
//' the naming conventions used in this package. If given, then
//' \code{surv_model}, \code{obs_model}, \code{size_model}, \code{sizeb_model},
//' \code{sizec_model}, \code{repst_model}, \code{fec_model}, \code{jsurv_model},
//' \code{jobs_model}, \code{jsize_model}, \code{jsizeb_model},
//' \code{jsizec_model}, \code{jrepst_model}, \code{jmatst_model},
//' \code{paramnames}, \code{yearcol}, and \code{patchcol} are not required.
//' Although this is optional input, it is recommended, and without it separate
//' vital rate model inputs (named \code{XX_model}) are required.
//' @param paramnames A data frame with three columns, the first describing all
//' terms used in linear modeling, the second (must be called \code{mainparams})
//' giving the general model terms that will be used in matrix creation, and the
//' third showing the equivalent terms used in modeling (must be named
//' \code{modelparams}). Function \code{\link{create_pm}()} can be used to
//' create a skeleton \code{paramnames} object, which can then be edited. Only
//' required if \code{modelsuite} is not supplied.
//' @param year Either a single integer value corresponding to the year to
//' project, or a vector of \code{times} elements with the year to use at each
//' time step. Defaults to \code{NA}, in which the first year in the set of years
//' in the dataset is projected. If a vector shorter than \code{times} is
//' supplied, then this vector will be cycled.
//' @param patch A value of \code{NA}, a single string value corresponding to the
//' patch to project, or a vector of \code{times} elements with the patch to use
//' at each time step. If a vector shorter than \code{times} is supplied, then
//' this vector will be cycled. Note that this function currently does not
//' handle multiple projections for different patches in the same run.
//' @param sp_density Either a single numeric value of spatial density to use in
//' vital rate models in all time steps, or a vector of \code{times} elements of
//' such numeric values. Defaults to \code{NA}.
//' @param ind_terms An optional data frame with 3 columns and \code{times} rows
//' giving the values of individual covariates a, b, and c, respectively, for
//' each projected time. Unused terms must be set to \code{0} (use of \code{NA}
//' will produce errors.)
//' @param dev_terms An optional data frame with 14 columns and \code{times}
//' rows showing the values of the deviation terms to be added to each linear
//' vital rate. The column order should be: 1: survival, 2: observation, 3:
//' primary size, 4: secondary size, 5: tertiary size, 6: reproduction, 7:
//' fecundity, 8: juvenile survival, 9: juvenile observation, 10: juvenile
//' primary size, 11: juvenile secondary size, 12: juvenile tertiary size, 13:
//' juvenile reproduction, and 14: juvenile maturity transition.  Unused terms
//' must be set to \code{0} (use of \code{NA} will produce errors.)
//' @param surv_model A linear model predicting survival probability. This can 
//' be a model of class \code{glm} or \code{glmer}, and requires a predicted
//' binomial variable under a logit link. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param obs_model A linear model predicting sprouting or observation
//' probability. This can be a model of class \code{glm} or \code{glmer}, and
//' requires a predicted binomial variable under a logit link. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param size_model A linear model predicting primary size. This can be a model
//' of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
//' \code{vglm}, \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param sizeb_model A linear model predicting secondary size. This can be a
//' model of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
//' \code{vglm}, \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param sizec_model A linear model predicting tertiary size. This can be a
//' model of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
//' \code{vglm}, \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param repst_model A linear model predicting reproduction probability. This 
//' can be a model of class \code{glm} or \code{glmer}, and requires a predicted
//' binomial variable under a logit link. Ignored if \code{modelsuite} is
//' provided. This model must have been developed in a modeling exercise testing
//' the impacts of occasions \emph{t} and \emph{t}-1.
//' @param fec_model A linear model predicting fecundity. This can be a model of
//' class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl}, \code{vglm},
//' \code{lm}, or \code{lmer}. Ignored if \code{modelsuite} is provided. This
//' model must have been developed in a modeling exercise testing the impacts of
//' occasions \emph{t} and \emph{t}-1.
//' @param jsurv_model A linear model predicting juvenile survival probability.
//' This can be a model of class \code{glm} or \code{glmer}, and requires a
//' predicted binomial variable under a logit link. Ignored if \code{modelsuite}
//' is provided. This model must have been developed in a modeling exercise
//' testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jobs_model A linear model predicting juvenile sprouting or observation
//' probability. This can be a model of class \code{glm} or \code{glmer}, and
//' requires a predicted binomial variable under a logit link. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jsize_model A linear model predicting juvenile primary size. This
//' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
//' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jsizeb_model A linear model predicting juvenile secondary size. This
//' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
//' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jsizec_model A linear model predicting juvenile tertiary size. This
//' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
//' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. Ignored if
//' \code{modelsuite} is provided. This model must have been developed in a
//' modeling exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
//' @param jrepst_model A linear model predicting reproduction probability of a 
//' mature individual that was immature in time \emph{t}. This can be a model
//' of class \code{glm} or \code{glmer}, and requires a predicted binomial
//' variable under a logit link. Ignored if \code{modelsuite} is provided. This
//' model must have been developed in a modeling exercise testing the impacts of
//' occasions \emph{t} and \emph{t}-1.
//' @param jmatst_model A linear model predicting maturity probability of an 
//' individual that was immature in time \emph{t}. This can be a model of class
//' \code{glm} or \code{glmer}, and requires a predicted binomial variable under
//' a logit link. Ignored if \code{modelsuite} is provided. This model must have
//' been developed in a modeling exercise testing the impacts of occasions
//' \emph{t} and \emph{t}-1.
//' @param start_vec An optional numeric vector denoting the starting stage
//' distribution for the projection. Defaults to a single individual of each
//' stage.
//' @param start_frame An optional data frame characterizing stages, age-stages,
//' or stage-pairs that should be set to non-zero values in the starting vector,
//' and what those values should be. Can only be used with \code{lefkoMat}
//' objects.
//' @param tweights An optional numeric vector denoting the probabilistic
//' weightings of year terms. Defaults to equal weighting among occasions.
//' @param density An optional data frame describing the matrix elements that
//' will be subject to density dependence, and the exact kind of density
//' dependence that they will be subject to. The data frame used should be an
//' object of class \code{lefkoDens}, which is the output from function
//' \code{\link{density_input}()}.
//' 
//' @return A list of class \code{lefkoProj}, which always includes the first
//' three elements of the following, and also includes the remaining elements
//' below when a \code{lefkoMat} object is used as input:
//' \item{projection}{A list of lists of matrices showing the total number of
//' individuals per stage per occasion. The first list corresponds to each
//' pop-patch followed by each population. The inner list corresponds to
//' replicates within each pop-patch or population.}
//' \item{stage_dist}{A list of lists of the actual stage distribution in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{rep_value}{A list of lists of the actual reproductive value in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{pop_size}{A list of matrices showing the total population size in each
//' occasion per replicate (row within data frame) per pop-patch or population
//' (list element).}
//' \item{labels}{A data frame showing the order of populations and patches in
//' item \code{projection}.}
//' \item{ahstages}{The original stageframe used in the study.}
//' \item{hstages}{A data frame showing the order of historical stage pairs.}
//' \item{agestages}{A data frame showing the order of age-stage pairs.}
//' \item{labels}{A short data frame indicating the population (always \code{1},
//' and patch (either the numeric index of the single chosen patch, or \code{1}
//' in all other cases).)}
//' \item{control}{A short vector indicating the number of replicates and the
//' number of occasions projected per replicate.}
//' \item{density}{The data frame input under the density option. Only provided
//' if input by the user.}
//' 
//' @section Notes:
//' Population projection can be a very time-consuming activity, and it is most
//' time-consuming when matrices need to be created at each time step. We have
//' created this function to be as quick as possible, but some options will slow
//' the analysis down. First, the \code{err_check} option should always be set
//' to \code{FALSE}, as the added created output will not only slow the analysis
//' down but also potentially crash the memory, if matrices are large enough.
//' Second, the \code{repvalue} option should be set to \code{FALSE} unless
//' reproductive values are genuinely needed, since this step requires
//' concurrent backward projection and so in some cases may double total run
//' time. Finally, if the only needed data is the total population size and
//' actual age/stage structure at each time step, then setting \code{growthonly
//' = TRUE} will yield the quickest possible run time.
//' 
//' Projections with large matrices may take a long time to run. To assess the
//' likely running time, try using a low number of iterations on a single
//' replicate first. For example, set \code{nreps = 1} and \code{times = 10} for
//' a trial run. If a full run is set and takes too long, press the STOP button
//' in RStudio to cancel the projection run.
//' 
//' Consistently positive population growth can quickly lead to population size
//' numbers larger than can be handled computationally. In that circumstance, a
//' continuously rising population size will suddenly become \code{NaN} for the
//' remainder of the projection.
//' 
//' This function does not reduce the dimensionality of matrices developed for
//' projection.
//' 
//' @seealso \code{\link{projection3}()}
//' @seealso \code{\link{flefko3}()}
//' @seealso \code{\link{flefko2}()}
//' @seealso \code{\link{aflefko2}()}
//' @seealso \code{\link{fleslie}()}
//' 
//' @examples
//' \donttest{
//' # Lathyrus projection example with historical matrices
//' data(lathyrus)
//' 
//' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
//'   9)
//' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
//'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
//'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
//' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
//'   0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
//'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
//' 
//' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
//'   propstatus = propvector)
//' 
//' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
//'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988", 
//'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
//' 
//' lathvertln$feca2 <- round(lathvertln$feca2)
//' lathvertln$feca1 <- round(lathvertln$feca1)
//' lathvertln$feca3 <- round(lathvertln$feca3)
//' 
//' lathmodelsln3 <- modelsearch(lathvertln, historical = TRUE, 
//'   approach = "mixed", suite = "main", 
//'   vitalrates = c("surv", "obs", "size", "repst", "fec"), juvestimate = "Sdl",
//'   bestfit = "AICc&k", sizedist = "gaussian", fecdist = "poisson", 
//'   indiv = "individ", patch = "patchid", year = "year2", year.as.random = TRUE,
//'   patch.as.random = TRUE, show.model.tables = TRUE, quiet = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "mat", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "Sdl", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "Sd", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, "mat", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, "Sdl", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, "Sdl", NA, NA),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframeln, historical = TRUE)
//' 
//' # While we do not use MPMs to initialize f_projections3(), we do use MPMs to
//' # initialize functions start_input() and density_input().
//' lathmat3ln <- flefko3(year = "all", patch = "all", stageframe = lathframeln, 
//'   modelsuite = lathmodelsln3, data = lathvertln, supplement = lathsupp3, 
//'   reduce = FALSE)
//' 
//' e3m_sv <- start_input(lathmat3ln, stage2 = "Sd", stage1 = "Sd", value = 1000)
//' 
//' e3d <- density_input(lathmat3ln, stage3 = c("Sd", "Sdl"),
//'   stage2 = c("rep", "rep"), stage1 = c("all", "all"), style = 1,
//'   time_delay = 1, alpha = 1, beta = 0, type = c(2, 2), type_t12 = c(1, 1))
//' 
//' trial7 <- f_projection3(format = 1, data = lathvertln,
//'   modelsuite = lathmodelsln3, stageframe = lathframeln, nreps = 2,
//'   times = 1000, stochastic = TRUE, standardize = FALSE, growthonly = TRUE,
//'   integeronly = FALSE, substoch = 0)
//' summary(trial7)
//' 
//' # Now with density dependence and a set start vector
//' trial7a <- f_projection3(format = 1, data = lathvertln,
//'   modelsuite = lathmodelsln3, stageframe = lathframeln, nreps = 2,
//'   times = 1000, stochastic = TRUE, standardize = FALSE, growthonly = TRUE,
//'   integeronly = FALSE, substoch = 0, sp_density = 0, start_frame = e3m_sv,
//'   density = e3d)
//' summary(trial7a)
//' }
//' 
//' @export f_projection3
// [[Rcpp::export]]
Rcpp::List f_projection3(DataFrame data, int format, bool prebreeding = true,
  int start_age = NA_INTEGER, int last_age = NA_INTEGER, int fecage_min = NA_INTEGER,
  int fecage_max = NA_INTEGER, bool cont = true, bool stochastic = false,
  bool standardize = false, bool growthonly = true, bool repvalue = false,
  bool integeronly = false, int substoch = 0, String ipm_method = "CDF",
  int nreps = 1, int times = 10000, double repmod = 1.0, double exp_tol = 700.0,
  double theta_tol = 100000000.0, bool random_inda = false, bool random_indb = false,
  bool random_indc = false, bool err_check = false,
  Nullable<DataFrame> stageframe = R_NilValue, Nullable<DataFrame> supplement = R_NilValue,
  Nullable<NumericMatrix> repmatrix = R_NilValue, Nullable<DataFrame> overwrite = R_NilValue,
  Nullable<List> modelsuite = R_NilValue, Nullable<DataFrame> paramnames = R_NilValue,
  Nullable<NumericVector> year = R_NilValue, Nullable<CharacterVector> patch = R_NilValue,
  Nullable<NumericVector> sp_density = R_NilValue, Nullable<RObject> ind_terms = R_NilValue,
  Nullable<RObject> dev_terms = R_NilValue, Nullable<RObject> surv_model = R_NilValue,
  Nullable<RObject> obs_model = R_NilValue, Nullable<RObject> size_model = R_NilValue,
  Nullable<RObject> sizeb_model = R_NilValue, Nullable<RObject> sizec_model = R_NilValue,
  Nullable<RObject> repst_model = R_NilValue, Nullable<RObject> fec_model = R_NilValue,
  Nullable<RObject> jsurv_model = R_NilValue, Nullable<RObject> jobs_model = R_NilValue,
  Nullable<RObject> jsize_model = R_NilValue, Nullable<RObject> jsizeb_model = R_NilValue,
  Nullable<RObject> jsizec_model = R_NilValue, Nullable<RObject> jrepst_model = R_NilValue,
  Nullable<RObject> jmatst_model = R_NilValue, Nullable<NumericVector> start_vec = R_NilValue,
  Nullable<RObject> start_frame = R_NilValue, Nullable<NumericVector> tweights = R_NilValue,
  Nullable<RObject> density = R_NilValue) {
  
  if (format < 1 || format > 5) {
    throw Rcpp::exception("Matrix format is not recognized.", false);
  }
  if (substoch <0 || substoch > 2) {
    throw Rcpp::exception("Option substoch must equal 0, 1, or 2.", false);
  }
  if (nreps < 1) {
    throw Rcpp::exception("Option nreps must equal a positive integer.", false);
  }
  if (times < 1) {
    throw Rcpp::exception("Option times must equal a positive integer.", false);
  }
  if (format < 4 && (!IntegerVector::is_na(start_age) || !IntegerVector::is_na(last_age))) {
    Rf_warningcall(R_NilValue, "Start and final ages cannot be used with matrix formats 1-3. Resetting these parameters....");
    start_age = 0;
    last_age = 0;
  }
  if (growthonly && repvalue) {
    Rf_warningcall(R_NilValue, "Option repvalue cannot be set to TRUE if growthonly is set to TRUE. Resetting repvalue to FALSE.");
    repvalue = false;
  }
  
  bool ipm_check_mid = stringcompare_simple(ipm_method, "mid", true);
  bool ipm_check_cdf = stringcompare_simple(ipm_method, "cdf", true);
  
  if (ipm_check_cdf) {
    ipm_method = "cdf";
  } else if (ipm_check_mid) {
    ipm_method = "midpoint";
  } else {
    Rf_warningcall(R_NilValue, "Option ipm_method is not understood. Will use cdf option.");
    ipm_method = "cdf";
  }

  // Vital rate models
  List msuite;
  RObject surmodl;
  RObject obsmodl;
  RObject sizmodl;
  RObject sibmodl;
  RObject sicmodl;
  RObject repmodl;
  RObject fecmodl;
  RObject jsurmodl;
  RObject jobsmodl;
  RObject jsizmodl;
  RObject jsibmodl;
  RObject jsicmodl;
  RObject jrepmodl;
  RObject jmatmodl;
  DataFrame pmnames;
  
  int modelcheck {0};
  NumericVector model1 = {1.0};
  NumericVector model0 = {0.0};
  CharacterVector modelnone {"none"};
  bool pmn_provided = false;
  
  if (modelsuite.isNotNull()) {
    Rcpp::List ms_intermediate(modelsuite);
    msuite = ms_intermediate;
    
    surmodl = msuite["survival_model"];
    obsmodl = msuite["observation_model"];
    sizmodl = msuite["size_model"];
    sibmodl = msuite["sizeb_model"];
    sicmodl = msuite["sizec_model"];
    repmodl = msuite["repstatus_model"];
    fecmodl = msuite["fecundity_model"];
    jsurmodl = msuite["juv_survival_model"];
    jobsmodl = msuite["juv_observation_model"];
    jsizmodl = msuite["juv_size_model"];
    jsibmodl = msuite["juv_sizeb_model"];
    jsicmodl = msuite["juv_sizec_model"];
    jrepmodl = msuite["juv_reproduction_model"];
    jmatmodl = msuite["juv_maturity_model"];
    
    pmnames = as<DataFrame>(msuite["paramnames"]);
    
    pmn_provided = true;
    modelcheck++;
    
  } else {
    if (surv_model.isNotNull()) {
      RObject sum_intermediate = RObject(surv_model);
      surmodl = sum_intermediate;
      modelcheck++;
    } else {
      surmodl = clone(model1);
    }
    if (obs_model.isNotNull()) {
      RObject obm_intermediate = RObject(obs_model);
      obsmodl = obm_intermediate;
      modelcheck++;
    } else {
      obsmodl = clone(model1);
    }
    if (size_model.isNotNull()) {
      RObject sim_intermediate = RObject(size_model);
      sizmodl = sim_intermediate;
      modelcheck++;
    } else {
      sizmodl = clone(model1);
    }
    if (sizeb_model.isNotNull()) {
      RObject sbm_intermediate = RObject(sizeb_model);
      sibmodl = sbm_intermediate;
      modelcheck++;
    } else {
      sibmodl = clone(model1);
    }
    if (sizec_model.isNotNull()) {
      RObject scm_intermediate = RObject(sizec_model);
      sicmodl = scm_intermediate;
      modelcheck++;
    } else {
      sicmodl = clone(model1);
    }
    if (repst_model.isNotNull()) {
      RObject rpm_intermediate = RObject(repst_model);
      repmodl = rpm_intermediate;
      modelcheck++;
    } else {
      repmodl = clone(model1);
    }
    if (fec_model.isNotNull()) {
      RObject fem_intermediate = RObject(fec_model);
      fecmodl = fem_intermediate;
      modelcheck++;
    } else {
      fecmodl = clone(model1);
    }
    
    if (jsurv_model.isNotNull()) {
      RObject jsum_intermediate = RObject(jsurv_model);
      jsurmodl = jsum_intermediate;
      modelcheck++;
    } else {
      jsurmodl = clone(model1);
    }
    if (jobs_model.isNotNull()) {
      RObject jobm_intermediate = RObject(jobs_model);
      jobsmodl = jobm_intermediate;
      modelcheck++;
    } else {
      jobsmodl = clone(model1);
    }
    if (jsize_model.isNotNull()) {
      RObject jsim_intermediate = RObject(jsize_model);
      jsizmodl = jsim_intermediate;
      modelcheck++;
    } else {
      jsizmodl = clone(model1);
    }
    if (jsizeb_model.isNotNull()) {
      RObject jsbm_intermediate = RObject(jsizeb_model);
      jsibmodl = jsbm_intermediate;
      modelcheck++;
    } else {
      jsibmodl = clone(model1);
    }
    if (jsizec_model.isNotNull()) {
      RObject jscm_intermediate = RObject(jsizec_model);
      jsicmodl = jscm_intermediate;
      modelcheck++;
    } else {
      jsicmodl = clone(model1);
    }
    if (jrepst_model.isNotNull()) {
      RObject jrpm_intermediate = RObject(jrepst_model);
      jrepmodl = jrpm_intermediate;
      modelcheck++;
    } else {
      jrepmodl = clone(model1);
    }
    if (jmatst_model.isNotNull()) {
      RObject jmat_intermediate = RObject(jmatst_model);
      jmatmodl = jmat_intermediate;
      modelcheck++;
    } else {
      jmatmodl = clone(model1);
    }
    
    if (paramnames.isNotNull()) {
      DataFrame pmn_intermediate = DataFrame(paramnames);
      pmnames = pmn_intermediate;
      
      pmn_provided = true;
    }
  }
  if (modelcheck == 0) {
    throw Rcpp::exception("This function requires a lefkoMod object or vital rate models.", false);
  }
  if (!pmn_provided) {
    throw Rcpp::exception("A paramnames object is required if a lefkoMod object is not supplied.", false);
  }
  
  // Stageframe
  DataFrame sframe;
  RObject maingroups;
  
  if (stageframe.isNotNull()) {
    DataFrame sf_intermediate = DataFrame(stageframe);
    sframe = sf_intermediate;
    
    CharacterVector maingroups_int = as<CharacterVector>(sframe["group"]);
    CharacterVector maingroups_int2 = unique(maingroups_int);
    
    CharacterVector maingroups_sorted = stringsort(maingroups_int2);
    maingroups = as<RObject>(maingroups_sorted);
    
  } else if (format < 5) {
    throw Rcpp::exception("A stageframe is required for all MPM formats except Leslie MPMs.", false);
  } else {
    sframe = R_NilValue;
    
    CharacterVector maingroups_sorted = {"0"};
    maingroups = as<RObject>(maingroups_sorted);
  }
  
  // Demographic data
  int no_vars = data.length();
  int used_yearcol {-1};
  int used_patchcol {-1};
  int used_agecol {-1};
  int used_indacol {-1};
  int used_indbcol {-1};
  int used_indccol {-1};
  
  CharacterVector modelparam_names = pmnames["modelparams"];
  std::string year2var = as<std::string>(modelparam_names(0));
  std::string patchvar = as<std::string>(modelparam_names(2));
  std::string agevar = as<std::string>(modelparam_names(21));
  std::string indcova2var = as<std::string>(modelparam_names(23));
  std::string indcova1var = as<std::string>(modelparam_names(24));
  std::string indcovb2var = as<std::string>(modelparam_names(25));
  std::string indcovb1var = as<std::string>(modelparam_names(26));
  std::string indcovc2var = as<std::string>(modelparam_names(27));
  std::string indcovc1var = as<std::string>(modelparam_names(28));
  std::string group2var = as<std::string>(modelparam_names(29));
  std::string group1var = as<std::string>(modelparam_names(30));
  
  CharacterVector data_names = data.attr("names");
  
  for (int i = 0; i < no_vars; i++) {
    if (stringcompare_hard(as<std::string>(data_names(i)), year2var)) {
      used_yearcol = i;
    }
    if (stringcompare_hard(as<std::string>(data_names(i)), patchvar)) {
      used_patchcol = i;
    }
    if (stringcompare_hard(as<std::string>(data_names(i)), agevar)) {
      used_agecol = i;
    }
    if (stringcompare_hard(as<std::string>(data_names(i)), indcova2var)) {
      used_indacol = i;
    }
    if (stringcompare_hard(as<std::string>(data_names(i)), indcovb2var)) {
      used_indbcol = i;
    }
    if (stringcompare_hard(as<std::string>(data_names(i)), indcovc2var)) {
      used_indccol = i;
    }
  }
  if (used_yearcol < 0 || used_yearcol > (no_vars - 1)) {
    throw Rcpp::exception("Correct variable in dataset for time t not found. Is paramnames object correct?", 
      false);
  }
  IntegerVector all_years = as<IntegerVector>(data[used_yearcol]);
  IntegerVector all_years_x = unique(all_years);
  IntegerVector mainyears_int = int_sort(all_years_x);
  NumericVector mainyears = as<NumericVector>(mainyears_int);
  
  CharacterVector mainpatches;
  
  if (used_patchcol > -1) {
    CharacterVector all_patches = as<CharacterVector>(data[used_patchcol]);
    CharacterVector mainpatches_int = unique(all_patches);
    
    mainpatches = stringsort(mainpatches_int);
  } else {
    mainpatches = {NA_STRING};
  }
  
  IntegerVector mainages;
  IntegerVector actualages;
  int age_limit {0};
  
  if (format > 3) { // Age and age-by-stage cases
    if (used_agecol > -1) {
      IntegerVector all_ages = as<IntegerVector>(data[used_agecol]);
      IntegerVector mainages_pre = unique(all_ages);
      mainages = int_sort(mainages_pre);
    }
    age_limit = max(mainages) + 1;
    
    if (IntegerVector::is_na(start_age)) {
      if (prebreeding) {
        start_age = 1;
      } else {
        start_age = 0;
      }
    }
    if (IntegerVector::is_na(last_age)) {
      last_age = max(mainages) + 1;
    }
    if (IntegerVector::is_na(fecage_min)) {
      fecage_min = min(mainages);
    }
    if (IntegerVector::is_na(fecage_max)) {
      fecage_max = last_age;
    }
    
    if (start_age > age_limit || last_age > age_limit) {
      Rf_warningcall(R_NilValue, "Entered start_age or last_age is beyond what is found in the dataset.");
    }
    if (fecage_min > age_limit || fecage_max > age_limit) {
      Rf_warningcall(R_NilValue, "Entered fecage_min or fecage_max is beyond what is found in the dataset.");
    }
    
    if (last_age < (start_age + 1)) {
      throw Rcpp::exception("Please set last_age to be greater than start_age.");
    }
    if (fecage_max < fecage_min) {
      throw Rcpp::exception("Please set fecage_max to be greater than or equal to fecage_min.");
    }
  }
  
  // Ind_names objects
  RObject inda_names;
  RObject indb_names;
  RObject indc_names;
  CharacterVector inda_names_ch;
  CharacterVector indb_names_ch;
  CharacterVector indc_names_ch;
  
  if (used_indacol > -1) {
    if (random_inda) {
      CharacterVector inda_data = as<CharacterVector>(data[used_indacol]);
      CharacterVector unique_a = unique(inda_data);
      CharacterVector sorted_a = stringsort(unique_a);
      inda_names_ch = sorted_a;
      
      inda_names = as<RObject>(sorted_a);
      
    } else {
      IntegerVector notrandom_a = {0};
      inda_names = as<RObject>(notrandom_a);
    }
  } else {
    IntegerVector notrandom_a = {0};
    inda_names = as<RObject>(notrandom_a);
    random_inda = false;
  }
  
  if (used_indbcol > -1) {
    if (random_indb) {
      CharacterVector indb_data = as<CharacterVector>(data[used_indbcol]);
      CharacterVector unique_b = unique(indb_data);
      CharacterVector sorted_b = stringsort(unique_b);
      indb_names_ch = sorted_b;
      
      indb_names = as<RObject>(sorted_b);
      
    } else {
      IntegerVector notrandom_b = {0};
      indb_names = as<RObject>(notrandom_b);
    }
  } else {
    IntegerVector notrandom_b = {0};
    indb_names = as<RObject>(notrandom_b);
    random_indb = false;
  }
  
  if (used_indccol > -1) {
    if (random_indc) {
      CharacterVector indc_data = as<CharacterVector>(data[used_indccol]);
      CharacterVector unique_c = unique(indc_data);
      CharacterVector sorted_c = stringsort(unique_c);
      indc_names_ch = sorted_c;
      
      indc_names = as<RObject>(sorted_c);
      
    } else {
      IntegerVector notrandom_c = {0};
      indc_names = as<RObject>(notrandom_c);
    }
  } else {
    IntegerVector notrandom_c = {0};
    indc_names = as<RObject>(notrandom_c);
    random_indc = false;
  }
  
  // Check years and patches entries
  List all_years_topull(nreps);
  NumericVector years_topull;
  NumericMatrix years_projected (times, nreps);
  int num_years = mainyears.length();
  if (year.isNotNull()) {
    NumericVector year_int (year);
    
    NumericVector years_unmatched = setdiff(year_int, mainyears);
    if (years_unmatched.length() > 0) {
      throw Rcpp::exception("Some input year values do not match years documented in the dataset.");
    }
    if (stochastic) throw Rcpp::exception("If option year is set, then projection cannot be stochastic.");
    years_topull = year_int;
  } else if (!stochastic) {
    NumericVector year_int = clone(mainyears);
    years_topull = year_int;
    Rf_warningcall(R_NilValue, "Option year not set, so will cycle through existing years.");
  }
  
  NumericVector twinput;
  if (tweights.isNotNull()) {
    twinput = as<NumericVector>(tweights);
    if (twinput.length() != num_years) {
      throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the dataset.", false);
    }
  } else {
    NumericVector tw_temp(num_years);
    tw_temp.fill(1);
    twinput = tw_temp;
  }
  
  if (stochastic) {
    years_topull = Rcpp::RcppArmadillo::sample(mainyears, times * nreps, true, twinput);
  }
  
  CharacterVector patches_topull;
  CharacterVector patches_projected (times);
  int chosenpatch {1};
  if (patch.isNotNull()) {
    CharacterVector patch_int (patch);
    
    CharacterVector patches_unmatched = setdiff(patch_int, mainpatches);
    if (patches_unmatched.length() > 0) {
      throw Rcpp::exception("Some input patch values do not match patches documented in the dataset.");
    }
    patches_topull = patch_int;
    int crazy_patch {-1};
    if (patch_int.length() == 1) {
      for (int i = 0; i < mainpatches.length(); i++) {
        if (stringcompare_hard(as<std::string>(patch_int(0)), as<std::string>(mainpatches(i)))) {
          crazy_patch = i;
        }
      }
      if (crazy_patch != -1) chosenpatch = crazy_patch + 1;
    }
  } else {
    CharacterVector patch_one (1);
    patch_one(0) = mainpatches(0);
    patches_topull = patch_one;
    Rf_warningcall(R_NilValue, "Option patch not set, so will set to first patch/population.");
  }
  
  // Handle spatial density vector
  NumericVector spdensity_topull;
  NumericVector spdensity_projected (times);
  if (sp_density.isNotNull()) {
    NumericVector spdensity_int (sp_density);
    spdensity_topull = spdensity_int;
  } else{
    NumericVector spdensity_int = {0.0};
    spdensity_topull = spdensity_int;
  }
  
  // Ind_terms data frame
  NumericVector f2_inda_values(times);
  NumericVector f2_indb_values(times);
  NumericVector f2_indc_values(times);
  NumericVector f1_inda_values(times);
  NumericVector f1_indb_values(times);
  NumericVector f1_indc_values(times);
  CharacterVector r2_inda_values(times);
  CharacterVector r2_indb_values(times);
  CharacterVector r2_indc_values(times);
  CharacterVector r1_inda_values(times);
  CharacterVector r1_indb_values(times);
  CharacterVector r1_indc_values(times);
  
  NumericVector f_inda_topull;
  NumericVector f_indb_topull;
  NumericVector f_indc_topull;
  CharacterVector r_inda_topull;
  CharacterVector r_indb_topull;
  CharacterVector r_indc_topull;
  
  bool lessthan_warning = false;
  bool greaterthan_warning = false;
  
  if (ind_terms.isNotNull()) {
    DataFrame it_ROint = RObject(ind_terms);
    
    RObject inda_whatever;
    RObject indb_whatever;
    RObject indc_whatever;
    
    if (is<DataFrame>(it_ROint)) {
      DataFrame it_intermediate = DataFrame(it_ROint);
      int it_size = it_intermediate.size();
      if (it_size != 3) {
        throw Rcpp::exception("Data frame ind_terms should have 3 columns and times rows.", false);
      }
      
      inda_whatever = as<RObject>(it_intermediate[0]);
      indb_whatever = as<RObject>(it_intermediate[1]);
      indc_whatever = as<RObject>(it_intermediate[2]);
      
    } else if (is<NumericMatrix>(it_ROint)) {
      NumericMatrix it_intermediate = NumericMatrix(it_ROint);
      int it_rows = it_intermediate.nrow();
      int it_cols = it_intermediate.ncol();
      
      if (it_rows != 3 && it_cols != 3) {
        throw Rcpp::exception("Numeric matrix ind_terms should have 3 columns and times rows.", false);
      }
      
      if (it_cols == 3) {
        NumericVector chuck1 = it_intermediate(_, 0);
        NumericVector chuck2 = it_intermediate(_, 1);
        NumericVector chuck3 = it_intermediate(_, 2);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
      } else if (it_rows == 3) {
        NumericVector chuck1 = it_intermediate(0, _);
        NumericVector chuck2 = it_intermediate(1, _);
        NumericVector chuck3 = it_intermediate(2, _);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
      }
    } else if (is<CharacterMatrix>(it_ROint)) {
      CharacterMatrix it_intermediate = CharacterMatrix(it_ROint);
      int it_rows = it_intermediate.nrow();
      int it_cols = it_intermediate.ncol();
      
      if (it_rows != 3 && it_cols != 3) {
        throw Rcpp::exception("String matrix ind_terms should have 3 columns and times rows.", false);
      }
      
      if (it_cols == 3) {
        CharacterVector chuck1 = it_intermediate(_, 0);
        CharacterVector chuck2 = it_intermediate(_, 1);
        CharacterVector chuck3 = it_intermediate(_, 2);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
       } else if (it_rows == 3) {
        CharacterVector chuck1 = it_intermediate(0, _);
        CharacterVector chuck2 = it_intermediate(1, _);
        CharacterVector chuck3 = it_intermediate(2, _);
        inda_whatever = as<RObject>(chuck1);
        indb_whatever = as<RObject>(chuck2);
        indc_whatever = as<RObject>(chuck3);
       }
    }
    
    if (is<NumericVector>(inda_whatever)) {
      if (random_inda) {
        Rf_warningcall(R_NilValue, "Indcov a appears to be numeric. Will assume that random_inda = FALSE. To alter this behavior, please convert indcov a into a character vector.");
        random_inda = false;
      }
      NumericVector inda_another_int = as<NumericVector>(inda_whatever);
      
      int check_len = inda_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_inda_topull = inda_another_int;
      r_inda_topull = {"none"};
      
    } else if (is<CharacterVector>(inda_whatever)) {
      if (!random_inda) {
        Rf_warningcall(R_NilValue, "Indcov a appears to be categorical. Will assume that random_inda = TRUE. To alter this behavior, please convert indcov a into a numeric vector.");
        random_inda = true;
      }
      CharacterVector inda_another_int = as<CharacterVector>(inda_whatever);
      
      int inda_inputlength = inda_another_int.length();
      int maininda_length = inda_names_ch.length();
      
      for (int i = 0; i < inda_inputlength; i++) {
        for (int j = 0; j < maininda_length; j++) {
          if (!stringcompare_hard(as<std::string>(inda_another_int(i)), as<std::string>(inda_names_ch(j)))) {
            throw Rcpp::exception("Some input ind cov a values do not match categories documented in the dataset.");
          }
        }
      }
      int check_len = inda_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      r_inda_topull = inda_another_int;
      f_inda_topull = {0};
    }
    
    if (is<NumericVector>(indb_whatever)) {
      if (random_indb) {
        Rf_warningcall(R_NilValue, "Indcov b appears to be numeric. Will assume that random_indb = FALSE. To alter this behavior, please convert indcov b into a character vector.");
        random_indb = false;
      }
      NumericVector indb_another_int = as<NumericVector>(indb_whatever);
      
      int check_len = indb_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_indb_topull = indb_another_int;
      r_indb_topull = {"none"};
      
    } else if (is<CharacterVector>(indb_whatever)) {
      if (!random_indb) {
        Rf_warningcall(R_NilValue, "Indcov b appears to be categorical. Will assume that random_indb = TRUE. To alter this behavior, please convert indcov b into a numeric vector.");
        random_indb = true;
      }
      CharacterVector indb_another_int = as<CharacterVector>(indb_whatever);
      
      int indb_inputlength = indb_another_int.length();
      int mainindb_length = indb_names_ch.length();
      
      for (int i = 0; i < indb_inputlength; i++) {
        for (int j = 0; j < mainindb_length; j++) {
          if (!stringcompare_hard(as<std::string>(indb_another_int(i)), as<std::string>(indb_names_ch(j)))) {
            throw Rcpp::exception("Some input ind cov b values do not match categories documented in the dataset.");
          }
        }
      }
      int check_len = indb_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      r_indb_topull = indb_another_int;
      f_indb_topull = {0};
    }
    
    if (is<NumericVector>(indc_whatever)) {
      if (random_indc) {
        Rf_warningcall(R_NilValue, "Indcov c appears to be numeric. Will assume that random_indc = FALSE. To alter this behavior, please convert indcov c into a character vector.");
        random_indc = false;
      }
      NumericVector indc_another_int = as<NumericVector>(indc_whatever);
      
      int check_len = indc_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      f_indc_topull = indc_another_int;
      r_indc_topull = {"none"};
      
    } else if (is<CharacterVector>(indc_whatever)) {
      if (!random_indc) {
        Rf_warningcall(R_NilValue, "Indcov c appears to be categorical. Will assume that random_indc = TRUE. To alter this behavior, please convert indcov c into a numeric vector.");
        random_indc = true;
      }
      CharacterVector indc_another_int = as<CharacterVector>(indc_whatever);
      
      int indc_inputlength = indc_another_int.length();
      int mainindc_length = indc_names_ch.length();
      
      for (int i = 0; i < indc_inputlength; i++) {
        for (int j = 0; j < mainindc_length; j++) {
          if (!stringcompare_hard(as<std::string>(indc_another_int(i)), as<std::string>(indc_names_ch(j)))) {
            throw Rcpp::exception("Some input ind cov c values do not match categories documented in the dataset.");
          }
        }
      }
      int check_len = indc_another_int.length();
      if (check_len > times) greaterthan_warning = true;
      if (check_len < times) lessthan_warning = true;
      
      r_indc_topull = indc_another_int;
      f_indc_topull = {0};
    }
    
  } else {
    f_inda_topull = clone(model0);
    r_inda_topull = clone(modelnone);
    f_indb_topull = clone(model0);
    r_indb_topull = clone(modelnone);
    f_indc_topull = clone(model0);
    r_indc_topull = clone(modelnone);
  }
  
  if (greaterthan_warning) {
    Rf_warningcall(R_NilValue, "More values of individual covariates have been supplied than required, so some will be cut.");
  }
  if (lessthan_warning) {
    Rf_warningcall(R_NilValue, "Fewer values of individual covariates have been supplied than required, so input values will be cycled.");
  }
  
  // dev_terms data frame or matrix
  NumericVector sur_dev_values(times);
  NumericVector obs_dev_values(times);
  NumericVector siz_dev_values(times);
  NumericVector sib_dev_values(times);
  NumericVector sic_dev_values(times);
  NumericVector rep_dev_values(times);
  NumericVector fec_dev_values(times);
  NumericVector jsur_dev_values(times);
  NumericVector jobs_dev_values(times);
  NumericVector jsiz_dev_values(times);
  NumericVector jsib_dev_values(times);
  NumericVector jsic_dev_values(times);
  NumericVector jrep_dev_values(times);
  NumericVector jmat_dev_values(times);
  
  NumericVector surv_dev_extracted;
  NumericVector obs_dev_extracted;
  NumericVector size_dev_extracted;
  NumericVector sizeb_dev_extracted;
  NumericVector sizec_dev_extracted;
  NumericVector repst_dev_extracted;
  NumericVector fec_dev_extracted;
  NumericVector jsurv_dev_extracted;
  NumericVector jobs_dev_extracted;
  NumericVector jsize_dev_extracted;
  NumericVector jsizeb_dev_extracted;
  NumericVector jsizec_dev_extracted;
  NumericVector jrepst_dev_extracted;
  NumericVector jmatst_dev_extracted;
  int veclimits {0};
  
  bool lessthan_warning_dev = false;
  bool greaterthan_warning_dev = false;
  
  if (dev_terms.isNotNull()) {
    RObject dt_intermediate = RObject(dev_terms);
    
    if (is<NumericMatrix>(dt_intermediate)) {
      NumericMatrix dt_mat = as<NumericMatrix>(dt_intermediate);
      int dt_rows = dt_mat.nrow();
      int dt_cols = dt_mat.ncol();
      
      if (dt_rows != 14 && dt_cols != 14) {
        throw Rcpp::exception("Deviation term matrix must have 14 columns.");
      }
      
      if (dt_rows == 14) {
        surv_dev_extracted = dt_mat(0, _);
        obs_dev_extracted = dt_mat(1, _);
        size_dev_extracted = dt_mat(2, _);
        sizeb_dev_extracted = dt_mat(3, _);
        sizec_dev_extracted = dt_mat(4, _);
        repst_dev_extracted = dt_mat(5, _);
        fec_dev_extracted = dt_mat(6, _);
        jsurv_dev_extracted = dt_mat(7, _);
        jobs_dev_extracted = dt_mat(8, _);
        jsize_dev_extracted = dt_mat(9, _);
        jsizeb_dev_extracted = dt_mat(10, _);
        jsizec_dev_extracted = dt_mat(11, _);
        jrepst_dev_extracted = dt_mat(12, _);
        jmatst_dev_extracted = dt_mat(13, _);
        
        veclimits = dt_cols;
        
        if (dt_cols < times) lessthan_warning_dev = true;
        if (dt_cols > times) greaterthan_warning_dev = true;
        
      } else {
        surv_dev_extracted = dt_mat(_, 0);
        obs_dev_extracted = dt_mat(_, 1);
        size_dev_extracted = dt_mat(_, 2);
        sizeb_dev_extracted = dt_mat(_, 3);
        sizec_dev_extracted = dt_mat(_, 4);
        repst_dev_extracted = dt_mat(_, 5);
        fec_dev_extracted = dt_mat(_, 6);
        jsurv_dev_extracted = dt_mat(_, 7);
        jobs_dev_extracted = dt_mat(_, 8);
        jsize_dev_extracted = dt_mat(_, 9);
        jsizeb_dev_extracted = dt_mat(_, 10);
        jsizec_dev_extracted = dt_mat(_, 11);
        jrepst_dev_extracted = dt_mat(_, 12);
        jmatst_dev_extracted = dt_mat(_, 13);
        
        veclimits = dt_rows;
        
        if (dt_rows < times) lessthan_warning_dev = true;
        if (dt_rows > times) greaterthan_warning_dev = true;
      }
      
    } else if (is<DataFrame>(dt_intermediate)) {
      DataFrame dt_frame = as<DataFrame>(dt_intermediate);
      int dt_vars = dt_frame.size();
      
      if (dt_vars != 14) {
        throw Rcpp::exception("Deviation term data frame must have 14 numeric columns.");
      }
      
      RObject surv_dev_a = as<RObject>(dt_frame[0]);
      RObject obs_dev_a = as<RObject>(dt_frame[1]);
      RObject size_dev_a = as<RObject>(dt_frame[2]);
      RObject sizeb_dev_a = as<RObject>(dt_frame[3]);
      RObject sizec_dev_a = as<RObject>(dt_frame[4]);
      RObject repst_dev_a = as<RObject>(dt_frame[5]);
      RObject fec_dev_a = as<RObject>(dt_frame[6]);
      RObject jsurv_dev_a = as<RObject>(dt_frame[7]);
      RObject jobs_dev_a = as<RObject>(dt_frame[8]);
      RObject jsize_dev_a = as<RObject>(dt_frame[9]);
      RObject jsizeb_dev_a = as<RObject>(dt_frame[10]);
      RObject jsizec_dev_a = as<RObject>(dt_frame[11]);
      RObject jrepst_dev_a = as<RObject>(dt_frame[12]);
      RObject jmatst_dev_a = as<RObject>(dt_frame[13]);
      
      if (is<NumericVector>(surv_dev_a)) {
        surv_dev_extracted = as<NumericVector>(surv_dev_a);
        veclimits = surv_dev_extracted.length();
        
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(obs_dev_a)) {
        obs_dev_extracted = as<NumericVector>(obs_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(size_dev_a)) {
        size_dev_extracted = as<NumericVector>(size_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(sizeb_dev_a)) {
        sizeb_dev_extracted = as<NumericVector>(sizeb_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(sizec_dev_a)) {
        sizec_dev_extracted = as<NumericVector>(sizec_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(repst_dev_a)) {
        repst_dev_extracted = as<NumericVector>(repst_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(fec_dev_a)) {
        fec_dev_extracted = as<NumericVector>(fec_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jsurv_dev_a)) {
        jsurv_dev_extracted = as<NumericVector>(jsurv_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jobs_dev_a)) {
        jobs_dev_extracted = as<NumericVector>(jobs_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jsize_dev_a)) {
        jsize_dev_extracted = as<NumericVector>(jsize_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jsizeb_dev_a)) {
        jsizeb_dev_extracted = as<NumericVector>(jsizeb_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jsizec_dev_a)) {
        jsizec_dev_extracted = as<NumericVector>(jsizec_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jrepst_dev_a)) {
        jrepst_dev_extracted = as<NumericVector>(jrepst_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      if (is<NumericVector>(jmatst_dev_a)) {
        jmatst_dev_extracted = as<NumericVector>(jmatst_dev_a);
      } else {
        throw Rcpp::exception("Some dev_terms columns appear not to be numeric.");
      }
      
    } else {
      throw Rcpp::exception("Option dev_terms must be input as a data frame of 14 numeric variables.");
    }
  } else {
    veclimits = 1;
    surv_dev_extracted = clone(model0);
    obs_dev_extracted = clone(model0);
    size_dev_extracted = clone(model0);
    sizeb_dev_extracted = clone(model0);
    sizec_dev_extracted = clone(model0);
    repst_dev_extracted = clone(model0);
    fec_dev_extracted = clone(model0);
    jsurv_dev_extracted = clone(model0);
    jobs_dev_extracted = clone(model0);
    jsize_dev_extracted = clone(model0);
    jsizeb_dev_extracted = clone(model0);
    jsizec_dev_extracted = clone(model0);
    jrepst_dev_extracted = clone(model0);
    jmatst_dev_extracted = clone(model0);
  }
  
  if (greaterthan_warning_dev) {
    Rf_warningcall(R_NilValue, "More values of intercept deviations have been supplied than required, so some will be cut.");
  }
  if (lessthan_warning_dev) {
    Rf_warningcall(R_NilValue, "Fewer values of intercept deviations have been supplied than required, so input values will be cycled.");
  }
  
  // Main for loop adjusting lengths of input vectors
  int year_counter {0};
  int patch_counter {0};
  int spdensity_counter {0};
  int finda_counter {0};
  int findb_counter {0};
  int findc_counter {0};
  int rinda_counter {0};
  int rindb_counter {0};
  int rindc_counter {0};
  int dev_counter {0};
  
  int year_limit = years_topull.length() / nreps;
  int patch_limit = patches_topull.length();
  int spdensity_limit = spdensity_topull.length();
  int finda_limit = f_inda_topull.length();
  int findb_limit = f_indb_topull.length();
  int findc_limit = f_indc_topull.length();
  int rinda_limit = r_inda_topull.length();
  int rindb_limit = r_indb_topull.length();
  int rindc_limit = r_indc_topull.length();
  
  for (int i = 0; i < times; i++) {
    if (year_counter == year_limit) year_counter = 0;
    if (patch_counter == patch_limit) patch_counter = 0;
    if (spdensity_counter == spdensity_limit) spdensity_counter = 0;
    if (finda_counter >= finda_limit) finda_counter = 0;
    if (rinda_counter >= rinda_limit) rinda_counter = 0;
    if (findb_counter >= findb_limit) findb_counter = 0;
    if (rindb_counter >= rindb_limit) rindb_counter = 0;
    if (findc_counter >= findc_limit) findc_counter = 0;
    if (rindc_counter >= rindc_limit) rindc_counter = 0;
    if (dev_counter >= veclimits) dev_counter = 0;
    
    for (int j = 0; j < nreps; j++) {
      years_projected(i, j) = years_topull(year_counter + times * j);
    }
    
    patches_projected(i) = patches_topull(patch_counter);
    spdensity_projected(i) = spdensity_topull(spdensity_counter);
    if (NumericVector::is_na(spdensity_projected(i))) spdensity_projected(i) = 0.0;
    f2_inda_values(i) = f_inda_topull(finda_counter);
    r2_inda_values(i) = r_inda_topull(rinda_counter);
    f2_indb_values(i) = f_indb_topull(findb_counter);
    r2_indb_values(i) = r_indb_topull(rindb_counter);
    f2_indc_values(i) = f_indc_topull(findc_counter);
    r2_indc_values(i) = r_indc_topull(rindc_counter);
    
    if (i > 0) {
      f1_inda_values(i) = f2_inda_values(i-1);
      r1_inda_values(i) = r2_inda_values(i-1);
      f1_indb_values(i) = f2_indb_values(i-1);
      r1_indb_values(i) = r2_indb_values(i-1);
      f1_indc_values(i) = f2_indc_values(i-1);
      r1_indc_values(i) = r2_indc_values(i-1);
      
    } else {
      r1_inda_values(i) = modelnone(0);
      r1_indb_values(i) = modelnone(0);
      r1_indc_values(i) = modelnone(0);
    }
    
    sur_dev_values(i) = surv_dev_extracted(dev_counter);
    obs_dev_values(i) = obs_dev_extracted(dev_counter);
    siz_dev_values(i) = size_dev_extracted(dev_counter);
    sib_dev_values(i) = sizeb_dev_extracted(dev_counter);
    sic_dev_values(i) = sizec_dev_extracted(dev_counter);
    rep_dev_values(i) = repst_dev_extracted(dev_counter);
    fec_dev_values(i) = fec_dev_extracted(dev_counter);
    jsur_dev_values(i) = jsurv_dev_extracted(dev_counter);
    jobs_dev_values(i) = jobs_dev_extracted(dev_counter);
    jsiz_dev_values(i) = jsize_dev_extracted(dev_counter);
    jsib_dev_values(i) = jsizeb_dev_extracted(dev_counter);
    jsic_dev_values(i) = jsizec_dev_extracted(dev_counter);
    jrep_dev_values(i) = jrepst_dev_extracted(dev_counter);
    jmat_dev_values(i) = jmatst_dev_extracted(dev_counter);
    
    year_counter++;
    patch_counter++;
    spdensity_counter++;
    finda_counter++;
    rinda_counter++;
    findb_counter++;
    rindb_counter++;
    findc_counter++;
    rindc_counter++;
    dev_counter++;
  }
  
  // Allstages
  DataFrame new_stageframe;
  arma::mat new_repmatrix;
  DataFrame new_ovtable;
  DataFrame allstages;
  
  if (format < 5) {
    bool agemat = false;
    bool historical = false;
    int ehrlen {1};
    int style {0};
    int filter {1};
    
    if (format == 2) ehrlen = 2;
    if (format == 3) style = 1;
    if (format == 4) {
      agemat = true;
      style = 2;
      filter = 2;
    }
    if (format < 3) historical = true;
    
    List melchett = sf_reassess(sframe, supplement, overwrite, repmatrix,
      agemat, historical, ehrlen);
    new_stageframe = as<DataFrame>(melchett["stageframe"]);
    new_repmatrix = as<arma::mat>(melchett["repmatrix"]);
    new_ovtable = as<DataFrame>(melchett["ovtable"]);
    
    // the old pizzle needs to be called
    DataFrame allstages_pre = theoldpizzle(new_stageframe, new_ovtable,
      new_repmatrix, start_age, last_age, ehrlen, style, cont, filter);
    allstages = allstages_pre;
    
  } else {
    DataFrame melchett = sf_leslie(start_age, last_age, fecage_min, fecage_max, cont);
    new_stageframe = melchett;
    allstages = melchett;
    
    CharacterVector maingroups_ch = {"0"};
    maingroups = as<RObject>(maingroups_ch);
    actualages = seq(start_age, last_age);
  }
  
  double maxsize {0.0};
  double maxsizeb {0.0};
  double maxsizec {0.0};
  
  if (format < 5) {
    NumericVector size3 = allstages["size3"];
    NumericVector size2n = allstages["size2n"];
    NumericVector size2o = allstages["size2o"];
    NumericVector sizeb3 = allstages["sizeb3"];
    NumericVector sizeb2n = allstages["sizeb2n"];
    NumericVector sizeb2o = allstages["sizeb2o"];
    NumericVector sizec3 = allstages["sizec3"];
    NumericVector sizec2n = allstages["sizec2n"];
    NumericVector sizec2o = allstages["sizec2o"];
    
    NumericVector maxveca = {max(size3), max(size2n), max(size2o)}; // What about NAs?
    NumericVector maxvecb = {max(sizeb3), max(sizeb2n), max(sizeb2o)};
    NumericVector maxvecc = {max(sizec3), max(sizec2n), max(sizec2o)};
    
    maxsize = max(maxveca); // What about NAs?
    maxsizeb = max(maxvecb);
    maxsizec = max(maxvecc);
  }
  
  // modelextract proxy lists
  List surv_proxy = modelextract(surmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List obs_proxy = modelextract(obsmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List size_proxy = modelextract(sizmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List sizeb_proxy = modelextract(sibmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List sizec_proxy = modelextract(sicmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List repst_proxy = modelextract(repmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List fec_proxy = modelextract(fecmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jsurv_proxy = modelextract(jsurmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jobs_proxy = modelextract(jobsmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jsize_proxy = modelextract(jsizmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jsizeb_proxy = modelextract(jsibmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jsizec_proxy = modelextract(jsicmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jrepst_proxy = modelextract(jrepmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  List jmatst_proxy = modelextract(jmatmodl, pmnames, mainyears, mainpatches,
    maingroups, inda_names, indb_names, indc_names);
  
  // Main projection set-up
  int yearnumber {0};
  int patchnumber {0};
  NumericVector used_devs;
  
  List madsexmadrigal_oneyear;
  List madsexmadrigal_forward;
  arma::mat Amat;
  arma::mat Umat;
  arma::mat Fmat;
  
  used_devs = {sur_dev_values(0), obs_dev_values(0), siz_dev_values(0), sib_dev_values(0),
    sic_dev_values(0), rep_dev_values(0), fec_dev_values(0), jsur_dev_values(0),
    jobs_dev_values(0), jsiz_dev_values(0), jsib_dev_values(0), jsic_dev_values(0),
    jrep_dev_values(0), jmat_dev_values(0)};
  
  if (format < 5) {
    madsexmadrigal_oneyear = jerzeibalowski(allstages, new_stageframe,
      format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
      repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
      jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
      f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
      r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
      r2_indc_values, r1_indc_values, used_devs, spdensity_projected(0),
      repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
      yearnumber, patchnumber, exp_tol, theta_tol, ipm_method, err_check,
      true);
    
  } else {
    madsexmadrigal_oneyear = motherbalowski(actualages, new_stageframe,
      surv_proxy, fec_proxy, f2_inda_values, f1_inda_values, f2_indb_values,
      f1_indb_values, f2_indc_values, f1_indc_values, r2_inda_values,
      r1_inda_values, r2_indb_values, r1_indb_values, r2_indc_values,
      r1_indc_values, sur_dev_values(0), fec_dev_values(0), spdensity_projected(0),
      repmod, last_age, false, yearnumber, patchnumber, exp_tol, theta_tol, true);
  }
  
  Umat = as<arma::mat>(madsexmadrigal_oneyear["U"]);
  int meanmatrows = Umat.n_rows;
  
  DataFrame ahstages;
  if (format < 5) {
    ahstages = sframe;
  } else {
    ahstages = new_stageframe;
  }
  
  DataFrame hstages;
  if (format < 3) {
    hstages = hst_maker(ahstages);
  } else {
    hstages = R_NilValue;
  }
  
  DataFrame agestages;
  if (format == 4) {
    agestages = age_maker(ahstages, start_age, last_age);
  } else {
    agestages = R_NilValue;
  }
  
  // Here we will check if the matrix is large and sparse
  int sparse_switch {0};
  int test_elems = Umat.n_elem;
  arma::uvec nonzero_elems = find(Umat);
  int all_nonzeros = nonzero_elems.n_elem;
  double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
  if (sparse_check <= 0.5 && meanmatrows > 100) {
    sparse_switch = 1;
  }
  
  // Start data frame and vector
  arma::vec startvec;
  if(start_frame.isNotNull()) {
    if (!is<DataFrame>(start_frame)) {
      throw Rcpp::exception("Option start_frame must be a data frame created with function start_input().");
    }
    Rcpp::DataFrame start_thru(start_frame);
    startvec.set_size(meanmatrows);
    startvec.zeros();
    arma::uvec start_elems = as<arma::uvec>(start_thru["row_num"]);
    start_elems = start_elems - 1;
    arma::vec start_values = as<arma::vec>(start_thru["value"]);
    
    if (start_elems.max() > (meanmatrows - 1)) {
      throw Rcpp::exception("Start vector input frame includes element indices too high for this MPM.",
        false);
    }
    for (int i = 0; i < start_elems.n_elem; i++) {
      startvec(start_elems(i)) = start_values(i);
    }
    
  } else if (start_vec.isNotNull()) {
    startvec = as<arma::vec>(start_vec);
    if (startvec.n_elem != meanmatrows) {
      throw Rcpp::exception("Start vector must be the same length as the number of rows in each matrix.",
        false);
    }
    
  } else {
    startvec.set_size(meanmatrows);
    startvec.ones();
  }
  
  // Density dependence inputs
  Rcpp::DataFrame dens_input;
  List dens_index;
  double pop_size {0};
  double changing_element_U {0.0};
  double changing_element_F {0.0};
  double changing_colsum {0.0};
  
  int time_delay {1};
  int dens_switch {0};
  bool warn_trigger_neg = false;
  bool warn_trigger_1 = false;
  
  arma::uvec dyn_index321;
  arma::uvec dyn_index_col;
  arma::uvec dyn_style;
  arma::vec dyn_alpha;
  arma::vec dyn_beta;
  arma::uvec dyn_delay;
  arma::uvec dyn_type;
  int n_dyn_elems {0};
  
  if (density.isNotNull()) {
    if (!is<DataFrame>(density)) {
      throw Rcpp::exception("Option density must be a data frame created with function density_input().");
    }
    Rcpp::DataFrame dens_thru(density);
    dens_input = dens_thru;
    dens_switch = 1;
    
    Rcpp::StringVector di_stage3 = as<StringVector>(dens_input["stage3"]);
    Rcpp::StringVector di_stage2 = as<StringVector>(dens_input["stage2"]);
    Rcpp::StringVector di_stage1 = as<StringVector>(dens_input["stage1"]);
    int di_size = di_stage3.length();
    
    if (format < 3) { // Historical matrices
      StringVector stage3 = as<StringVector>(hstages["stage_2"]);
      StringVector stage2r = as<StringVector>(hstages["stage_1"]);
      StringVector stage2c = as<StringVector>(hstages["stage_2"]);
      StringVector stage1 = as<StringVector>(hstages["stage_1"]);
      int hst_size = stage3.length();
      
      arma::uvec hst_3(hst_size);
      arma::uvec hst_2r(hst_size);
      arma::uvec hst_2c(hst_size);
      arma::uvec hst_1(hst_size);
      hst_3.zeros();
      hst_2r.zeros();
      hst_2c.zeros();
      hst_1.zeros();
      
      arma::uvec di_stage32_id(di_size);
      arma::uvec di_stage21_id(di_size);
      arma::uvec di_index(di_size);
      di_stage32_id.zeros();
      di_stage21_id.zeros();
      di_index.zeros();
      
      for (int i = 0; i < di_size; i++) {
        for (int j = 0; j < hst_size; j++) {
          if (di_stage3(i) == stage3(j)) {
            hst_3(j) = 1;
          } else {
            hst_3(j) = 0;
          }
        }
        
        for (int j = 0; j < hst_size; j++) {
          if (di_stage2(i) == stage2r(j)) {
            hst_2r(j) = 1;
          } else {
            hst_2r(j) = 0;
          }
        }
        
        for (int j = 0; j < hst_size; j++) {
          if (di_stage2(i) == stage2c(j)) {
            hst_2c(j) = 1;
          } else {
            hst_2c(j) = 0;
          }
        }
        
        for (int j = 0; j < hst_size; j++) {
          if (di_stage1(i) == stage1(j)) {
            hst_1(j) = 1;
          } else {
            hst_1(j) = 0;
          }
        }
        
        arma::uvec find_hst3 = find(hst_3);
        arma::uvec find_hst2r = find(hst_2r);
        arma::uvec find_hst2c = find(hst_2c);
        arma::uvec find_hst1 = find(hst_1);
        
        arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
        arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
        
        di_stage32_id(i) = pop_32(0);
        di_stage21_id(i) = pop_21(0);
        di_index(i) = pop_32(0) + (pop_21(0) * hst_size);
        
        hst_3.zeros();
        hst_2r.zeros();
        hst_2c.zeros();
        hst_1.zeros();
      }
      
      dens_index = Rcpp::List::create(_["index32"] = di_stage32_id,
        _["index21"] = di_stage21_id, _["index321"] = di_index);
      
    } else { // Ahistorical and age-based matrices
      StringVector stage3 = as<StringVector>(ahstages["stage"]);
      StringVector stage2 = as<StringVector>(ahstages["stage"]);
      int ahst_size = stage3.length();
      
      arma::uvec ahst_3(ahst_size);
      arma::uvec ahst_2(ahst_size);
      ahst_3.zeros();
      ahst_2.zeros();
      
      arma::uvec di_stage32_id(di_size);
      arma::uvec di_stage21_id(di_size);
      arma::uvec di_index(di_size);
      di_stage32_id.zeros();
      di_stage21_id.zeros();
      di_index.zeros();
      
      for (int i = 0; i < di_size; i++) {
        for (int j = 0; j < ahst_size; j++) {
          if (di_stage3(i) == stage3(j)) {
            ahst_3(j) = 1;
          } else {
            ahst_3(j) = 0;
          }
        }
        
        for (int j = 0; j < ahst_size; j++) {
          if (di_stage2(i) == stage2(j)) {
            ahst_2(j) = 1;
          } else {
            ahst_2(j) = 0;
          }
        }
        
        arma::uvec find_ahst3 = find(ahst_3);
        arma::uvec find_ahst2 = find(ahst_2);
        di_stage32_id(i) = find_ahst3(0);
        di_stage21_id(i) = find_ahst2(0);
        di_index(i) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
        
        ahst_3.zeros();
        ahst_2.zeros();
      }
      
      dens_index = Rcpp::List::create(_["index3"] = di_stage32_id,
        _["index2"] = di_stage21_id, _["index321"] = di_index);
    }
    dyn_index321 = as<arma::uvec>(dens_index["index321"]);
    dyn_index_col = as<arma::uvec>(dens_index[1]);
    dyn_style = as<arma::uvec>(dens_input["style"]);
    dyn_alpha = as<arma::vec>(dens_input["alpha"]);
    dyn_beta = as<arma::vec>(dens_input["beta"]);
    dyn_delay = as<arma::uvec>(dens_input["time_delay"]);
    dyn_type = as<arma::uvec>(dens_input["type"]);
    n_dyn_elems = dyn_index321.n_elem;
  }
  
  // Main projection loop
  List A_all(nreps);
  List F_all(nreps);
  List U_all(nreps);
  List out_all(nreps);
  List A_mats(times);
  List F_mats(times);
  List U_mats(times);
  List out_mats(times);
  
  IntegerMatrix years_int = refsort_num(years_projected, mainyears);
  IntegerVector patches_int = refsort_str(patches_projected, mainpatches);
  
  arma::vec theseventhson = startvec;
  arma::rowvec theseventhgrandson = startvec.as_row();
  arma::mat popproj(meanmatrows, (times + 1), fill::zeros); // Population vector
  arma::mat wpopproj(meanmatrows, (times + 1), fill::zeros); // Population w vector
  arma::mat vpopproj(meanmatrows, (times + 1)); // Population v vector
  arma::mat Rvecmat(1, (times + 1), fill::zeros);
  arma::mat thesecondprophecy;
  
  popproj.col(0) = startvec;
  Rvecmat(0) = sum(startvec);
  if (!growthonly) {
    wpopproj.col(0) = startvec / sum(startvec);
  }
  
  List all_projections (nreps);
  List all_stagedist (nreps);
  List all_repvalues (nreps);
  arma::mat all_R (nreps, times+1);
  
  if (sparse_switch == 0 || format == 5) {
    for (int rep = 0; rep < nreps; rep++) {
      
      theseventhson = startvec;
      arma::rowvec theseventhgrandson = startvec.as_row();
      
      for (int i = 0; i < times; i++) {
        if (i % 50 == 0) Rcpp::checkUserInterrupt();
        
        yearnumber = years_int(i, rep) - 1;
        patchnumber = patches_int(i) - 1;
        
        used_devs = {sur_dev_values(i), obs_dev_values(i), siz_dev_values(i), sib_dev_values(i),
          sic_dev_values(i), rep_dev_values(i), fec_dev_values(i), jsur_dev_values(i),
          jobs_dev_values(i), jsiz_dev_values(i), jsib_dev_values(i), jsic_dev_values(i),
          jrep_dev_values(i), jmat_dev_values(i)};
        
        if (format < 5) {
          madsexmadrigal_oneyear = jerzeibalowski(allstages, new_stageframe,
            format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
            repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
            jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
            f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
            r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
            r2_indc_values, r1_indc_values, used_devs, spdensity_projected(i),
            repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
            yearnumber, patchnumber, exp_tol, theta_tol, ipm_method, err_check,
            true);
        } else {
          madsexmadrigal_oneyear = motherbalowski(actualages, new_stageframe,
            surv_proxy, fec_proxy, f2_inda_values, f1_inda_values, f2_indb_values,
            f1_indb_values, f2_indc_values, f1_indc_values, r2_inda_values,
            r1_inda_values, r2_indb_values, r1_indb_values, r2_indc_values,
            r1_indc_values, sur_dev_values(i), fec_dev_values(i), spdensity_projected(i),
            repmod, last_age, false, yearnumber, patchnumber, exp_tol, theta_tol, true);
        }
        
        Umat = as<arma::mat>(madsexmadrigal_oneyear["U"]);
        Fmat = as<arma::mat>(madsexmadrigal_oneyear["F"]);
        
        for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
          time_delay = dyn_delay(j);
          if (time_delay > 0) time_delay = time_delay - 1;
          
          if (i >= time_delay) {
            pop_size = sum(popproj.col(i - time_delay));
            
            if (dyn_style(j) == 1) { // Ricker
              changing_element_U = Umat(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              changing_element_F = Fmat(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              
            } else if (dyn_style(j) == 2) { // Beverton-Holt
              changing_element_U = Umat(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              changing_element_F = Fmat(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              
            } else if (dyn_style(j) == 3) { // Usher function
              changing_element_U = Umat(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              changing_element_F = Fmat(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              
            } else if (dyn_style(j) == 4) { // Logistic function
              double used_popsize = pop_size;
              if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
                used_popsize = dyn_alpha(j);
              }
              changing_element_U = Umat(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
              changing_element_F = Fmat(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
            }
            
            if (substoch == 1) {
              if (changing_element_U > 1.0 && dyn_type(j) == 1) {
                changing_element_U = 1.0;
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
            } else if (substoch == 2 && dyn_type(j) == 1) {
              double barnyard_antics {0.0};
              arma::vec given_col = Umat.col(dyn_index_col(j));
              arma::uvec gc_negs = find(given_col < 0.0);
              if (gc_negs.n_elem > 0) {
                barnyard_antics = sum(given_col(gc_negs));
              }
              changing_colsum = sum(given_col) - Umat(dyn_index321(j)) - barnyard_antics;
              
              if (changing_element_U > (1.0 - changing_colsum)) {
                changing_element_U = (1.0 - changing_colsum);
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
            } else if (substoch > 0 && dyn_type(j) == 2) {
              if (changing_element_F < 0.0) {
                changing_element_F = 0.0;
              }
            }
            Umat(dyn_index321(j)) = changing_element_U;
            Fmat(dyn_index321(j)) = changing_element_F;
            
            if (dyn_type(j) == 1 && Umat(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              Rf_warningcall(R_NilValue, "Some probabilities with value > 1.0 produced during density adjustment.");
            } else if ((Umat(dyn_index321(j)) < 0.0 || Fmat(dyn_index321(j)) < 0.0) && !warn_trigger_neg) {
              warn_trigger_neg = true;
              Rf_warningcall(R_NilValue, "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
        
        Amat = Umat + Fmat;
        
        if (err_check) {
          F_mats(i) = Fmat;
          U_mats(i) = Umat;
          A_mats(i) = Amat;
          if (format < 5) out_mats(i) = as<DataFrame>(madsexmadrigal_oneyear["out"]);
        }
        
        theseventhson = Amat * theseventhson;
        if (integeronly) {
          theseventhson = floor(theseventhson);
        }
        popproj.col(i+1) = theseventhson;
        Rvecmat(i+1) = sum(theseventhson);
        
        if (standardize) {
          theseventhson = theseventhson / sum(theseventhson);
        }
        if (!growthonly) {
          wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
          
          if (repvalue) {
            NumericVector second_devs = {sur_dev_values(times - (i+1)), obs_dev_values(times - (i+1)),
              siz_dev_values(times - (i+1)), sib_dev_values(times - (i+1)), sic_dev_values(times - (i+1)),
              rep_dev_values(times - (i+1)), fec_dev_values(times - (i+1)), jsur_dev_values(times - (i+1)),
              jobs_dev_values(times - (i+1)), jsiz_dev_values(times - (i+1)), jsib_dev_values(times - (i+1)),
              jsic_dev_values(times - (i+1)), jrep_dev_values(times - (i+1)), jmat_dev_values(times - (i+1))};
            
            if (format < 5) {
              madsexmadrigal_forward = jerzeibalowski(allstages, new_stageframe,
                format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
                repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
                jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
                f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
                r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
                r2_indc_values, r1_indc_values, used_devs, spdensity_projected(times - (i+1)),
                repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
                yearnumber, patchnumber, exp_tol, theta_tol, ipm_method, err_check,
                true);
            } else {
              madsexmadrigal_forward = motherbalowski(actualages, new_stageframe,
                surv_proxy, fec_proxy, f2_inda_values, f1_inda_values, f2_indb_values,
                f1_indb_values, f2_indc_values, f1_indc_values, r2_inda_values,
                r1_inda_values, r2_indb_values, r1_indb_values, r2_indc_values,
                r1_indc_values, sur_dev_values(times - (i+1)), fec_dev_values(times - (i+1)),
                spdensity_projected(times - (i+1)), repmod, last_age, false, yearnumber,
                patchnumber, exp_tol, theta_tol, true);
            }
            arma::mat second_U = as<arma::mat>(madsexmadrigal_forward["U"]);
            arma::mat second_F = as<arma::mat>(madsexmadrigal_forward["F"]);
            thesecondprophecy = second_U + second_F;
            theseventhgrandson = theseventhgrandson * thesecondprophecy;
            
            double seventhgrandsum = sum(theseventhgrandson);
            arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
            theseventhgrandson = theseventhgrandson / seventhgrandsum;
            vpopproj.col(times - (i+1)) = midwife;
          } // if repvalue
        } // if growthonly
      } // times for loop
      
      all_projections(rep) = popproj;
      all_R.row(rep) = Rvecmat;
      
      if (!growthonly) {
        all_stagedist(rep) = wpopproj;
        if (repvalue) all_repvalues(rep) = vpopproj;
      }
      
      if (err_check) {
        A_all(rep) = A_mats;
        F_all(rep) = F_mats;
        U_all(rep) = U_mats;
        if (format < 5) out_all(rep) = out_mats;
      }
    } // nreps for loop
  } else {
    arma::sp_mat Umat_sp;
    arma::sp_mat Fmat_sp;
    arma::sp_mat Amat_sp;
    arma::sp_mat thesecondprophecy_sp;
    arma::sp_mat theseventhson_sp(theseventhson);
    arma::sp_mat theseventhgrandson_sp(theseventhgrandson);
    
    for (int rep = 0; rep < nreps; rep++) {
      for (int i = 0; i < times; i++) {
        if (i % 50 == 0) Rcpp::checkUserInterrupt();
        
        yearnumber = years_int(i, rep) - 1;
        patchnumber = patches_int(i) - 1;
        
        used_devs = {sur_dev_values(i), obs_dev_values(i), siz_dev_values(i), sib_dev_values(i),
          sic_dev_values(i), rep_dev_values(i), fec_dev_values(i), jsur_dev_values(i),
          jobs_dev_values(i), jsiz_dev_values(i), jsib_dev_values(i), jsic_dev_values(i),
          jrep_dev_values(i), jmat_dev_values(i)};
        
        madsexmadrigal_oneyear = jerzeibalowski_sp(allstages, new_stageframe,
          format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
          repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
          jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
          f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
          r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
          r2_indc_values, r1_indc_values, used_devs, spdensity_projected(i),
          repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
          yearnumber, patchnumber, exp_tol, theta_tol, ipm_method, err_check,
          true);
        
        Umat_sp = as<arma::sp_mat>(madsexmadrigal_oneyear["U"]);
        Fmat_sp = as<arma::sp_mat>(madsexmadrigal_oneyear["F"]);
        
        for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
          time_delay = dyn_delay(j);
          if (time_delay > 0) time_delay = time_delay - 1;
          
          if (i >= time_delay) {
            pop_size = sum(popproj.col(i - time_delay));
            
            if (dyn_style(j) == 1) { // Ricker
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
              
            } else if (dyn_style(j) == 2) { // Beverton-Holt
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
              
            } else if (dyn_style(j) == 3) { // Usher function
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
              
            } else if (dyn_style(j) == 4) { // Logistic function
              double used_popsize = pop_size;
              if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
                used_popsize = dyn_alpha(j);
              }
              changing_element_U = Umat_sp(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
              changing_element_F = Fmat_sp(dyn_index321(j)) * 
                (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
            }
            
            if (substoch == 1) {
              if (changing_element_U > 1.0 && dyn_type(j) == 1) {
                changing_element_U = 1.0;
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
              
            } else if (substoch == 2 && dyn_type(j) == 1) {
              double barnyard_antics {0.0};
              arma::vec given_col = arma::vec(Umat_sp.col(dyn_index_col(j)));
              arma::uvec gc_negs = find(given_col < 0.0);
              if (gc_negs.n_elem > 0) {
                barnyard_antics = sum(given_col(gc_negs));
              }
              changing_colsum = sum(given_col) - Umat_sp(dyn_index321(j)) - barnyard_antics;
              
              if (changing_element_U > (1.0 - changing_colsum)) {
                changing_element_U = (1.0 - changing_colsum);
              } else if (changing_element_U < 0.0) {
                changing_element_U = 0.0;
              }
              
            } else if (substoch > 0 && dyn_type(j) == 2) {
              if (changing_element_F < 0.0) {
                changing_element_F = 0.0;
              }
            }
            Umat_sp(dyn_index321(j)) = changing_element_U;
            Fmat_sp(dyn_index321(j)) = changing_element_F;
            
            if (dyn_type(j) == 1 && Umat_sp(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              Rf_warningcall(R_NilValue, "Some probabilities with value > 1.0 produced during density adjustment.");
            } else if ((Umat_sp(dyn_index321(j)) < 0.0 || Fmat_sp(dyn_index321(j)) < 0.0) &&
              !warn_trigger_neg) {
              warn_trigger_neg = true;
              Rf_warningcall(R_NilValue, "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
        
        Amat_sp = Umat_sp + Fmat_sp;
        
        if (err_check) {
          F_mats(i) = Fmat_sp;
          U_mats(i) = Umat_sp;
          A_mats(i) = Amat_sp;
          if (format < 5) out_mats(i) = as<DataFrame>(madsexmadrigal_oneyear["out"]);
        }
        
        theseventhson_sp = Amat_sp * theseventhson_sp;
        if (integeronly) {
          theseventhson_sp = floor(theseventhson_sp);
        }
        popproj.col(i+1) = arma::vec(theseventhson_sp);
        Rvecmat(i+1) = sum(popproj.col(i+1));
        
        if (standardize) {
          theseventhson_sp = theseventhson_sp / Rvecmat(i+1);
        }
        if (!growthonly) {
          wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
          
          if (repvalue) {
            NumericVector second_devs = {sur_dev_values(times - (i+1)), obs_dev_values(times - (i+1)),
              siz_dev_values(times - (i+1)), sib_dev_values(times - (i+1)), sic_dev_values(times - (i+1)),
              rep_dev_values(times - (i+1)), fec_dev_values(times - (i+1)), jsur_dev_values(times - (i+1)),
              jobs_dev_values(times - (i+1)), jsiz_dev_values(times - (i+1)), jsib_dev_values(times - (i+1)),
              jsic_dev_values(times - (i+1)), jrep_dev_values(times - (i+1)), jmat_dev_values(times - (i+1))};
            
            madsexmadrigal_forward = jerzeibalowski_sp(allstages, new_stageframe,
              format, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy,
              repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy,
              jsizec_proxy, jrepst_proxy, jmatst_proxy, f2_inda_values, f1_inda_values,
              f2_indb_values, f1_indb_values, f2_indc_values, f1_indc_values,
              r2_inda_values, r1_inda_values, r2_indb_values, r1_indb_values,
              r2_indc_values, r1_indc_values, used_devs, spdensity_projected(times - (i+1)),
              repmod, maxsize, maxsizeb, maxsizec, start_age, last_age, false,
              yearnumber, patchnumber, exp_tol, theta_tol, ipm_method, err_check,
              true);
            arma::sp_mat second_U = as<arma::sp_mat>(madsexmadrigal_forward["U"]);
            arma::sp_mat second_F = as<arma::sp_mat>(madsexmadrigal_forward["F"]);
            thesecondprophecy_sp = second_U + second_F;
            theseventhgrandson_sp = theseventhgrandson_sp * thesecondprophecy_sp;
            
            arma::rowvec partway(theseventhgrandson_sp);
            double seventhgrandsum = sum(partway);
            arma::sp_mat midwife = theseventhgrandson_sp / seventhgrandsum;
            theseventhgrandson_sp = midwife;
            vpopproj.col(times - (i+1)) = midwife;
          } // if repvalue
        } // if growthonly
      } // times for loop
      
      all_projections(rep) = popproj;
      all_R.row(rep) = Rvecmat;
      
      if (!growthonly) {
        all_stagedist(rep) = wpopproj;
        if (repvalue) all_repvalues(rep) = vpopproj;
      }
      
      if (err_check) {
        A_all(rep) = A_mats;
        F_all(rep) = F_mats;
        U_all(rep) = U_mats;
        if (format < 5) out_all(rep) = out_mats;
      }
    } // nreps for loop
  }
  
  // Final output prep
  IntegerVector control = {nreps, times};
  DataFrame newlabels = DataFrame::create(_["pop"] = 1, _["patch"] = chosenpatch);
  CharacterVector output_class = {"lefkoProj"};
  List output;
  
  List all_R_list(1);
  all_R_list(0) = all_R;
  
  if (err_check) {
    List output_err(37);
    
    output_err(0) = all_projections;
    output_err(1) = all_stagedist;
    output_err(2) = all_repvalues;
    output_err(3) = all_R_list;
    output_err(4) = ahstages;
    output_err(5) = hstages;
    output_err(6) = agestages;
    output_err(7) = newlabels;
    output_err(8) = control;
    output_err(9) = dens_input;
    output_err(10) = pmnames;
    output_err(11) = mainyears;
    output_err(12) = mainpatches;
    output_err(13) = maingroups;
    output_err(14) = mainages;
    output_err(15) = allstages;
    output_err(16) = surv_proxy;
    output_err(17) = obs_proxy;
    output_err(18) = size_proxy;
    output_err(19) = sizeb_proxy;
    output_err(20) = sizec_proxy;
    output_err(21) = repst_proxy;
    output_err(22) = fec_proxy;
    output_err(23) = jsurv_proxy;
    output_err(24) = jobs_proxy;
    output_err(25) = jsize_proxy;
    output_err(26) = jsizeb_proxy;
    output_err(27) = jsizec_proxy;
    output_err(28) = jrepst_proxy;
    output_err(29) = jmatst_proxy;
    output_err(30) = years_projected;
    output_err(31) = patches_projected;
    output_err(32) = spdensity_projected;
    output_err(33) = A_all;
    output_err(34) = U_all;
    output_err(35) = F_all;
    output_err(36) = out_all;
    
    CharacterVector output_err_names = {"projection", "stage_dist", "rep_value",
      "pop_size", "ahstages", "hstages", "agestages", "labels", "control", "density",
      "paramnames", "mainyears", "mainpatches", "maingroups", "mainages",
      "allstages", "surv_proxy", "obs_proxy", "size_proxy", "sizeb_proxy",
      "sizec_proxy", "repst_proxy", "fec_proxy", "jsurv_proxy", "jobs_proxy",
      "jsize_proxy", "jsizeb_proxy", "jsizec_proxy", "jrepst_proxy", "jmatst_proxy",
      "years_projected", "patches_projected", "spdensity_projected", "A_all",
      "U_all", "F_all", "out_all"};
    output_err.attr("names") = output_err_names;
    output_err.attr("class") = output_class;
    
    output = output_err;
  } else {
    List output_noerr(10);
    
    output_noerr(0) = all_projections;
    output_noerr(1) = all_stagedist;
    output_noerr(2) = all_repvalues;
    output_noerr(3) = all_R_list;
    output_noerr(4) = ahstages;
    output_noerr(5) = hstages;
    output_noerr(6) = agestages;
    output_noerr(7) = newlabels;
    output_noerr(8) = control;
    output_noerr(9) = dens_input;
    
    CharacterVector output_noerr_names = {"projection", "stage_dist", "rep_value",
      "pop_size", "ahstages", "hstages", "agestages", "labels", "control", "density"};
    output_noerr.attr("names") = output_noerr_names;
    output_noerr.attr("class") = output_class;
    
    output = output_noerr;
  }
  
  return output;
}

// pop dynamics

//' Vectorize Matrix for Historical Mean Matrix Estimation
//' 
//' Function \code{flagrantcrap()} vectorizes core indices of matrices
//' input as list elements.
//' 
//' @name flagrantcrap
//' 
//' @param Xmat A matrix originally a part of a list object.
//' @param allindices A vector of indices to remove from the matrix
//' 
//' @return A column vector of specifically called elements from the input
//' matrix.
//' 
//' @keywords internal
//' @noRd
arma::vec flagrantcrap(arma::mat Xmat, arma::uvec allindices) {
  
  arma::vec newcol = Xmat.elem(allindices);
  
  return newcol;
}

//' Vectorize Matrix for Ahistorical Mean Matrix Estimation
//' 
//' Function \code{moreflagrantcrap()} vectorizes matrices input as list
//' elements.
//' 
//' @name moreflagrantcrap
//' 
//' @param Xmat A matrix originally a part of a list object.
//' 
//' @return A column vector of the input matrix.
//' 
//' @keywords internal
//' @noRd
arma::vec moreflagrantcrap(arma::mat Xmat) {
  
  arma::vec newcol = arma::vectorise(Xmat);
  
  return newcol;
}

//' Calculate Logarithms of Non-Zero Elements of Sparse Matrix
//' 
//' Function \code{spmat_log} finds the non-zero elements in a sparse matrix,
//' calculates their logs, and inserts them back into the matrix and returns it.
//' Based on code developed by Coatless Professor and posted by him on
//' StackOverflow.
//' 
//' @name spmat_log
//' 
//' @param B A sparse matrix. Note that this is assumed to be a population
//' projection matrix, meaning that all values are either 0 or positive.
//' 
//' @return A sparse matrix with non-zero values as logs of the elements in the
//' input matrix.
//' 
//' @keywords internal
//' @noRd
arma::sp_mat spmat_log(arma::sp_mat coremat)
{
  arma::sp_mat::const_iterator start = coremat.begin();
  arma::sp_mat::const_iterator end   = coremat.end();
  arma::sp_mat::const_iterator it = start; 
  
  int n = std::distance(start, end);
  
  if (n > 0) {
    arma::umat locs(2, n);
    arma::uvec temp(2);
    arma::vec vals(n);
    arma::vec logvals(n);
    locs.zeros();
    temp.zeros();
    vals.zeros();
    logvals.zeros();
    
    for(int i = 0; i < n; ++i) {
      temp(0) = it.row();
      temp(1) = it.col();
      locs.col(i) = temp;
      vals(i) = coremat(temp(0), temp(1));
      logvals(i) = log(vals(i));
      ++it; // increment
    }
    
    coremat = arma::sp_mat(locs, logvals, coremat.n_rows, coremat.n_cols);
  }
  
  return coremat;
}

//' Estimates Mean LefkoMat Object for Historical MPM
//' 
//' Function \code{turbogeodiesel()} estimates mean historical population
//' projection matrices, treating the mean as element-wise arithmetic.
//' 
//' @name turbogeodiesel
//' 
//' @param loy A data frame denoting the population, patch, and occasion
//' designation for each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' @param agestages This is the \code{agestages} object held by \code{mats}.
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param patchmats A logical value stating whether to estimate patch-level
//' means.
//' @param popmats A logical value stating whether to estimate population-level
//' means.
//' 
//' @return A list using the structure of a lefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.turbogeodiesel)]]
List turbogeodiesel(DataFrame loy, List Umats, List Fmats, DataFrame hstages, 
  DataFrame agestages, DataFrame stages, bool patchmats, bool popmats) {
  
  StringVector pops = as<StringVector>(loy["pop"]);
  arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
  StringVector patches = as<StringVector>(loy["patch"]);
  arma::uvec year2 = as<arma::uvec>(loy["year2"]);
  arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
  arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
  arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  int loydim = pops.length();
  int numofpops = uniquepops.n_elem;
  int numofpatches = uniquepoppatches.n_elem;
  
  if (numofpatches == 1) popmats = 0;
  
  StringVector poporderlong(loydim);
  arma::uvec poporderlong_num(loydim);
  StringVector patchorderlong(loydim);
  arma::uvec annmatriceslong(loydim);
  arma::uvec meanassign(loydim);
  poporderlong_num.zeros();
  annmatriceslong.zeros();
  meanassign.zeros();
  
  pop_num = pop_num + 1;
  poppatchc = poppatchc + 1;
  
  poporderlong(0) = pops(0);
  poporderlong_num(0) = pop_num(0);
  patchorderlong(0) = patches(0);
  annmatriceslong(0) = 1;
  meanassign(0) = 1;
  
  int counter {0};
  
  StringVector uniquepops_str(numofpops);
  uniquepops_str(0) = pops(0);
  int popcounter {0};
  
  // Here we assess how many mean matrices we need, and their overall order
  if (loydim > 1) {
    for (int i = 1; i < loydim; i++) {
      if (poppatchc(i) != poppatchc(i-1)) {
        counter++;
        poporderlong(counter) = pops(i);
        poporderlong_num(counter) = pop_num(i);
        patchorderlong(counter) = patches(i);
        annmatriceslong(counter) = 1;
        meanassign(i) = meanassign(i-1) + 1;
        
        if (pop_num(i) != pop_num(i-1)) {
          popcounter += 1;
          uniquepops_str(popcounter) = pops(i);
        }
        
      } else {
        annmatriceslong(counter) = annmatriceslong(counter) + 1;
        meanassign(i) = meanassign(i-1);
      }
    }
  }
  
  arma::uvec toestimate = find(poporderlong_num);
  int popcount = toestimate.n_elem;
  
  int totalmatrices = toestimate.n_elem + numofpops;
  
  if (patchmats == 1 && popmats == 0) {
    totalmatrices = toestimate.n_elem;
  } else if (patchmats == 0 && popmats == 1) {
    totalmatrices = numofpops;
  }
  
  arma::uvec poporder = poporderlong_num.elem(toestimate);
  arma::uvec patchorder = poppatchc.elem(toestimate);
  arma::uvec annmatrices = annmatriceslong.elem(toestimate);
  
  StringVector poporder_str(popcount);
  StringVector patchorder_str(popcount);
  
  for (int i = 0; i < popcount; i++) {
    poporder_str(i) = pops(toestimate(i));
    patchorder_str(i) = patches(toestimate(i));
  }
  
  // This next chunk predicts which elements will be targeted for arithmetic mean estimation
  int format_int {0};
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  StringVector stagenames = as<StringVector>(stages["stage"]);
  int numstages = astages.n_elem;
  
  if (stagenames(numstages - 1) == "AlmostBorn") format_int = 1;
  
  arma::uvec hstage3in = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec hstage2nin = as<arma::uvec>(hstages["stage_id_1"]);
  int numhstages = hstage3in.n_elem;
  
  int predictedsize = 2 * numstages * numstages * numstages;
  
  arma::uvec hsindexl(predictedsize);
  hsindexl.zeros();
  
  counter = 0;
  
  if (format_int == 0) { // Ehrlen format
    for (int i1 = 0; i1 < numhstages; i1++) {
      for (int i2 = 0; i2 < numhstages; i2++) {
        if (hstage3in(i1) == hstage2nin(i2)) {
          hsindexl(counter) = (i1 * numhstages) + i2;
          counter++;
        }
      }
    }
  } else { // deVries format
    for (int i1 = 0; i1 < numhstages; i1++) {
      for (int i2 = 0; i2 < numhstages; i2++) {
        if (hstage3in(i1) == hstage2nin(i2)) {
          hsindexl(counter) = (i1 * numhstages) + i2;
          counter++;
          
        } else if (hstage2nin(i2) == numstages || hstage3in(i1) == numstages) {
          hsindexl(counter) = (i1 * numhstages) + i2;
          counter++;
        }
      }
    }
  }
  
  arma::uvec hsgood = find(hsindexl);
  arma::uvec hsindex = hsindexl.elem(hsgood);
  arma::uvec zerovec(1);
  zerovec.zeros();
  arma::uvec allindices = join_cols(zerovec, hsindex);
  
  // Now we build U and F matrices of element-wise arithmetic means, where
  // each column corresponds to the predicted non-zero elements of each mean
  // matrix, and each matrix is presented as a column vector within the 
  // overall matrix. The A matrix is the sum of U and F.
  int core_elem = counter;
  
  arma::mat umatvec(core_elem, totalmatrices);
  arma::mat fmatvec(core_elem, totalmatrices);
  umatvec.zeros();
  fmatvec.zeros();
  
  int patchchoice {0};
  int popchoice {0};
  
  pop_num = pop_num - 1;
  poppatchc = poppatchc - 1;
  
  for (int i = 0; i < loydim; i++) {
    if (patchmats == 1) {
      patchchoice = poppatchc(i);
      
      umatvec.col(patchchoice) = umatvec.col(patchchoice) +
        (flagrantcrap(as<arma::mat>(Umats[i]), allindices) / yearsinpatch(i));
      fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
        (flagrantcrap(as<arma::mat>(Fmats[i]), allindices) / yearsinpatch(i));
    }
    
    if (popmats == 1) {
      if (patchmats == 1) {
        popchoice = numofpatches + pop_num(i);
      } else {
        popchoice = pop_num(i);
      }
      
      umatvec.col(popchoice) = umatvec.col(popchoice) +
        (flagrantcrap(as<arma::mat>(Umats[i]), allindices) / (yearsinpatch(i) * patchesinpop(i)));
      fmatvec.col(popchoice) = fmatvec.col(popchoice) +
        (flagrantcrap(as<arma::mat>(Fmats[i]), allindices) / (yearsinpatch(i) * patchesinpop(i)));
    }
  }
  arma::mat amatvec = umatvec + fmatvec;
  
  // Here we create the cheat sheet algorithm
  int cheatsheetlength {1};
  if (numofpatches > 1) cheatsheetlength = numofpops + numofpatches;
  StringVector poporder_redone(cheatsheetlength);
  StringVector patchorder_redone(cheatsheetlength);
  
  if (numofpatches > 1) {
    for (int i = 0; i < numofpatches; i++) {
      poporder_redone(i) = poporderlong(i);
      patchorder_redone(i) = patchorderlong(i);
    }
    
    for (int i = 0; i < numofpops; i++) {
      poporder_redone(i+numofpatches) = uniquepops_str(i);
      patchorder_redone(i+numofpatches) = "0";
    }
    
  } else {
    poporder_redone(0) = poporderlong(0);
    patchorder_redone(0) = patchorderlong(0);
  }
  DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone,
    _["patch"] = patchorder_redone);
  
  // Now we will create the main list objects holding the matrices
  List U(totalmatrices);
  List F(totalmatrices);
  List A(totalmatrices);
  
  arma::mat umat_base(numhstages, numhstages);
  arma::mat fmat_base(numhstages, numhstages);
  arma::mat amat_base(numhstages, numhstages);
  
  for (int i = 0; i < totalmatrices; i++) {
    umat_base.zeros();
    fmat_base.zeros();
    amat_base.zeros();
    
    umat_base.elem(allindices) = umatvec.col(i);
    fmat_base.elem(allindices) = fmatvec.col(i);
    amat_base.elem(allindices) = amatvec.col(i);
    
    U(i) = umat_base;
    F(i) = fmat_base;
    A(i) = amat_base;
  }

  // Matrix QC output
  arma::uvec utrans = find(umatvec);
  arma::uvec ftrans = find(fmatvec);
  int totalutrans = utrans.n_elem;
  int totalftrans = ftrans.n_elem;
  
  NumericVector matrixqc(3);
  matrixqc(0) = totalutrans; // summed number of non-zero u transitions
  matrixqc(1) = totalftrans; // summed number of non-zero f transitions
  matrixqc(2) = totalmatrices;
  
  // Final output
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F,
    _["hstages"] = hstages, _["agestages"] = agestages, _["ahstages"] = stages,
    _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
  
  return output;
}

//' Estimates Mean LefkoMat Object for Ahistorical MPM
//' 
//' Function \code{geodiesel()} estimates mean ahistorical population
//' projection matrices, treating the mean as element-wise arithmetic. The
//' function can handle both normal ahistorical MPMs and age x stage ahistorical
//' MPMs.
//' 
//' @name geodiesel
//' 
//' @param loy A data frame denoting the population, patch, and occasion
//' designation of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param agestages This is the \code{agestages} object held by \code{mats}.
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param patchmats A logical value stating whether to estimate patch-level
//' means.
//' @param popmats A logical value stating whether to estimate population-level
//' means.
//' 
//' @return A list using the structure of a LefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.geodiesel)]]
List geodiesel(DataFrame loy, List Umats, List Fmats, DataFrame agestages,
  DataFrame stages, bool patchmats, bool popmats) {
  
  StringVector pops = as<StringVector>(loy["pop"]);
  arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
  StringVector patches = as<StringVector>(loy["patch"]);
  arma::uvec year2 = as<arma::uvec>(loy["year2"]);
  arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
  arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
  arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  int loydim = pops.length();
  int numofpops = uniquepops.n_elem;
  int numofpatches = uniquepoppatches.n_elem;
  
  if (numofpatches == 1) popmats = 0;
  
  StringVector poporderlong(loydim);
  arma::uvec poporderlong_num(loydim);
  StringVector patchorderlong(loydim);
  arma::uvec annmatriceslong(loydim);
  arma::uvec meanassign(loydim);
  poporderlong_num.zeros();
  annmatriceslong.zeros();
  meanassign.zeros();
  
  pop_num = pop_num + 1;
  poppatchc = poppatchc + 1;
  
  poporderlong(0) = pops(0);
  poporderlong_num(0) = pop_num(0);
  patchorderlong(0) = patches(0);
  annmatriceslong(0) = 1;
  meanassign(0) = 1;
  
  int counter {0};
  
  StringVector uniquepops_str(numofpops);
  uniquepops_str(0) = pops(0);
  int popcounter {0};
  
  // Here we assess how many mean matrices we need, and their overall order
  if (loydim > 1) {
    for (int i = 1; i < loydim; i++) {
      if (poppatchc(i) != poppatchc(i-1)) {
        counter++;
        poporderlong(counter) = pops(i);
        poporderlong_num(counter) = pop_num(i);
        patchorderlong(counter) = patches(i);
        annmatriceslong(counter) = 1;
        meanassign(i) = meanassign(i-1) + 1;
        
        if (pop_num(i) != pop_num(i-1)) {
          popcounter += 1;
          uniquepops_str(popcounter) = pops(i);
        }
        
      } else {
        annmatriceslong(counter) = annmatriceslong(counter) + 1;
        meanassign(i) = meanassign(i-1);
      }
    }
  }
  arma::uvec toestimate = find(poporderlong_num);
  int popcount = toestimate.n_elem;
  
  int totalmatrices = toestimate.n_elem + numofpops;
  
  if (patchmats == 1 && popmats == 0) {
    totalmatrices = toestimate.n_elem;
  } else if (patchmats == 0 && popmats == 1) {
    totalmatrices = numofpops;
  }
  
  arma::uvec poporder = poporderlong_num.elem(toestimate);
  arma::uvec patchorder = poppatchc.elem(toestimate);
  arma::uvec annmatrices = annmatriceslong.elem(toestimate);
  
  StringVector poporder_str(popcount);
  StringVector patchorder_str(popcount);
  
  for (int i = 0; i < popcount; i++) {
    poporder_str(i) = pops(toestimate(i));
    patchorder_str(i) = patches(toestimate(i));
  }
  
  // This next chunk predicts which elements will be targeted for arithmetic mean estimation
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  int initialstages = astages.n_elem;
  
  // Now we will tet for the presence of ages, and determine the matrix dimensions required
  arma::mat initUmat = Umats(0);
  int colsused = initUmat.n_cols;
  int agemultiplier = colsused / initialstages;
  
  int numstages = astages.n_elem * agemultiplier;
  
  // Now we build U and F matrices of element-wise arithmetic means, where
  // each column corresponds to the predicted non-zero elements of each mean
  // matrix, and each matrix is presented as a column vector within the 
  // overall matrix. The A matrix is the sum of U and F.
  int core_elem = numstages * numstages;
  
  arma::mat umatvec(core_elem, totalmatrices);
  arma::mat fmatvec(core_elem, totalmatrices);
  umatvec.zeros();
  fmatvec.zeros();
  
  int patchchoice {0};
  int popchoice {0};
  
  pop_num = pop_num - 1;
  poppatchc = poppatchc - 1;
  
  for (int i = 0; i < loydim; i ++) {
    if (patchmats == 1) {
      patchchoice = poppatchc(i);
      
      umatvec.col(patchchoice) = umatvec.col(patchchoice) +
        (moreflagrantcrap(as<arma::mat>(Umats[i])) / yearsinpatch(i));
      fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
        (moreflagrantcrap(as<arma::mat>(Fmats[i])) / yearsinpatch(i));
    }
    
    if (popmats == 1) {
      if (patchmats == 1) {
        popchoice = numofpatches + pop_num(i);
      } else {
        popchoice = pop_num(i);
      }
      
      umatvec.col(popchoice) = umatvec.col(popchoice) +
        (moreflagrantcrap(as<arma::mat>(Umats[i])) / (yearsinpatch(i) * patchesinpop(i)));
      fmatvec.col(popchoice) = fmatvec.col(popchoice) +
        (moreflagrantcrap(as<arma::mat>(Fmats[i])) / (yearsinpatch(i) * patchesinpop(i)));
    }
  }
  arma::mat amatvec = umatvec + fmatvec;
  
  // Here we create the cheat sheet algorithm
  int cheatsheetlength {1};
  if (numofpatches > 1) cheatsheetlength = numofpops + numofpatches;
  StringVector poporder_redone(cheatsheetlength);
  StringVector patchorder_redone(cheatsheetlength);
  
  if (numofpatches > 1) {
    for (int i = 0; i < numofpatches; i++) {
      poporder_redone(i) = poporderlong(i);
      patchorder_redone(i) = patchorderlong(i);
    }
    
    for (int i = 0; i < numofpops; i++) {
      poporder_redone(i+numofpatches) = uniquepops_str(i);
      patchorder_redone(i+numofpatches) = "0";
    }
    
  } else {
    poporder_redone(0) = poporderlong(0);
    patchorder_redone(0) = patchorderlong(0);
  }
  
  DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, 
    _["patch"] = patchorder_redone);
  
  // Now we will create the main list objects to hold the matrices
  List U(totalmatrices);
  List F(totalmatrices);
  List A(totalmatrices);
  
  arma::mat umat_base = umatvec.col(0);
  arma::mat fmat_base = fmatvec.col(0);
  arma::mat amat_base = amatvec.col(0);
  
  umat_base.reshape(numstages, numstages);
  fmat_base.reshape(numstages, numstages);
  amat_base.reshape(numstages, numstages);
  
  U(0) = umat_base;
  F(0) = fmat_base;
  A(0) = amat_base;
  
  if (totalmatrices > 1) {
    for (int i = 1; i < totalmatrices; i++) {
      umat_base.zeros();
      fmat_base.zeros();
      amat_base.zeros();
      
      umat_base = umatvec.col(i);
      fmat_base = fmatvec.col(i);
      amat_base = amatvec.col(i);
      
      umat_base.reshape(numstages, numstages);
      fmat_base.reshape(numstages, numstages);
      amat_base.reshape(numstages, numstages);
      
      U(i) = umat_base;
      F(i) = fmat_base;
      A(i) = amat_base;
    }
  }
  
  // Matrix QC output
  arma::uvec utrans = find(umatvec);
  arma::uvec ftrans = find(fmatvec);
  int totalutrans = utrans.n_elem;
  int totalftrans = ftrans.n_elem;
  
  NumericVector matrixqc(3);
  matrixqc(0) = totalutrans; // summed number of U transitions
  matrixqc(1) = totalftrans; // summed number of F transitions
  matrixqc(2) = totalmatrices;
  
  // Final output
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F, 
    _["hstages"] = R_NilValue, _["agestages"] = agestages, _["ahstages"] = stages,
    _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
  
  return output;
}

//' Full Eigen Analysis of a Single Dense Matrix
//' 
//' Function \code{decomp3()} returns all eigenvalues, right eigenvectors, and
//' left eigenvectors estimated for a matrix by the \code{eig_gen}() function
//' in the C++ Armadillo library. Works with dense matrices.
//' 
//' @name decomp3
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.decomp3)]]
List decomp3(arma::mat Amat) {
  arma::cx_vec Aeigval;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
  
  List output = List::create(Named("eigenvalues") = Aeigval,
    _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Full Eigen Analysis of a Single Sparse Matrix
//' 
//' Function \code{decomp3sp()} returns all eigenvalues, right eigenvectors, and
//' left eigenvectors estimated for a matrix by the \code{eigs_gen}() function
//' in the C++ Armadillo library. Works with sparse matrices.
//' 
//' @name decomp3sp
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.decomp3sp)]]
List decomp3sp(arma::mat Amat) {
  arma::sp_mat spAmat(Amat);
  arma::sp_mat t_spAmat = spAmat.t();
  arma::cx_vec Aeigval;
  arma::cx_vec Aeigvall;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eigs_gen(Aeigval, Aeigvecr, spAmat, 1);
  eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
  
  List output = List::create(Named("eigenvalues") = Aeigval,
    _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Full Eigen Analysis of a Single Sparse Matrix, with Sparse Input
//' 
//' \code{decomp3sp_inp()} returns all eigenvalues, right eigenvectors, and left
//' eigenvectors estimated for a matrix by the \code{eigs_gen}() function
//' in the C++ Armadillo library. Works with sparse matrices.
//' 
//' @name decomp3sp_inp
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @section Notes:
//' This function works slightly differently from function \code{decomp3sp()} in
//' that the latter function requires a sparse matrix provided in dense format,
//' while this function requires a sparse matrix in sparse format.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.decomp3sp_inp)]]
List decomp3sp_inp(arma::sp_mat spAmat) {
  arma::sp_mat t_spAmat = spAmat.t();
  arma::cx_vec Aeigval;
  arma::cx_vec Aeigvall;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eigs_gen(Aeigval, Aeigvecr, spAmat, 1);
  eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
  
  List output = List::create(Named("eigenvalues") = Aeigval,
    _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Estimate Deterministic Population Growth Rate of Any Matrix
//' 
//' \code{lambda3matrix()} returns the dominant eigenvalue of a single
//' dense or sparse projection matrix, provided in dense matrix format.
//' 
//' @name lambda3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse matrix
//' format.
//'
//' @return This function returns the dominant eigenvalue of the matrix. This
//' is given as the largest real part of all eigenvalues estimated via the 
//' \code{eig_gen}() and \code{eigs_gen}() functions in the C++ Armadillo
//' library.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.lambda3matrix)]]
double lambda3matrix(arma::mat Amat, bool sparse) {
  double lambda {0};
  
  if (!sparse) {
    List eigenstuff = decomp3(Amat);
    arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
    
    lambda = max(realeigenvals);
    
  } else {
    List eigenstuff = decomp3sp(Amat);
    arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
    
    lambda = max(realeigenvals);
  }
  
  return lambda;
}

//' Estimate Stable Stage Distribution of Any Population Matrix
//' 
//' \code{ss3matrix()} returns the stable stage distribution for a 
//' dense or sparse population matrix.
//' 
//' @name ss3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix.
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ss3matrix)]]
arma::vec ss3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.0000000001); // This line replaces all numbers lower than 1 x 10-10 with 0
  
  double rvsum = sum(realrightvec);
  realrightvec = realrightvec / rvsum;
  
  return realrightvec;
}

//' Estimate Reproductive Value of Any Population Matrix
//' 
//' \code{rv3matrix()} returns the reproductive values for stages in a
//' dense or sparse population matrix (both provided in dense matrix format).
//' The function provides standard reproductive values, meaning that the overall
//' reproductive values of basic life history stages in a historical matrix are
//' not provided (the \code{\link{repvalue3.lefkoMat}()} function estimates
//' these on the basis of stage description information provided in the
//' \code{lefkoMat} object used as input in that function).
//' 
//' @name rv3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a vector characterizing the reproductive
//' values for stages of a population projection matrix.
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.rv3matrix)]]
arma::vec rv3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.0000000001); // This line replaces all numbers lower than 1 x 10-10 with 0

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  
  realleftvec = realleftvec / rlvmin;
  
  return realleftvec;
}

//' Estimate Deterministic Sensitivities of Any Population Matrix
//' 
//' \code{sens3matrix()} returns the sensitivity of lambda with respect
//' to each element in a dense or sparse matrix (provided in dense matrix
//' format). This is accomplished via the \code{eig_gen}() and \code{eigs_gen}()
//' functions in the C++ Armadillo library.
//' 
//' @name sens3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a matrix of deterministic sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3matrix)]]
arma::mat sens3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;
  
  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of A Spars Matrixe
//' 
//' \code{sens3sp_matrix()} returns the sensitivity of lambda with respect
//' to each element in a sparse matrix, provided in sparse matrix format. This
//' is accomplished via the \code{eigs_gen}() function in the C++ Armadillo
//' library.
//' 
//' @name sens3sp_matrix
//' 
//' @param Aspmat A population projection matrix in sparse matrix format.
//' @param refmat A sparse matrix used for reference to create associated 0s in
//' the sensitivity matrix.
//' 
//' @return This function returns a sparse matrix of deterministic
//' sensitivities. Zeroes are derived from the reference matrix, and replace
//' non-zero entries that will be zeroed out in the following math. Currently
//' used in LTRE estimation.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3sp_matrix)]]
arma::sp_mat sens3sp_matrix(arma::sp_mat Aspmat, arma::sp_mat refmat) {
  List eigenstuff = decomp3sp_inp(Aspmat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;
  
  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel);
  arma::sp_mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      if (refmat(i, j) != 0) {
        smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
      }
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of a Historical LefkoMat Object
//' 
//' \code{sens3hlefko()} returns the sensitivity of lambda with respect
//' to each historical stage-pair in the matrix, and the associated
//' sensitivity for each life history stage. This is accomplished via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @name sens3hlefko
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic sensitivity
//' matrices:
//' \item{h_smat}{Matrix of sensitivities corresponding to the historical
//' matrix.}
//' \item{ah_smat}{Matrix of sensitivities corresponding to the ahistorical
//' matrix.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3hlefko)]]
List sens3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  arma::uvec stage_id = as<arma::uvec>(ahstages["stage_id"]);
  arma::uvec h_stage_2 = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec h_stage_1 = as<arma::uvec>(hstages["stage_id_1"]);
  
  List eigenstuff = decomp3sp(Amat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // Using a lower threshold than previously here
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;
  
  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // Using a lower threshold than previously here

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  int ahstagelength = stage_id.n_elem;
  
  arma::vec wcorrah (ahstagelength);
  arma::vec vcorrah (ahstagelength);
  arma::vec vwprodah (ahstagelength);
  arma::mat ahsens(ahstagelength, ahstagelength);
  wcorrah.zeros();
  vcorrah.zeros();
  vwprodah.zeros();
  ahsens.zeros();
  
  // Here we create the scalar product vw and the ahistorical stable stage distribution w
  int ahrows {0};
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
    ahrows = h_stage_2(i) - 1;
    wcorrah(ahrows) = wcorrah(ahrows) + realrightvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop creates a corrected reproductive value vector
  for (int i = 0; i < rvel; i++) {
    ahrows = h_stage_2(i) - 1;
    
    if (wcorrah(ahrows) != 0) {
      vcorrah(ahrows) = vwprod(i) / wcorrah(ahrows) + vcorrah(ahrows);
    } else {
      // Some stages have expected corrected stable stage proportions of 0
      vcorrah(ahrows) = 0.0 + vcorrah(ahrows);
    }
  }
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  // This next section creates the ahistorical sensitivity matrix
  for (int i = 0; i < ahstagelength; i++) {
    vwprodah(i) = wcorrah(i) * vcorrah(i);
  }
  double vwscalarah = sum(vwprodah);
  
  for (int i = 0; i < ahstagelength; i++) {
    for (int j = 0; j < ahstagelength; j++) {
      ahsens(i, j) = vcorrah(i) * wcorrah(j) / vwscalarah;
    }
  }
  
  List output = List::create(Named("h_smat") = smat, _["ah_smat"] = ahsens);
  
  return output;
}

//' Estimate Deterministic Elasticities of Any Population Matrix
//' 
//' \code{elas3matrix()} returns the elasticity of lambda with respect
//' to each element in a dense or sparse matrix, both provided in dense matrix
//' format. This is accomplished via the \code{eig_gen}() and \code{eigs_gen}()
//' functions in the C++ Armadillo library.
//' 
//' @name elas3matrix
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a matrix of deterministic elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3matrix)]]
arma::mat elas3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);

  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // Lower threshold than used in w and v
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;

  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // Lower threshold than used in w and v

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Deterministic Elasticities of a Historical LefkoMat Object
//' 
//' \code{elas3hlefko()} returns the elasticity of lambda with respect
//' to each historical stage-pair in the matrix, and the summed elasticities
//' for each life history stage. This is accomplished via the \code{eigs_gen}()
//' function in the C++ Armadillo library.
//' 
//' @name elas3hlefko
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic elasticity
//' matrices:
//' \item{h_emat}{Matrix of elasticities corresponding to the historical matrix.}
//' \item{ah_emat}{Matrix of elasticities corresponding to the ahistorical
//' matrix, but using summed historical elasticities as the basis of estimation.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3hlefko)]]
List elas3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  arma::uvec stage_id = as<arma::uvec>(ahstages["stage_id"]);
  arma::uvec h_stage_2 = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec h_stage_1 = as<arma::uvec>(hstages["stage_id_1"]);
  
  List eigenstuff = decomp3sp(Amat);
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;

  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // The next few lines set up the empty ahistorical matrix
  int ahstagelength = stage_id.n_elem;
  arma::mat ahelas(ahstagelength, ahstagelength);
  ahelas.zeros();
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
      ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) =
        ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) + emat(i, j);
    }
  }
  
  List output = List::create(Named("h_emat") = emat, _["ah_emat"] = ahelas);
  
  return output;
}

//' Core Time-based Population Matrix Projection Function
//' 
//' Function \code{proj3()} runs the matrix projections used in other functions
//' in package \code{lefko3}.
//' 
//' @name proj3
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{$A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3)]]
arma::mat proj3(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool standardize, bool growthonly, bool integeronly) {
  int sparse_switch {0};
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  arma::mat popproj(nostages, (theclairvoyant + 1)); // Population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1)); // Population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1)); // Population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1));
  popproj.zeros();
  wpopproj.zeros();
  vpopproj.zeros();
  Rvecmat.zeros();
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  arma::mat finaloutput;
  
  // Here we will check if the matrix is large and sparse
  int test_elems = as<arma::mat>(core_list(0)).n_elem;
  arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
  int all_nonzeros = nonzero_elems.n_elem;
  double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
  if (sparse_check <= 0.5 && start_vec.n_elem > 100) {
    sparse_switch = 1;
  } else sparse_switch = 0;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  if (sparse_switch == 0) {
    // Dense matrix projection
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
      
      theseventhson = theprophecy * theseventhson;
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      popproj.col(i+1) = theseventhson;
      Rvecmat(i+1) = sum(theseventhson);
      
      if (standardize) {
        theseventhson = theseventhson / sum(theseventhson);
      }
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        thesecondprophecy = as<arma::mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * thesecondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  } else {
    // Sparse matrix projection
    arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    
    int matlist_length = core_list.size();
    arma::mat first_mat = core_list(0);
    arma::sp_mat new_sparse = arma::sp_mat(first_mat);
    Rcpp::List sparse_list = List::create(_["1"] = new_sparse);
    if(matlist_length > 1) {
      for (int i = 1; i < matlist_length; i++) {
        first_mat = as<arma::mat>(core_list(i));
        new_sparse = arma::sp_mat(first_mat);
        sparse_list.push_back(new_sparse);
      }
    }
    arma::sp_mat sparse_prophecy;
    arma::sp_mat sparse_secondprophecy;
    
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      sparse_prophecy = as<arma::sp_mat>(sparse_list[(mat_order(i))]);
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
      Rvecmat(i+1) = sum(popproj.col(i+1));
      
      if (standardize) {
        sparse_seventhson = sparse_seventhson / sum(popproj.col(i+1));
      }
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        sparse_secondprophecy = as<arma::sp_mat>(sparse_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Slimmed-down Time-based Population Sparse Matrix Projection Function
//' 
//' Function \code{proj3sp()} runs the matrix projections used in some other
//' functions in package \code{lefko3}, but only when the input is sparse. This
//' is a slimmed down version of function \code{proj3()}
//' 
//' @name proj3sp
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to
//' the \code{$A} list within a \code{lefkoMat} object. Matrices must be in
//' \code{arma::sp_mat} format.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3sp)]]
arma::mat proj3sp(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool standardize, bool growthonly, bool integeronly) {
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  arma::mat popproj(nostages, (theclairvoyant + 1)); // Population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1)); // Population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1)); // Population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1));
  popproj.zeros();
  wpopproj.zeros();
  vpopproj.zeros();
  Rvecmat.zeros();
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
  arma::mat finaloutput;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  // Sparse matrix projection
  arma::sp_mat sparse_prophecy;
  arma::sp_mat sparse_secondprophecy;
    
  for (int i = 0; i < theclairvoyant; i++) {
    if (i % 50 == 0) Rcpp::checkUserInterrupt();
    
    sparse_prophecy = as<arma::sp_mat>(core_list[(mat_order(i))]);
    sparse_seventhson = sparse_prophecy * sparse_seventhson;
    if (integeronly) {
      sparse_seventhson = floor(sparse_seventhson);
    }
    popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
    Rvecmat(i+1) = sum(popproj.col(i+1));
    
    if (standardize) {
      sparse_seventhson = sparse_seventhson / sum(popproj.col(i+1));
    }
    
    if (!growthonly) {
      wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
      sparse_secondprophecy = as<arma::sp_mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
      theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
      
      double seventhgrandsum = sum(theseventhgrandson);
      arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
      theseventhgrandson = theseventhgrandson / seventhgrandsum;
      vpopproj.col(theclairvoyant - (i+1)) = midwife;
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Core Time-based Density-Dependent Population Matrix Projection Function
//' 
//' Function \code{proj3dens()} runs density-dependent matrix projections.
//' 
//' @name proj3dens
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0
//' and fecundity to be non-negative, and \code{2} forces all column rows to
//' total no more than 1.0.
//' @param dens_input The original \code{lefkoDens} data frame supplied through
//' the \code{\link{density_input}()} function.
//' @param dens_index A list giving the indices of elements in object
//' \code{dens_input}.
//' @param allow_warnings A logical value indicating whether the function should
//' send warnings if estimated values fall outside of the realm of possibility.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @section Notes:
//' There is no option to standardize population vectors here, because density
//' dependence requires the full population size to be tracked.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3dens)]]
arma::mat proj3dens(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool growthonly, bool integeronly, int substoch, Rcpp::DataFrame dens_input,
  Rcpp::List dens_index, bool allow_warnings = false) {
  int sparse_switch {0};
  int time_delay {1};
  double pop_size {0};
  bool warn_trigger_neg = false;
  bool warn_trigger_1 = false;
  
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  // Density dependence arguments
  arma::uvec dyn_index321 = as<arma::uvec>(dens_index["index321"]);
  arma::uvec dyn_index_col = as<arma::uvec>(dens_index[1]);
  arma::uvec dyn_style = as<arma::uvec>(dens_input["style"]);
  arma::vec dyn_alpha = as<arma::vec>(dens_input["alpha"]);
  arma::vec dyn_beta = as<arma::vec>(dens_input["beta"]);
  arma::uvec dyn_delay = as<arma::uvec>(dens_input["time_delay"]);
  arma::uvec dyn_type = as<arma::uvec>(dens_input["type"]);
  int n_dyn_elems = dyn_index321.n_elem;
  
  // Matrices and vectors for projection results
  arma::mat popproj(nostages, (theclairvoyant + 1)); // Population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1)); // Population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1)); // Population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1));
  popproj.zeros();
  wpopproj.zeros();
  vpopproj.zeros();
  Rvecmat.zeros();
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  arma::mat finaloutput;
  
  // Here we will check if the matrix is large and sparse
  int test_elems = as<arma::mat>(core_list(0)).n_elem;
  arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
  int all_nonzeros = nonzero_elems.n_elem;
  double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
  if (sparse_check <= 0.5 && theseventhson.n_elem > 100) {
    sparse_switch = 1;
  } else sparse_switch = 0;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  double changing_element {0.0};
  double changing_colsum {0.0};
  
  if (sparse_switch == 0) {
    // Dense matrix projection
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        if (i >= time_delay) {
          pop_size = sum(popproj.col(i - time_delay));
          
          if (dyn_style(j) == 1) { // Ricker
            changing_element = theprophecy(dyn_index321(j)) * 
              dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
            
          } else if (dyn_style(j) == 2) { // Beverton-Holt
            changing_element = theprophecy(dyn_index321(j)) * 
              dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
            
          } else if (dyn_style(j) == 3) { // Usher function
            changing_element = theprophecy(dyn_index321(j)) * 
              (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
            
          } else if (dyn_style(j) == 4) { // Logistic function
            double used_popsize = pop_size;
            if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
              used_popsize = dyn_alpha(j);
            }
            changing_element = theprophecy(dyn_index321(j)) * 
              (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
          }
          
          if (substoch == 1) {
            if (changing_element > 1.0 && dyn_type(j) == 1) {
              changing_element = 1.0;
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch == 2 && dyn_type(j) == 1) {
            double barnyard_antics {0.0};
            arma::vec given_col = theprophecy.col(dyn_index_col(j));
            arma::uvec gc_negs = find(given_col < 0.0);
            if (gc_negs.n_elem > 0) {
              barnyard_antics = sum(given_col(gc_negs));
            }
            changing_colsum = sum(given_col) - theprophecy(dyn_index321(j)) - barnyard_antics;
            
            if (changing_element > (1.0 - changing_colsum)) {
              changing_element = (1.0 - changing_colsum);
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch > 0 && dyn_type(j) == 2) {
            if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          }
          theprophecy(dyn_index321(j)) = changing_element;
          
          if (allow_warnings) {
            if (dyn_type(j) == 1 && theprophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              Rf_warningcall(R_NilValue, "Some probabilities with value > 1.0 produced during density adjustment.");
            } else if (theprophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
              warn_trigger_neg = true;
              Rf_warningcall(R_NilValue, "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
      }
      theseventhson = theprophecy * theseventhson; // thechosenone
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      popproj.col(i+1) = theseventhson;
      Rvecmat(i+1) = sum(theseventhson);
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        thesecondprophecy = as<arma::mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * thesecondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  } else {
    // Sparse matrix projection
    arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    int matlist_length = core_list.size();
    Rcpp::List sparse_list(matlist_length);
    
    arma::mat first_mat = core_list(0);
    arma::sp_mat new_sparse = arma::sp_mat(first_mat);
    sparse_list(0) = new_sparse;
    
    if(matlist_length > 1) {
      for (int i = 1; i < matlist_length; i++) {
        first_mat = as<arma::mat>(core_list(i));
        new_sparse = arma::sp_mat(first_mat);
        sparse_list(i) = new_sparse;
      }
    }
    
    arma::sp_mat sparse_prophecy;
    arma::sp_mat sparse_secondprophecy;
    
    for (int i = 0; i < theclairvoyant; i++) {
      if (i % 50 == 0) Rcpp::checkUserInterrupt();
      
      sparse_prophecy = as<arma::sp_mat>(sparse_list[(mat_order(i))]);
      
      for (int j = 0; j < n_dyn_elems; j++) { // Density dependence
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        if (i >= time_delay) {
          pop_size = sum(popproj.col(i - time_delay));
          
          if (dyn_style(j) == 1) { // Ricker
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
            
          } else if (dyn_style(j) == 2) { // Beverton-Holt
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
            
          } else if (dyn_style(j) == 3) { // Usher function
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
            
          } else if (dyn_style(j) == 4) { // Logistic function
            double used_popsize = pop_size;
            if (dyn_beta(j) > 0.0 && pop_size > dyn_alpha(j)) {
              used_popsize = dyn_alpha(j);
            }
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              (1 - used_popsize / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
          }
          
          if (substoch == 1 && dyn_type(j) == 1) {
            if (changing_element > 1.0) {
              changing_element = 1.0;
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch == 2 && dyn_type(j) == 1) {
            double barnyard_antics {0.0};
            arma::vec given_col = arma::vec(sparse_prophecy.col(dyn_index_col(j)));
            arma::uvec gc_negs = find(given_col < 0.0);
            if (gc_negs.n_elem > 0) {
              barnyard_antics = sum(given_col(gc_negs));
            }
            changing_colsum = sum(given_col) - sparse_prophecy(dyn_index321(j)) - barnyard_antics;
            
            if (changing_element > (1.0 - changing_colsum)) {
              changing_element = (1.0 - changing_colsum);
            } else if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          } else if (substoch > 0 && dyn_type(j) == 2) {
            if (changing_element < 0.0) {
              changing_element = 0.0;
            }
          }
          sparse_prophecy(dyn_index321(j)) = changing_element;
          
          if (allow_warnings) {
            if (dyn_type(j) == 1 && sparse_prophecy(dyn_index321(j)) > 1.0 && !warn_trigger_1) {
              warn_trigger_1 = true;
              Rf_warningcall(R_NilValue, "Some probabilities with value > 1.0 produced during density adjustment.");
            } else if (sparse_prophecy(dyn_index321(j)) < 0.0 && !warn_trigger_neg) {
              warn_trigger_neg = true;
              Rf_warningcall(R_NilValue, "Some matrix elements with value < 0.0 produced during density adjustment.");
            }
          }
        }
      }
      
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
      Rvecmat(i+1) = sum(popproj.col(i+1));
      
      if (!growthonly) {
        wpopproj.col(i+1) = popproj.col(i+1) / Rvecmat(i+1);
        sparse_secondprophecy = as<arma::sp_mat>(sparse_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Conduct Population Projection Simulations
//' 
//' Function \code{projection3()} runs projection simulations. It projects the
//' population and patches forward in time by a user-defined number of
//' occasions. A given set of matrices is utilized and not recreated, although
//' elements may be altered if density dependence is set. Projections may be
//' deterministic or stochastic, and may be density dependent in either case. If
//' deterministic, then projections will be cyclical if matrices exist covering
//' multiple occasions for each population or patch. If stochastic, then annual
//' matrices will be shuffled within patches and populations. Also produces
//' replicates if set.
//' 
//' @name projection3
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param nreps The number of replicate projections.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' 10,000.
//' @param historical An optional logical value only used if object \code{mpm}
//' is a list of matrices, rather than a \code{lefkoMat} object. Defaults to
//' \code{FALSE} for the former case, and overridden by information supplied in
//' the \code{lefkoMat} object for the latter case.
//' @param stochastic A logical value denoting whether to conduct a stochastic
//' projection or a deterministic / cyclical projection.
//' @param standardize A logical value denoting whether to re-standardize the
//' population size to 1.0 at each occasion. Defaults to \code{FALSE}.
//' @param growthonly A logical value indicating whether to produce only the
//' projected population size at each occasion, or a vector showing the stage
//' distribution followed by the reproductive value vector followed by the full
//' population size at each occasion. Defaults to \code{TRUE}.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion to the nearest
//' integer. Defaults to \code{FALSE}.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0,
//' and forces fecundity to be non-negative; and \code{2} forces all column rows
//' in the survival-transition matrices to total no more than 1.0, in addition
//' to the actions outlined for option \code{1}.
//' @param sub_warnings A logical value indicating whether to warn the user if
//' density dependence yields matrix values outside of the realm of possibility.
//' Generally, this means that survival-transition elements altered to values
//' outside of the interval [0, 1], and negative fecundity values, will both
//' yield warnings. Defaults to \code{TRUE}.
//' @param year Either a single integer value corresponding to the year to
//' project, or a vector of \code{times} elements with the year to use at each
//' time step. If a vector shorter than \code{times} is supplied, then this
//' vector will be cycled. If not provided, then all annual matrices will be
//' cycled within patches or populations.
//' @param start_vec An optional numeric vector denoting the starting stage
//' distribution for the projection. Defaults to a single individual of each
//' stage.
//' @param start_frame An optional data frame characterizing stages, age-stages,
//' or stage-pairs that should be set to non-zero values in the starting vector,
//' and what those values should be. Can only be used with \code{lefkoMat}
//' objects.
//' @param tweights An optional numeric vector denoting the probabilistic
//' weightings of annual matrices. Defaults to equal weighting among occasions.
//' @param density An optional data frame describing the matrix elements that
//' will be subject to density dependence, and the exact kind of density
//' dependence that they will be subject to. The data frame used should be an
//' object of class \code{lefkoDens}, which is the output from function
//' \code{\link{density_input}()}.
//' 
//' @return A list of class \code{lefkoProj}, which always includes the first
//' three elements of the following, and also includes the remaining elements
//' below when a \code{lefkoMat} object is used as input:
//' \item{projection}{A list of lists of matrices showing the total number of
//' individuals per stage per occasion. The first list corresponds to each
//' pop-patch followed by each population. The inner list corresponds to
//' replicates within each pop-patch or population.}
//' \item{stage_dist}{A list of lists of the actual stage distribution in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{rep_value}{A list of lists of the actual reproductive value in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{pop_size}{A list of data frames showing the total population size in
//' each occasion per replicate (row within data frame) per pop-patch or
//' population (list element).}
//' \item{labels}{A data frame showing the order of populations and patches in
//' item \code{projection}.}
//' \item{ahstages}{The original stageframe used in the study.}
//' \item{hstages}{A data frame showing the order of historical stage pairs.}
//' \item{agestages}{A data frame showing the order of age-stage pairs.}
//' \item{control}{A short vector indicating the number of replicates and the
//' number of occasions projected per replicate.}
//' \item{density}{The data frame input under the density option. Only provided
//' if input by the user.}
//' 
//' @section Notes:
//' Projections are run both at the patch level and at the population level.
//' Population level estimates will be noted at the end of the
//' data frame with 0 entries for patch designation.
//' 
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//' 
//' Starting vectors can be input in one of two ways: 1) as \code{start_vec}
//' input, which is a vector of numbers of the numbers of individuals in each
//' stage, stage pair, or age-stage, with the length of the vector necessarily
//' as long as there are rows in the matrices of the MPM; or 2) as
//' \code{start_frame} input, which is a data frame showing only those stages,
//' stage pairs, or age-stages that should begin with more than 0 individuals,
//' and the numbers of individuals that those stages should start with (this
//' object is created using the \code{\link{start_input}()} function). If both
//' are provided, then \code{start_frame} takes precedence and \code{start_vec}
//' is ignored. If neither is provided, then \code{projection3()} automatically
//' assumes that each stage, stage pair, or age-stage begins with a single
//' individual. Importantly, if a \code{lefkoMat} object is not used, and a list
//' of matrices is provided instead, then \code{start_frame} cannot be utilized
//' and a full \code{start_vec} must be provided to conduct a simulation with
//' starting numbers of individuals other than 1 per stage.
//' 
//' The resulting data frames in element \code{projection} are separated by
//' pop-patch according to the order provided in element \code{labels}, but the
//' matrices for each element of \code{projection} have the result of each
//' replicate stacked in order on top of one another without any break or
//' indication. Results for each replicate must be separated using the
//' information provided in elements \code{control} and the 3 stage
//' descriptor elements.
//' 
//' Density dependent projections are automatically set up if object
//' \code{density} is input. If this object is not included, then density
//' independent projections will be set up. Note that currently, density
//' dependent projections can only be performed with \code{lefkoMat} objects.
//' 
//' The stage distributions and reproductive values produced are not the
//' asymptotic values as would be given by the standardized right and left
//' eigenvectors associated with the dominant eigenvalue of a matrix, but are
//' vectors describing these values at the specific points in time projected.
//' See equations 14.86 and 14.88 and section 14.4 on Sensitivity and Elasticity
//' Analysis under Environmental Stochasticity in Caswell (2001, Matrix
//' Population Models, Sinauer Associates) for more details.
//' 
//' Consistently positive population growth can quickly lead to population size
//' numbers larger than can be handled computationally. In that circumstance, a
//' continuously rising population size will suddenly become \code{NaN} for the
//' remainder of the projection.
//' 
//' Users wishing to run a projection of a single patch in a \code{lefkoMat}
//' object with multiple patches should subset the MPM first to contain only
//' the patch needed. This can be accomplished with the
//' \code{\link{subset_lM}()} function.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{density_input}()}
//' @seealso \code{\link{f_projection3}()}
//' 
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
//'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
//'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
//' 
//' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
//' 
//' # Cypripedium example
//' data(cypdata)
//'  
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4, 
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04", 
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE, 
//'   NRasRep = TRUE)
//' 
//' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
//'     "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
//'     "SL", "SL", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
//'     "SL", "SL", "SL", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
//'     NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
//'     NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"), 
//'   supplement = cypsupp3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- projection3(cypmatrix3r, nreps = 5, stochastic = TRUE)
//' 
//' @export projection3
// [[Rcpp::export]]
Rcpp::List projection3(List mpm, int nreps = 1, int times = 10000,
  bool historical = false, bool stochastic = false, bool standardize = false,
  bool growthonly = true, bool integeronly = false, int substoch = 0,
  bool sub_warnings = true, Nullable<IntegerVector> year = R_NilValue,
  Nullable<NumericVector> start_vec = R_NilValue, Nullable<DataFrame> start_frame = R_NilValue,
  Nullable<NumericVector> tweights = R_NilValue, Nullable<DataFrame> density = R_NilValue) {
  
  Rcpp::List dens_index;
  Rcpp::DataFrame dens_input;
  
  int theclairvoyant = times;
  int dens_switch {0};
  int used_matsize {0};
  int total_projrows {0};
  bool year_override = false;
  
  if (theclairvoyant < 1) {
    throw Rcpp::exception("Option times must be a positive integer.", false);
  }
  if (nreps < 1) {
    throw Rcpp::exception("Option nreps must be a positive integer.", false);
  }
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Option substoch must be set to 0, 1, or 2.", false);
  }
  
  arma::uvec theprophecy(theclairvoyant);
  theprophecy.zeros();
  
  arma::vec startvec;
  arma::mat projection;
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = as<List>(mpm["A"]);
    List umats = as<List>(mpm["U"]);
    List fmats = as<List>(mpm["F"]);
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    historical = false;
    bool agebystage = false;
    
    if (hstages.length() > 1) {
      historical = true;
    }
    if (agestages.length() > 1) {
      agebystage = true;
    }
    
    if (density.isNotNull()) { 
      Rcpp::DataFrame dens_thru(density);
      dens_input = dens_thru;
      dens_switch = 1;
      
      Rcpp::StringVector di_stage3 = as<StringVector>(dens_input["stage3"]);
      Rcpp::StringVector di_stage2 = as<StringVector>(dens_input["stage2"]);
      Rcpp::StringVector di_stage1 = as<StringVector>(dens_input["stage1"]);
      int di_size = di_stage3.length();
      
      if (historical) {
        StringVector stage3 = as<StringVector>(hstages["stage_2"]);
        StringVector stage2r = as<StringVector>(hstages["stage_1"]);
        StringVector stage2c = as<StringVector>(hstages["stage_2"]);
        StringVector stage1 = as<StringVector>(hstages["stage_1"]);
        int hst_size = stage3.length();
        
        arma::uvec hst_3(hst_size);
        arma::uvec hst_2r(hst_size);
        arma::uvec hst_2c(hst_size);
        arma::uvec hst_1(hst_size);
        hst_3.zeros();
        hst_2r.zeros();
        hst_2c.zeros();
        hst_1.zeros();
        
        arma::uvec di_stage32_id(di_size);
        arma::uvec di_stage21_id(di_size);
        arma::uvec di_index(di_size);
        di_stage32_id.zeros();
        di_stage21_id.zeros();
        di_index.zeros();
        
        for (int i = 0; i < di_size; i++) { // This loop runs through each density_input line
          for (int j = 0; j < hst_size; j++) {
            if (di_stage3(i) == stage3(j)) {
              hst_3(j) = 1;
            } else {
              hst_3(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage2(i) == stage2r(j)) {
              hst_2r(j) = 1;
            } else {
              hst_2r(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage2(i) == stage2c(j)) {
              hst_2c(j) = 1;
            } else {
              hst_2c(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage1(i) == stage1(j)) {
              hst_1(j) = 1;
            } else {
              hst_1(j) = 0;
            }
          }
          
          arma::uvec find_hst3 = find(hst_3);
          arma::uvec find_hst2r = find(hst_2r);
          arma::uvec find_hst2c = find(hst_2c);
          arma::uvec find_hst1 = find(hst_1);
          
          arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
          arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
          
          di_stage32_id(i) = pop_32(0);
          di_stage21_id(i) = pop_21(0);
          di_index(i) = pop_32(0) + (pop_21(0) * hst_size);
          
          hst_3.zeros();
          hst_2r.zeros();
          hst_2c.zeros();
          hst_1.zeros();
        }
        
        dens_index = Rcpp::List::create(_["index32"] = di_stage32_id,
          _["index21"] = di_stage21_id, _["index321"] = di_index);
        
      } else {
        StringVector stage3 = as<StringVector>(stageframe["stage"]);
        StringVector stage2 = as<StringVector>(stageframe["stage"]);
        int ahst_size = stage3.length();
        
        arma::uvec ahst_3(ahst_size);
        arma::uvec ahst_2(ahst_size);
        ahst_3.zeros();
        ahst_2.zeros();

        arma::uvec di_stage32_id(di_size);
        arma::uvec di_stage21_id(di_size);
        arma::uvec di_index(di_size);
        di_stage32_id.zeros();
        di_stage21_id.zeros();
        di_index.zeros();
        
        for (int i = 0; i < di_size; i++) { // This loop runs through each density_input line
          for (int j = 0; j < ahst_size; j++) {
            if (di_stage3(i) == stage3(j)) {
              ahst_3(j) = 1;
            } else {
              ahst_3(j) = 0;
            }
          }
          
          for (int j = 0; j < ahst_size; j++) {
            if (di_stage2(i) == stage2(j)) {
              ahst_2(j) = 1;
            } else {
              ahst_2(j) = 0;
            }
          }
          
          arma::uvec find_ahst3 = find(ahst_3);
          arma::uvec find_ahst2 = find(ahst_2);
          di_stage32_id(i) = find_ahst3(0);
          di_stage21_id(i) = find_ahst2(0);
          di_index(i) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
          
          ahst_3.zeros();
          ahst_2.zeros();
        }
        
        dens_index = Rcpp::List::create(_["index3"] = di_stage32_id,
          _["index2"] = di_stage21_id, _["index321"] = di_index);
      }
    }
    
    IntegerVector yearorder;
    StringVector patchorder;
    if (labels.length() < 3) {
      StringVector label_elements = labels.attr("names");
      std::string patch_named = "patch";
      
      for (int i = 0; i < label_elements.length(); i++) {
        if (stringcompare_hard(as<std::string>(label_elements(i)), "patch")) {
          Rf_warningcall(R_NilValue, "This function takes annual matrices as input. This lefkoMat object appears to be a set of mean matrices, and may lack annual matrices. Will project only the mean.");
        }
      }
      
      StringVector patch_projected = as<StringVector>(labels["patch"]);
      IntegerVector years_projected(patch_projected.length());
      for (int i = 0; i < patch_projected.length(); i++) {
        years_projected(i) = 1;
      }
      
      patchorder = patch_projected;
      yearorder = years_projected;
    } else {
      patchorder = as<StringVector>(labels["patch"]);
      yearorder = as<IntegerVector>(labels["year2"]);
    }
    StringVector poporder = as<StringVector>(labels["pop"]);
    arma::uvec armayearorder = as<arma::uvec>(yearorder);
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    IntegerVector uniqueyears = sort_unique(yearorder);
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    IntegerVector years_forward;
    if (year.isNotNull()) {
      if (stochastic) throw Rcpp::exception("Options year cannot be used when stochastic = TRUE.", false);
      
      IntegerVector years_ = as<IntegerVector>(year);
      
      int member_sum {0};
      for (int i = 0; i < years_.length(); i++) {
        for (int j = 0; j < yl; j++) {
          if (years_[i] == uniqueyears[j]) member_sum++;
        }
        if (member_sum == 0) {
          throw Rcpp::exception("Option year includes time indices that do not exist in the input lefkoMat object.", false);
        }
        member_sum = 0;
      }
      
      IntegerVector years_pre (times);
      
      int rampant_exigence {0};
      for (int i = 0; i < times; i++) {
        years_pre(i) = years_(rampant_exigence);
        rampant_exigence++;
        
        if (rampant_exigence >= years_.length()) {
          rampant_exigence = 0;
        }
      }
      years_forward = years_pre; // This is the programmed order of matrices for all times, if years is input
      year_override = true; // This variable decides whether to use years or the defaults matrix vectors
    }
    
    arma::vec twinput;
    if (tweights.isNotNull()) {
      twinput = as<arma::vec>(tweights);
      if (twinput.n_elem != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      if (!stochastic) throw Rcpp::exception("Option tweights can only be used when stochastic = TRUE.", false);
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize);
    arma::uvec yearsinpatch(loysize);
    patchesinpop.zeros();
    yearsinpatch.zeros();
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = summervacation.n_elem;
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = motorhead.n_elem;
      
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    List mean_lefkomat;
    
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the
    // simulation, and estimate all descriptive metrics
    List meanamats = as<List>(mean_lefkomat["A"]);
    List mmlabels = as<List>(mean_lefkomat["labels"]);
    StringVector mmpops = as<StringVector>(mmlabels["pop"]);
    StringVector mmpatches = as<StringVector>(mmlabels["patch"]);
    
    arma::mat thechosenone = as<arma::mat>(meanamats[0]);
    int meanmatsize = thechosenone.n_elem;
    int meanmatrows = thechosenone.n_rows;
    arma::vec startvec;
    int trials = meanamats.length();
    used_matsize = meanmatrows;
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = allppcs.n_elem;
    List plist_hold(allppcsnem);
    int pop_est {1};
    if (allppcsnem > 1) {
      pop_est = trials - allppcsnem;    
    }
    List projection_list(trials);
    
    if(start_frame.isNotNull()) {
      Rcpp::DataFrame start_thru(start_frame);
      startvec.set_size(meanmatrows);
      startvec.zeros();
      arma::uvec start_elems = as<arma::uvec>(start_thru["row_num"]);
      start_elems = start_elems - 1;
      arma::vec start_values = as<arma::vec>(start_thru["value"]);
      
      if (start_elems.max() > (meanmatrows - 1)) {
        throw Rcpp::exception("Start vector input frame includes element indices too high for this MPM.",
          false);
      }
      for (int i = 0; i < start_elems.n_elem; i++) {
        startvec(start_elems(i)) = start_values(i);
      }
      
    } else if (start_vec.isNotNull()) {
      startvec = as<arma::vec>(start_vec);
      if (startvec.n_elem != meanmatrows) {
        throw Rcpp::exception("Start vector must be the same length as the number of rows in each matrix.",
          false);
      }
      
    } else {
      startvec.set_size(meanmatrows);
      startvec.ones();
    }
    
    twinput = twinput / sum(twinput);
    
    for (int i= 0; i < allppcsnem; i++) {
      thechosenone = as<arma::mat>(meanamats[i]);
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      int chosen_yl = thenumbersofthebeast.n_elem;
      
      arma::uvec pre_prophecy (theclairvoyant, fill::zeros);
      if (year_override) {
        for (int j = 0; j < theclairvoyant; j++) {
          arma::uvec tnb_year_indices = find(armayearorder == years_forward(j));
          arma::uvec year_patch_intersect = intersect(thenumbersofthebeast, tnb_year_indices);
          
          pre_prophecy(j) = year_patch_intersect(0);
        }
      }
      // This loop takes care of multiple replicates, creating the final data frame
      // of results for each pop-patch
      for (int rep = 0; rep < nreps; rep++) {
        if (stochastic) {
          theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
          
        } else if (year_override) {
          theprophecy = pre_prophecy;
          
        } else {
          theprophecy.set_size(theclairvoyant);
          theprophecy.zeros();
          
          for (int j = 0; j < theclairvoyant; j++) {
            theprophecy(j) = thenumbersofthebeast(j % chosen_yl);
          }
        }
        
        if (dens_switch) {
          if (rep == 0) {
            projection = proj3dens(startvec, amats, theprophecy, growthonly,
              integeronly, substoch, dens_input, dens_index, sub_warnings);
          } else {
            arma::mat nextproj = proj3dens(startvec, amats, theprophecy,
              growthonly, integeronly, substoch, dens_input, dens_index,
              sub_warnings);
            projection = arma::join_cols(projection, nextproj);
          }
        } else {
          if (rep == 0) {
            projection = proj3(startvec, amats, theprophecy, standardize, growthonly,
              integeronly);
          } else {
            arma::mat nextproj = proj3(startvec, amats, theprophecy, standardize,
              growthonly, integeronly);
            projection = arma::join_cols(projection, nextproj);
          }
        }
      }
      
      projection_list(i) = projection;
    }
    
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize);
    arma::uvec yearmatch(loysize);
    popmatch.zeros();
    yearmatch.zeros();
    List meanmatyearlist(uniqueyears.length());
    
    if (allppcsnem > 1) { // Checks for pop-mean matrices separate from the the patch means
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // This loop goes through each population
        thechosenone = as<arma::mat>(meanamats[allppcsnem + i]);
        
        for (int j = 0; j < loysize; j++) { // Checks which A matrices match current population
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        arma::uvec neededmatspop = find(popmatch);
        
        for (int j = 0; j < yl; j++) { // Checks each year and develops matrix mean across patches
          for (int k = 0; k < loysize; k++) { // Develops vector to find all matrices for current year
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          
          // Catches matrix indices matching current year and pop
          int crankybankynem = crankybanky.n_elem;
          arma::mat crossmat(meanmatsize, crankybankynem);
          crossmat.zeros();
          for (int j = 0; j < crankybankynem; j++) {
            crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
          }
          
          arma::vec happymedium(meanmatsize);
          happymedium.zeros();
          for (int j = 0; j < meanmatsize; j++) {
            for (int k = 0; k < crankybankynem; k++) {
              happymedium(j) = happymedium(j) + crossmat(j, k) / (crankybankynem);
            }
          }
          arma::mat finalyearmat = happymedium;
          finalyearmat.reshape(meanmatrows, meanmatrows);
          meanmatyearlist(j) = finalyearmat;
        }
        
        int numyearsused = meanmatyearlist.length();
        arma::uvec choicevec = linspace<arma::uvec>(0, (numyearsused - 1), numyearsused);
        int chosen_yl = choicevec.n_elem;
      
        // This loop takes care of multiple replicates, creating the final data frame
        // of results for the pop mean(s)
        for (int rep = 0; rep < nreps; rep++) {
          if (stochastic) {
            theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant, true, twinput);
          } else {
            theprophecy.zeros();
            for (int j = 0; j < theclairvoyant; j++) {
              theprophecy(j) = choicevec(j % chosen_yl);
            }
          }
          
          if (dens_switch) {
            if (rep == 0) {
              projection = proj3dens(startvec, meanmatyearlist, theprophecy,
                growthonly, integeronly, substoch, dens_input, dens_index,
                sub_warnings);
            } else {
              arma::mat nextproj = proj3dens(startvec, meanmatyearlist, theprophecy,
                growthonly, integeronly, substoch, dens_input, dens_index,
                sub_warnings);
              projection = arma::join_cols(projection, nextproj);
            }
          } else {
            if (rep == 0) {
              projection = proj3(startvec, meanmatyearlist, theprophecy,
                standardize, growthonly, integeronly);
            } else {
              arma::mat nextproj = proj3(startvec, meanmatyearlist, theprophecy,
                standardize, growthonly, integeronly);
              projection = arma::join_cols(projection, nextproj);
            }
          }
        }
        projection_list(allppcsnem + i) = projection;
      }
    }
    
    // The final output will have a projection list with # elements = nreps, nested
    // within a list with # elements = # pop-patches
    List projection_set(nreps);
    List ss_set(nreps);
    List rv_set(nreps);
    arma::mat total_sizes_set(nreps, (times+1), fill::zeros);
    
    int length_ppy = projection_list.length();
    List final_projection(length_ppy);
    List final_ss(length_ppy);
    List final_rv(length_ppy);
    List final_ns(length_ppy);
    
    List output(9);
    
    if (!growthonly) {
      arma::mat list_proj(total_projrows, (times+1), fill::zeros);
      arma::mat extracted_proj(used_matsize, used_matsize, fill::zeros);
      int diversion = used_matsize * 3 + 1;
      
      for (int j = 0; j < length_ppy; j++) {
        list_proj = as<arma::mat>(projection_list[j]);
        
        for (int i = 0; i < nreps; i++) {
          extracted_proj = list_proj.rows((diversion * i), (diversion * i + used_matsize - 1));
          projection_set(i) = extracted_proj;
          
          extracted_proj = list_proj.rows((diversion * i + used_matsize),
            (diversion * i + (2 * used_matsize) - 1));
          ss_set(i) = extracted_proj;
          
          extracted_proj = list_proj.rows((diversion * i + (2 * used_matsize)),
            (diversion * i + (3 * used_matsize) - 1));
          rv_set(i) = extracted_proj;
          
          total_sizes_set.row(i) = list_proj.row(diversion * (i+1) - 1);
        }
        final_projection(j) = clone(projection_set);
        final_ss(j) = clone(ss_set);
        final_rv(j) = clone(rv_set);
        final_ns(j) = total_sizes_set;
      }
      
    } else {
      arma::mat list_proj(total_projrows, (times+1), fill::zeros);
      arma::mat extracted_proj(used_matsize, used_matsize, fill::zeros);
      int diversion = used_matsize;
      
      for (int j = 0; j < length_ppy; j++) {
        list_proj = as<arma::mat>(projection_list[j]);
        
        for (int i = 0; i < nreps; i++) {
          extracted_proj = list_proj.rows((diversion * i), (diversion * i + used_matsize - 1));
          projection_set(i) = extracted_proj;
          
          total_sizes_set.row(i) = sum(extracted_proj, 0);
        }
        final_projection(j) = clone(projection_set);
        final_ss(j) = NULL;
        final_rv(j) = NULL;
        final_ns(j) = total_sizes_set;
      }
    }
    
    DataFrame newlabels = DataFrame::create(_["pop"] = mmpops,
      _["patch"] = mmpatches);
    Rcpp::IntegerVector control = {nreps, times};
    
    output(0) = final_projection;
    output(1) = final_ss;
    output(2) = final_rv;
    output(3) = final_ns;
    output(4) = newlabels;
    output(5) = stageframe;
    output(6) = hstages;
    output(7) = agestages;
    output(8) = control;
    
    if (dens_switch) {
      output.push_back(dens_input);
      
      CharacterVector namevec = {"projection", "stage_dist", "rep_value", "pop_size",
        "labels", "ahstages", "hstages", "agestages", "control", "density"};
      output.attr("names") = namevec;
    } else {
      CharacterVector namevec = {"projection", "stage_dist", "rep_value", "pop_size",
        "labels", "ahstages", "hstages", "agestages", "control"};
      output.attr("names") = namevec;
    }
    output.attr("class") = "lefkoProj";
    
    return output;
    
  } else { // When not a lefkoMat object...
    List projection_list (1);
    List amats = mpm;
    int yl = amats.length();
    arma::mat firstmat = as<arma::mat>(amats[0]);
    int matrows = firstmat.n_rows;
    int matcols = firstmat.n_cols;
    used_matsize = matrows;
    
    arma::uvec uniqueyears(yl);
    for (int i = 0; i < yl; i++) {
      uniqueyears(i) = i;
    }
    arma::vec twinput;
    
    if (matrows != matcols) {
      throw Rcpp::exception("Supplied matrices must be square. Please check matrix dimensions.",
        false);
    }
    
    if (tweights.isNotNull()) {
      twinput = as<arma::vec>(tweights);
      if (twinput.n_elem != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    if (start_vec.isNotNull()) {
      startvec = as<arma::vec>(start_vec);
      if (startvec.n_elem != matrows) {
        throw Rcpp::exception("Start vector must be the same length as the number of rows in each matrix.", 
          false);
      }
      
    } else {
      startvec.set_size(matrows);
      startvec.ones();
    }
    
    IntegerVector years_forward;
    if (year.isNotNull()) {
      if (stochastic) throw Rcpp::exception("Options year cannot be used when stochastic = TRUE.", false);
      
      IntegerVector years_ = as<IntegerVector>(year);
      years_ = years_ - 1;
      
      int member_sum {0};
      for (int i = 0; i < years_.length(); i++) {
        for (int j = 0; j < yl; j++) {
          if (years_[i] == uniqueyears[j]) member_sum++;
        }
        if (member_sum == 0) {
          throw Rcpp::exception("Option year includes time indices that do not exist in the input lefkoMat object.", false);
        }
        member_sum = 0;
      }
      
      IntegerVector years_pre (times);
      
      int rampant_exigence {0};
      for (int i = 0; i < times; i++) {
        years_pre(i) = years_(rampant_exigence);
        rampant_exigence++;
        
        if (rampant_exigence >= years_.length()) {
          rampant_exigence = 0;
        }
      }
      years_forward = years_pre; // This is the programmed order of matrices for all times, if years is input
      year_override = true; // This variable decides whether to use years or the defaults matrix vectors
    }
    
    // Now we create the mean matrix
    arma::mat thechosenone(matrows, matcols);
    thechosenone.zeros();
    
    for (int i = 0; i < yl; i++) {
      arma::mat columnified = as<arma::mat>(amats[i]);
      thechosenone = thechosenone + (columnified / yl);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the simulation, and
    // estimate all descriptive metrics
    twinput = twinput / sum(twinput);
    arma::uvec thenumbersofthebeast = uniqueyears;
    
    // Here we loop multiple replicates, creating a data frame of results
    for (int rep = 0; rep < nreps; rep++) {
      if (stochastic) {
        theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
        
      } else if (year_override) {
        theprophecy = as<arma::uvec>(years_forward);
        
      } else {
        theprophecy.zeros();
        for (int i = 0; i < theclairvoyant; i++) {
          theprophecy(i) = thenumbersofthebeast(i % yl);
        }
      }
      
      if (rep == 0) {
        projection = proj3(startvec, amats, theprophecy, standardize, growthonly, integeronly);
        
      } else {
        arma::mat nextproj = proj3(startvec, amats, theprophecy, standardize, growthonly, integeronly);
        projection = arma::join_cols(projection, nextproj);
      }
    }
    projection_list(0) = projection;
    DataFrame newlabels = DataFrame::create(_["pop"] = 1, _["patch"] = 1);
    Rcpp::IntegerVector control = {nreps, times};
    
    Rcpp::List output = List::create(_["projection"] = projection_list, _["labels"] = newlabels,
      _["control"] = control);
    output.attr("class") = "lefkoProj";
    
    return output;
  }
}

//' Estimate Stochastic Population Growth Rate
//' 
//' Function \code{slambda3()} estimates the stochastic population growth rate,
//' \eqn{a}, defined as the long-term arithmetic mean of the log population 
//' growth rate estimated per simulated occasion. This function can handle both
//' lefkoMat objects and lists of full A matrices as input. 
//' 
//' @name slambda3
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param times Number of occasions to iterate. Defaults to \code{10000}.
//' @param historical An optional logical value only used if object \code{mpm}
//' is a list of matrices, rather than a \code{lefkoMat} object. Defaults to
//' \code{FALSE} for the former case, and overridden by information supplied in
//' the \code{lefkoMat} object for the latter case.
//' @param dense_only A logical value indicating whether to force matrices to be
//' run in dense format. Defaults to \code{FALSE}, and should only be used if
//' errors occur when running under default conditions.
//' @param tweights Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' 
//' @return A data frame with the following variables:
//' 
//' \item{pop}{The identity of the population.}
//' \item{patch}{The identity of the patch.}
//' \item{a}{Estimate of stochastic growth rate, estimated as the arithmetic
//' mean of the log population growth rate across simulated occasions.}
//' \item{var}{The estimated variance of a.}
//' \item{sd}{The standard deviation of a.}
//' \item{se}{The standard error of a.}
//'
//' @section Notes:
//' The log stochastic population growth rate, \eqn{a}, is as given in equation
//' 2 of Tuljapurkar, Horvitz, and Pascarella 2003. This term is estimated via
//' projection of randomly sampled matrices, similarly to the procedure outlined
//' in Box 7.4 of Morris and Doak (2002).
//'  
//' Stochastic growth rate is estimated both at the patch level and at the
//' population level. Population level estimates will be noted at the end of the
//' data frame with 0 entries for patch designation.
//' 
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//'
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
//'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
//'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
//' 
//' slambda3(ehrlen3)
//' 
//' # Cypripedium example
//' data(cypdata)
//'  
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4, 
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04", 
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE, 
//'   NRasRep = TRUE)
//' 
//' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
//'     "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
//'     "SL", "SL", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
//'     "SL", "SL", "SL", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
//'     NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
//'     NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"), 
//'   supplement = cypsupp3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- slambda3(cypmatrix3r, dense_only = TRUE)
//' cypstoch
//' 
//' @export slambda3
// [[Rcpp::export]]
DataFrame slambda3(List mpm, int times = 10000, bool historical = false,
  bool dense_only = false, Nullable<NumericVector> tweights = R_NilValue) {
  
  int theclairvoyant {0};
  int sparse_switch {0};
  theclairvoyant = times;
  
  if (theclairvoyant < 1) {
    throw Rcpp::exception("Option must equal a positive integer.", false);
  }
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = as<List>(mpm["A"]);
    List umats = as<List>(mpm["U"]);
    List fmats = as<List>(mpm["F"]);
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    if (hstages.length() > 1) {
      historical = true;
    } else {
      historical = false;
    }
    
    // Here we will check if the matrix is large and sparse
    int test_elems = as<arma::mat>(amats(0)).n_elem;
    int Amatrows = as<arma::mat>(amats(0)).n_rows;
    arma::uvec nonzero_elems = find(as<arma::mat>(amats(0)));
    int all_nonzeros = nonzero_elems.n_elem;
    double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
    if (sparse_check <= 0.5 && Amatrows > 100) {
      sparse_switch = 1;
    } else sparse_switch = 0;
    if (dense_only) sparse_switch = 0;
    
    if (labels.length() < 3) {
      throw Rcpp::exception("Function 'slambda3' requires annual matrices. This lefkoMat object appears to be a set of mean matrices, and lacks annual matrices.", false);
    }
    
    StringVector poporder = as<StringVector>(labels["pop"]);
    StringVector patchorder = as<StringVector>(labels["patch"]);
    IntegerVector yearorder = as<IntegerVector>(labels["year2"]);
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    IntegerVector uniqueyears = sort_unique(yearorder);
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    arma::vec twinput;
    if (tweights.isNotNull()) {
      twinput = as<arma::vec>(tweights);
      if (twinput.n_elem != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize);
    arma::uvec yearsinpatch(loysize);
    patchesinpop.zeros();
    yearsinpatch.zeros();
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = summervacation.n_elem;
      
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = motorhead.n_elem;
      
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    List mean_lefkomat;
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1);
    }
    
    // Here we take matrices for each patch, run the simulation, and estimate metrics
    List meanamats = as<List>(mean_lefkomat["A"]);
    List mmlabels = as<List>(mean_lefkomat["labels"]);
    StringVector mmpops = as<StringVector>(mmlabels["pop"]);
    StringVector mmpatches = as<StringVector>(mmlabels["patch"]);
    
    int meanmatsize = as<arma::mat>(meanamats[0]).n_elem; // thechosenone
    int meanmatrows = as<arma::mat>(meanamats[0]).n_rows; // thechosenone
    arma::vec startvec;
    int trials = meanamats.length();
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = allppcs.n_elem;
    
    arma::mat slmat(theclairvoyant, trials);
    slmat.zeros();
    arma::vec sl_mean(trials);
    arma::vec sl_var(trials);
    arma::vec sl_sd(trials);
    arma::vec sl_se(trials);
    sl_mean.zeros();
    sl_var.zeros();
    sl_sd.zeros();
    sl_se.zeros();
    
    twinput = twinput / sum(twinput);
    
    for (int i= 0; i < allppcsnem; i++) {
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
      startvec = ss3matrix(as<arma::mat>(meanamats[i]), sparse_switch); // thechosenone
      arma::mat projection = proj3(startvec, amats, theprophecy, 1, 1, 0);
      
      for (int j = 0; j < theclairvoyant; j++) {
        double madness = sum(projection.col(j+1));
        slmat(j,i) = log(madness);
      }
      
      sl_mean(i) = mean(slmat.col(i));
      sl_var(i) = var(slmat.col(i));
      sl_sd(i) = stddev(slmat.col(i));
      sl_se(i) = sl_sd(i) / sqrt(static_cast<double>(theclairvoyant));
    }
    
    int pop_est {1};
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize);
    arma::uvec yearmatch(loysize);
    popmatch.zeros();
    yearmatch.zeros();
    List meanmatyearlist(uniqueyears.length());
    
    if (allppcsnem > 1) { // Checks for pop-mean matrices separate from patch means
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // This loop goes through each population
        startvec = ss3matrix(as<arma::mat>(meanamats[allppcsnem + i]), sparse_switch); // thechosenone
        
        for (int j = 0; j < loysize; j++) { // Checks which A matrices match current pop
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        arma::uvec neededmatspop = find(popmatch);
        
        for (int j = 0; j < yl; j++) { // Checks each year and develops matrix mean across patches
          for (int k = 0; k < loysize; k++) { // Develops vector to find all matrices for current year
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          
          // This vector catches the indices of matrices that match the current year and pop
          int crankybankynem = crankybanky.n_elem;
          arma::mat crossmat(meanmatsize, crankybankynem);
          crossmat.zeros();
          for (int j = 0; j < crankybankynem; j++) {
            crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
          }
          
          arma::vec happymedium(meanmatsize);
          happymedium.zeros();
          for (int j = 0; j < meanmatsize; j++) {
            for (int k = 0; k < crankybankynem; k++) {
              happymedium(j) = happymedium(j) + crossmat(j, k) / (crankybankynem);
            }
          }
          
          arma::mat finalyearmat = happymedium;
          finalyearmat.reshape(meanmatrows, meanmatrows);
          meanmatyearlist(j) = finalyearmat;
        }
        
        int numyearsused = meanmatyearlist.length();
        arma::uvec choicevec = linspace<arma::uvec>(0, (numyearsused - 1), numyearsused);
        arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant, true, twinput);
        
        arma::mat projection = proj3(startvec, meanmatyearlist, theprophecy, 1, 1, 0);
        
        for (int j = 0; j < theclairvoyant; j++) {
          double madness = sum(projection.col(j+1));
          slmat(j,(allppcsnem +i)) = log(madness);
        }
        
        sl_mean((allppcsnem +i)) = mean(slmat.col((allppcsnem +i)));
        sl_var((allppcsnem +i)) = var(slmat.col((allppcsnem +i)));
        sl_sd((allppcsnem +i)) = stddev(slmat.col((allppcsnem +i)));
        sl_se((allppcsnem +i)) = sl_sd((allppcsnem +i)) / sqrt(static_cast<double>(theclairvoyant));
      }
    }
    return DataFrame::create(_["pop"] = mmpops, _["patch"] = mmpatches,
      _["a"] = sl_mean, _["var"] = sl_var, _["sd"] = sl_sd, _["se"] = sl_se);
    
  } else {
    List amats = mpm;
    
    // Here we will check if the matrix is large and sparse
    int test_elems = as<arma::mat>(amats(0)).n_elem;
    int Amatrows = as<arma::mat>(amats(0)).n_rows;
    arma::uvec nonzero_elems = find(as<arma::mat>(amats(0)));
    int all_nonzeros = nonzero_elems.n_elem;
    double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
    if (sparse_check <= 0.5 && Amatrows > 100) {
      sparse_switch = 1;
    } else sparse_switch = 0;
    
    if (dense_only) sparse_switch = 0;
    
    int yl = amats.length();
    arma::mat firstmat = as<arma::mat>(amats[0]);
    int matrows = firstmat.n_rows;
    int matcols = firstmat.n_cols;
    
    arma::uvec uniqueyears(yl);
    for (int i = 0; i < yl; i++) {
      uniqueyears(i) = i;
    }
    
    arma::vec twinput;
    if (matrows != matcols) {
      throw Rcpp::exception("Supplied matrices must be square. Please check matrix dimensions.",
        false);
    }
    
    if (tweights.isNotNull()) {
      twinput = as<arma::vec>(tweights);
      if (twinput.n_elem != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    // Now we create the mean matrix
    arma::mat thechosenone(matrows, matcols);
    thechosenone.zeros();
    
    for (int i = 0; i < yl; i++) {
      thechosenone = thechosenone + (as<arma::mat>(amats[i]) / yl);
    }
    
    // Here we take each matrix, run the simulation, and estimate all metrics
    arma::vec startvec;
    int trials {1};
    
    arma::mat slmat(theclairvoyant, trials);
    slmat.zeros();
    arma::vec sl_mean(trials);
    arma::vec sl_var(trials);
    arma::vec sl_sd(trials);
    arma::vec sl_se(trials);
    sl_mean.zeros();
    sl_var.zeros();
    sl_sd.zeros();
    sl_se.zeros();
    
    twinput = twinput / sum(twinput);
    
    arma::uvec thenumbersofthebeast = uniqueyears;
    arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
    startvec = ss3matrix(thechosenone, sparse_switch);
    arma::mat projection = proj3(startvec, amats, theprophecy, 1, 1, 0);
    
    for (int j = 0; j < theclairvoyant; j++) {
      slmat(j,0) = sum(projection.col(j+1));
    }
    
    sl_mean(0) = mean(slmat.col(0));
    sl_var(0) = var(slmat.col(0));
    sl_sd(0) = stddev(slmat.col(0));
    sl_se(0) = sl_sd(0) / sqrt(static_cast<double>(theclairvoyant));
    
    CharacterVector mmpops(1);
    CharacterVector mmpatches(1);
    mmpops(0) = "1";
    mmpatches(0) = "0";
    
    return DataFrame::create(_["pop"] = mmpops, _["patch"] = mmpatches,
      _["a"] = sl_mean, _["var"] = sl_var, _["sd"] = sl_sd, _["se"] = sl_se);
  }
}

//' Estimate Stochastic Sensitivity or Elasticity of Matrix Set
//' 
//' Function \code{stoch_senselas()} estimates the sensitivity and elasticity to
//' matrix elements of \eqn{a}, defined as the long-term arithmetic mean of the
//' log population growth estimated per simulated occasion (as given in equation 2
//' in Tuljapurkar, Horvitz, and Pascarella 2003). 
//' 
//' @name stoch_senselas
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param times Number of occasions to iterate. Defaults to 10,000.
//' @param historical An optional logical value only used if object \code{mpm}
//' is a list of matrices, rather than a \code{lefkoMat} object. Defaults to
//' \code{FALSE} for the former case, and overridden by information supplied in
//' the \code{lefkoMat} object for the latter case.
//' @param style An integer designating whether to estimate sensitivity matrices
//' (\code{1}) or elasticity matrices (\code{2}). Defaults to \code{1}.
//' @param tweights Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' 
//' @return A list of one or two cubes (3d array) where each slice corresponds
//' to a sensitivity or elasticity matrix for a specific pop-patch, followed by
//' the sensitivity or elasticity matrices of all populations (only if multiple
//' pop-patches occur in the input). Two such cubes are only provided when a
//' historical lefkoMat object is used as input, in which case the first
//' element is the historical sensitivity/elasticity matrix, and the second is
//' the ahistorical sensitivity/elasticity matrix.
//' 
//' @section Notes:
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//' 
//' This function currently requires all patches to have the same occasions, if
//' a \code{lefkoMat} object is used as input. Asymmetry in the number of
//' occasions across patches and/or populations will likely cause errors.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.stoch_senselas)]]
Rcpp::List stoch_senselas(List mpm, int times = 10000, bool historical = false,
  int style = 1, Nullable<NumericVector> tweights = R_NilValue) {
  int theclairvoyant = times;
  if (theclairvoyant < 1) {
    throw Rcpp::exception("Option times must be a positive integer.", false);
  }
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = as<List>(mpm["A"]);
    List umats = as<List>(mpm["U"]);
    List fmats = as<List>(mpm["F"]);
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    if (labels.length() < 3) {
      Rf_warningcall(R_NilValue, "This function takes annual matrices as input. This lefkoMat object appears to be a set of mean matrices, and may lack annual matrices.");
    }
    
    // Here we assess ahistorical versions of historical
    // sensitivities and elasticities
    arma::uvec ahstages_id = as<arma::uvec>(stageframe["stage_id"]);
    StringVector ahstages_name = as<StringVector>(stageframe["stage"]);
    int ahstages_num = ahstages_id.n_elem;
    arma::uvec hstages_id2(ahstages_num * ahstages_num);
    hstages_id2.zeros();
    int hstages_num {0};
    
    if (hstages.length() > 1) {
      historical = true;
      arma::uvec hstages_id = as<arma::uvec>(hstages["stage_id_2"]);
      hstages_num = hstages_id.n_elem;
      
      for (int i = 0; i < hstages_num; i++) {
        hstages_id2(i) = hstages_id(i);
      }
      
    } else {
      historical = false;
    }
    
    StringVector poporder = as<StringVector>(labels["pop"]);
    StringVector patchorder = as<StringVector>(labels["patch"]);
    IntegerVector yearorder = as<IntegerVector>(labels["year2"]);
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    IntegerVector uniqueyears = sort_unique(yearorder);
    arma::uvec uniqueyears_arma = as<arma::uvec>(uniqueyears);
    
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    arma::vec twinput;
    if (tweights.isNotNull()) {
      twinput = as<arma::vec>(tweights);
      if (twinput.n_elem != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    arma::vec twinput_corr = twinput / sum(twinput);
    arma::uvec theprophecy_allyears = Rcpp::RcppArmadillo::sample(uniqueyears_arma,
      theclairvoyant, true, twinput_corr);
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize);
    arma::uvec yearsinpatch(loysize);
    patchesinpop.zeros();
    yearsinpatch.zeros();
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = summervacation.n_elem;
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = motorhead.n_elem;
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    // Now we will create a set of means for patches and pops
    List mean_lefkomat;
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1);
    }
    
    // Now we will set up the preliminaries for the stochastic simulations
    List meanamats = as<List>(mean_lefkomat["A"]);
    List mmlabels = as<List>(mean_lefkomat["labels"]);
    StringVector mmpops = as<StringVector>(mmlabels["pop"]);
    StringVector mmpatches = as<StringVector>(mmlabels["patch"]);
    
    int meanmatsize = as<arma::mat>(meanamats[0]).n_elem; // thechosenone
    int meanmatrows = as<arma::mat>(meanamats[0]).n_rows; // thechosenone
    arma::vec startvec(meanmatrows);
    startvec.ones();
    startvec = startvec / meanmatrows; // Start vector for w and v calculations
    int trials = meanamats.length();
    
    // Here we are two cubes to hold sensitivity/elasticity matrices, the first for
    // general use while the second is for ahistorical versions of historical matrices
    arma::cube senscube(meanmatrows, meanmatrows, trials);
    arma::cube senscube_ah(ahstages_num, ahstages_num, trials);
    senscube.zeros();
    senscube_ah.zeros();
    
    // This next matrix will hold the year values for each run
    arma::umat yearspulled(trials, theclairvoyant);
    yearspulled.zeros();
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = allppcs.n_elem;
    arma::uvec year2arma = as<arma::uvec>(yearorder);
    
    // These matrices and vectors will hold R values
    arma::mat Rvecmat(trials, theclairvoyant);
    Rvecmat.zeros();
    
    for (int i= 0; i < allppcsnem; i++) { // This loop goes through each pop-patch
      arma::uvec theprophecy = theprophecy_allyears;
      theprophecy.zeros();
      
      arma::uvec tnotb_patch = find(ppcindex == allppcs(i));
      
      for (int j = 0; j < yl; j++) { // Creates main index marking matrices to use
        // Needs to be modified for situations in which patches do not have the same years
        arma::uvec tnotb_years = find(year2arma == uniqueyears(j));
        arma::uvec thenumbersofthebeast = intersect(tnotb_patch, tnotb_years);
        
        if (thenumbersofthebeast.n_elem > 0) {
          arma::uvec prophetic_yearindices = find(theprophecy_allyears == uniqueyears(j));
          
          if (prophetic_yearindices.n_elem > 0) {
            int replacement = thenumbersofthebeast(0);
            theprophecy.elem(prophetic_yearindices).fill(replacement);
          }
        }
      }
      yearspulled.row(i) = theprophecy.t();
      
      // The next section creates stable stage and rep value vectors arranged in
      // matrix format. The first two are general for whatever has been input,
      // whether historical or ahistorical, while the next two are specifically
      // for ahistorical versions of historical inputs
      arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
      arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
      arma::vec wprojection_ah(ahstages_num);
      arma::vec vprojection_ah(ahstages_num);
      wprojection.zeros();
      vprojection.zeros();
      wprojection_ah.zeros();
      vprojection_ah.zeros();
      
      // Control loop to develop w and v values
      arma::mat crazy_prophet = proj3(startvec, amats, theprophecy, 1, 0, 0);
      wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1),
        theclairvoyant);
      vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0,
        ((startvec.n_elem * 3) - 1), theclairvoyant);
      Rvecmat.row(i) = crazy_prophet.submat((startvec.n_elem * 3), 1,
        (startvec.n_elem * 3), theclairvoyant); // Rvec
      
      // All references should go to senscube, a 3d array to hold sensitivity matrices
      for (int j = 0; j < theclairvoyant; j++) {
        // Main loop for sensitivity matrices, adding each occasion to the
        // respective matrix for each pop-patch
        if (j % 50 == 0) Rcpp::checkUserInterrupt();
        
        arma::vec vtplus1 = vprojection.col(j+1);
        arma::vec wtplus1 = wprojection.col(j+1);
        arma::vec wt = wprojection.col(j);
        arma::mat currentsens_num = vtplus1 * wt.as_row(); // Numerator of key matrix equation
        arma::mat currentsens_den = (Rvecmat(i, j) * vtplus1.as_row() * wtplus1); // Denominator of equation
        double cd_double = currentsens_den(0,0);
        arma::mat currentsens = currentsens_num / (cd_double * theclairvoyant);
        
        // Creates ensitivity matrices
        if (style == 1) {
          senscube.slice(i) += currentsens;
          
          if (historical) {
            wprojection_ah.zeros();
            vprojection_ah.zeros();
            
            // Loop creates ahistorical stable stage dist for projected occasion j+1
            for (int k1 = 0; k1 < hstages_num; k1++) {
              int current_stage2 = hstages_id2(k1);
              wprojection_ah(current_stage2 - 1) = wprojection_ah(current_stage2 - 1)  +
                wtplus1(k1);
            } // k1 loop
            
            // Now the ahistorical reproductive value vector for occasion j+1
            for (int k2 = 0; k2 < hstages_num; k2++) {
              int current_stage2 = hstages_id2(k2);
              
              if (wprojection_ah(current_stage2 - 1) > 0) {
                vprojection_ah(current_stage2 - 1) = vprojection_ah(current_stage2 - 1) +
                  (vtplus1(k2) * wtplus1(k2) / wprojection_ah(current_stage2 - 1));
              }
            } // k2 loop
            
            // Now to propagate the projection sensitivity matrix, and add it to
            // the main sensitivity matrix
            arma::rowvec wtah_tpose = wprojection_ah.as_row();
            arma::rowvec vtah_tpose = vprojection_ah.as_row();
            arma::mat csah_num = vprojection_ah * wtah_tpose;
            arma::mat csah_den = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
            double cdah_double = csah_den(0,0);
            arma::mat csah = csah_num / (cdah_double * theclairvoyant);
            senscube_ah.slice(i) += csah;
          } // if historical statement
        } else {
          // This creates the elasticity matrices
          senscube.slice(i) += currentsens % as<arma::mat>(amats[(theprophecy(j))]);
        }
      }
    }
    
    // This section works on the pop means
    int pop_est {1};
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize);
    arma::uvec yearmatch(loysize);
    popmatch.zeros();
    yearmatch.zeros();
    List meanmatyearlist(yl);
    
    IntegerVector tnotb_all = seq(0, (yl - 1));
    arma::uvec theprophecy = theprophecy_allyears;
    theprophecy.zeros();
    
    for (int j = 0; j < yl; j++) { // Creates main index marking matrices to use
      arma::uvec prophetic_yearindices = find(theprophecy_allyears == uniqueyears(j));
      if (prophetic_yearindices.n_elem > 0) {
        theprophecy.elem(prophetic_yearindices).fill(j);
      }
    }
    
    if (allppcsnem > 1) { // Checks for pop-mean matrices, only occurring with >1 pop-patches
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // This loop goes through each population
        for (int j = 0; j < loysize; j++) { // Checks which A matrices match current pop
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        arma::uvec neededmatspop = find(popmatch == 1);
        
        for (int j = 0; j < yl; j++) { // Checks each year and develops matrix mean across patches
          for (int k = 0; k < loysize; k++) { // Develops a vector to find all matrices for current year
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          
          // This vector catches matrix indices that match the current year and pop
          int crankybankynem = crankybanky.n_elem;
          arma::mat crossmat(meanmatsize, crankybankynem);
          crossmat.zeros();
          
          for (int j = 0; j < crankybankynem; j++) {
            crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
          }
          
          arma::vec happymedium(meanmatsize);
          happymedium.zeros();
          for (int j = 0; j < meanmatsize; j++) {
            for (int k = 0; k < crankybankynem; k++) {
              happymedium(j) = happymedium(j) + crossmat(j, k) / (crankybankynem);
            }
          }
          
          arma::mat finalyearmat = happymedium;
          finalyearmat.reshape(meanmatrows, meanmatrows);
          meanmatyearlist(j) = finalyearmat;
        }
        yearspulled.row(allppcsnem + i) = theprophecy.t();
        
        // Here we use meanmatyearlist in place of amats
        arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
        arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
        arma::vec wprojection_ah(ahstages_num);
        arma::vec vprojection_ah(ahstages_num);
        wprojection.zeros();
        vprojection.zeros();
        wprojection_ah.zeros();
        vprojection_ah.zeros();
        
        // Here we run the control loop to decvelop w and v values
        arma::mat crazy_prophet = proj3(startvec, meanmatyearlist, theprophecy, 1, 0, 0);
        wprojection = crazy_prophet.submat(startvec.n_elem, 0,
          ((startvec.n_elem * 2) - 1), theclairvoyant);
        vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0,
          ((startvec.n_elem * 3) - 1), theclairvoyant);
        Rvecmat.row(allppcsnem + i) = crazy_prophet.submat((startvec.n_elem * 3), 1,
          (startvec.n_elem * 3), theclairvoyant); // Rvec
        
        // All references should go to senscube, a 3d array holding sensitivity matrices
        
        // Next is the main time loop for the sensitivity matrices, adding each
        // occasion to the respective matrix for each pop-patch
        for (int j = 0; j < theclairvoyant; j++) {  
          arma::vec vtplus1 = vprojection.col(j+1);
          arma::vec wtplus1 = wprojection.col(j+1);
          arma::vec wt = wprojection.col(j);
          
          arma::mat currentsens_num = vtplus1 * wt.as_row(); // Numerator of key matrix equation
          arma::mat currentsens_den = (Rvecmat((allppcsnem + i), j) *
            vtplus1.as_row() * wtplus1); // Denominator of key matrix equation
          double cd_double = currentsens_den(0,0);
          arma::mat currentsens = currentsens_num / (cd_double * theclairvoyant);
          
          if (style == 1) {
            // Sensitivity matrix
            senscube.slice(allppcsnem + i) += currentsens; 
            
            if (historical) {
              wprojection_ah.zeros();
              vprojection_ah.zeros();
              
              // Creates ahistorical stable stage distribution for projected occasion j+1
              for (int k1 = 0; k1 < hstages_num; k1++) {
                int current_stage2 = hstages_id2(k1);
                wprojection_ah(current_stage2 - 1) = wprojection_ah(current_stage2 - 1)  +
                  wtplus1(k1);
              } // k1 loop
              
              // Now the ahistorical reproductive value vector for occasion j+1
              for (int k2 = 0; k2 < hstages_num; k2++) {
                int current_stage2 = hstages_id2(k2);
                
                if (wprojection_ah(current_stage2 - 1) > 0) {
                  vprojection_ah(current_stage2 - 1) = vprojection_ah(current_stage2 - 1) +
                    (vtplus1(k2) * wtplus1(k2) / wprojection_ah(current_stage2 - 1));
                }
              } // k2 loop
              
              // Here we propagate the projection sensitivity matrix, and add it to
              // the main sensitivity matrix
              arma::rowvec wtah_tpose = wprojection_ah.as_row();
              arma::rowvec vtah_tpose = vprojection_ah.as_row();
              arma::mat csah_num = vprojection_ah * wtah_tpose;
              arma::mat csah_den = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
              double cdah_double = csah_den(0,0);
              arma::mat csah = csah_num / (cdah_double * theclairvoyant);
              senscube_ah.slice(i) += csah;
              
            } // if historical statement
          } else {
            // This is the elasticity matrix
            senscube.slice(allppcsnem + i) += currentsens % as<arma::mat>(meanmatyearlist[(theprophecy(j))]); 
          }
        }
      } // for loop i, for populations
    } // if statement, checking if more than one patch and determining if pop means need to be dealt with
    
    if (historical && style == 2) {
      for (int k = 0; k < trials; k++) {
        arma::uvec hstages_id1 = as<arma::uvec>(hstages["stage_id_1"]);
        arma::uvec hstages_id2 = as<arma::uvec>(hstages["stage_id_2"]);
        arma::mat elasah(ahstages_num, ahstages_num);
        elasah.zeros();
        
        arma::mat hslice = senscube.slice(k);
        for (int i = 0; i < hstages_num; i++) {
          for (int j = 0; j < hstages_num; j++) {
            elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) =
              elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) + hslice(j, i);
          }
        }
        senscube_ah.slice(k) = elasah;
      }
    }
    
    return Rcpp::List::create(_["maincube"] = senscube, _["ahcube"] = senscube_ah);
    
  } else { 
    // Single list of A matrices as input
    List amats = mpm;
    
    int yl = amats.length();
    arma::mat firstmat = as<arma::mat>(amats[0]);
    int matrows = firstmat.n_rows;
    int matcols = firstmat.n_cols;
    
    IntegerVector uniqueyears = seq(0, (yl - 1));
    arma::uvec uniqueyears_arma = as<arma::uvec>(uniqueyears);
    arma::vec twinput;
    
    if (matrows != matcols) {
      throw Rcpp::exception("Supplied matrices must be square. Please check matrix dimensions.",
        false);
    }
    
    if (tweights.isNotNull()) {
      twinput = as<arma::vec>(tweights);
      if (twinput.n_elem != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      
    } else {
      twinput.resize(yl);
      twinput.ones();
    }
    
    // Here we set up the vector of chosen occasions, sampled from all possible occasions
    arma::vec twinput_corr = twinput / sum(twinput);
    arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears_arma, theclairvoyant,
      true, twinput_corr);
    
    // Here we initialize an empty matrix and start vector for w and v
    // The matrix will be updated at each occasion
    arma::vec startvec(matrows);
    startvec.ones();
    startvec = startvec / matrows; // The is the start vector for w and v calculations
    int trials {1};
    
    // Here we initialize a flat cube to hold the sensitivity or elasticity matrix
    arma::cube senscube(matrows, matrows, trials);
    senscube.zeros();
    
    // These matrices and vectors will hold R values
    arma::mat Rvecmat(trials, theclairvoyant);
    Rvecmat.zeros();
    
    arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
    wprojection.zeros();
    vprojection.zeros();
    
    // Here we run the control loop to develop w and v values
    arma::mat crazy_prophet = proj3(startvec, amats, theprophecy, 1, 0, 0);
    wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1),
      theclairvoyant);
    vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0,
      ((startvec.n_elem * 3) - 1), theclairvoyant);
    Rvecmat.row(0) = crazy_prophet.submat((startvec.n_elem * 3), 1,
      (startvec.n_elem * 3), theclairvoyant); // Rvec
    
    // All references should go to senscube, a 3d array holding sensitivity matrices
    for (int j = 0; j < theclairvoyant; j++) { // Main loop for sensitivity matrices
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // Numerator of key matrix equation
      arma::mat currentsens_den = (Rvecmat(0, j) * vtplus1.as_row() * 
        wprojection.col(j+1)); // Denominator of key matrix equation
      double cd_double = currentsens_den(0,0);
      arma::mat currentsens = currentsens_num / (cd_double * theclairvoyant);
      
      if (style == 1) {
        senscube.slice(0) += currentsens; // Sensitivity matrix
      } else {
        senscube.slice(0) += currentsens % as<arma::mat>(amats[(theprophecy(j))]); // Elasticity matrix
      }
    }
    
  return Rcpp::List::create(_["maincube"] = senscube);    
  }
}

//' Creates Size Index for Elasticity Summaries of hMPMs
//' 
//' Function \code{bambi3()} creates an index of estimable elements in
//' historical matrices, and details the kind of transition that it is.
//' 
//' @name bambi3
//' 
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' 
//' @return A data frame with the following elements:
//' \item{index}{Vector index of matrix element in C++ terms.}
//' \item{transition}{Category of transition.}
//' \item{size3}{Size in occasion \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in occasion \emph{t}+1.}
//' \item{entrystatus3}{Entry status in occasion \emph{t}+1.}
//' \item{size2}{Size in occasion \emph{t}.}
//' \item{repstatus2}{Reproductive status in occasion \emph{t}.}
//' \item{entrystatus2}{Entry status in occasion \emph{t}.}
//' \item{size1}{Size in occasion \emph{t}-1.}
//' \item{repstatus1}{Reproductive status in occasion \emph{t}11.}
//' \item{entrystatus1}{Entry status in occasion \emph{t}-1.}
//'
//' The kind of transitions conforms to the following code: \code{10}: full
//' stasis, \code{11}: stasis to growth, \code{12}: full growth, \code{13}:
//' growth to stasis, \code{14}: stasis to shrinkage, \code{15}: full shrinkage,
//' \code{16}: shrinkage to stasis, \code{17}: growth to shrinkage, \code{18}:
//' shrinkage to growth, \code{20}: stasis to fecundity, \code{21}: growth to
//' fecundity, \code{22}: shrinkage to fecundity, \code{23}: fecundity to
//' stasis, \code{24}: fecundity to growth, \code{25}: fecundity to shrinkage,
//' \code{26}: fecundity to fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.bambi3)]]
DataFrame bambi3(DataFrame stages, DataFrame hstages) {
  StringVector stagenames = as<StringVector>(stages["stage"]);
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  arma::vec sizes = as<arma::vec>(stages["original_size"]);
  arma::uvec repstatus = as<arma::uvec>(stages["repstatus"]);
  arma::uvec entrystage = as<arma::uvec>(stages["entrystage"]);
  int numstages = astages.n_elem;
  
  arma::uvec hstage3in = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec hstage2nin = as<arma::uvec>(hstages["stage_id_1"]);
  int numhstages = hstage3in.n_elem;
  
  hstage3in = hstage3in - 1;
  hstage2nin = hstage2nin - 1;
  int predictedsize = numstages * numstages * numstages;
  
  arma::ivec hsindexl(predictedsize);
  arma::uvec transition_type(predictedsize);
  hsindexl.fill(-1);
  transition_type.zeros();
  
  arma::vec size1(predictedsize);
  arma::vec size2(predictedsize);
  arma::vec size3(predictedsize);
  size1.fill(-1);
  size2.fill(-1);
  size3.fill(-1);
  
  arma::uvec repstatus1(predictedsize);
  arma::uvec repstatus2(predictedsize);
  arma::uvec repstatus3(predictedsize);
  repstatus1.zeros();
  repstatus2.zeros();
  repstatus3.zeros();
  
  arma::uvec entrystatus1(predictedsize);
  arma::uvec entrystatus2(predictedsize);
  arma::uvec entrystatus3(predictedsize);
  entrystatus1.zeros();
  entrystatus2.zeros();
  entrystatus3.zeros();
  
  StringVector longnames3(predictedsize);
  StringVector longnames2(predictedsize);
  StringVector longnames1(predictedsize);
  
  int counter = 0;
  
  for (int i1 = 0; i1 < numhstages; i1++) {
    for (int i2 = 0; i2 < numhstages; i2++) {
      if (hstage3in(i1) == (hstage2nin(i2))) {
        hsindexl(counter) = (i1 * numhstages) + i2;
        
        int stage1 = hstage2nin(i2);
        longnames1(counter) = stagenames(stage1);
        size1(counter) = sizes(stage1);
        repstatus1(counter) = repstatus(stage1);
        entrystatus1(counter) = entrystage(stage1);
        
        int stage2 = hstage2nin(i1);
        longnames2(counter) = stagenames(stage2);
        size2(counter) = sizes(stage2);
        repstatus2(counter) = repstatus(stage2);
        entrystatus2(counter) = entrystage(stage2);
        
        int stage3 = hstage3in(i2);
        longnames3(counter) = stagenames(stage3);
        size3(counter) = sizes(stage3);
        repstatus3(counter) = repstatus(stage3);
        entrystatus3(counter) = entrystage(stage3);
        
        if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
          if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
            transition_type(counter) = 26; // Fecundity to fecundity
          } else if (size2(counter) == size1(counter)) {
            if (repstatus2(counter) > repstatus1(counter) || entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 21; // Growth to fecundity
            } else if (repstatus2(counter) < repstatus1(counter) ||
              entrystatus2(counter) > entrystatus1(counter)) {
              transition_type(counter) = 22; // Shrinkage to fecundity
            } else {
              transition_type(counter) = 20; // Stasis to fecundity
            }
          } else if (size2(counter) > size1(counter)) {
            transition_type(counter) = 21; // Growth to fecundity
          } else if (size2(counter) < size1(counter)) {
            transition_type(counter) = 22; // Shrinkage to fecundity
          }
        } else if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
          if (size3(counter) == size2(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 24; // Fecundity to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 25; // Fecundity to shrinkage
            } else {
              transition_type(counter) = 23; // Fecundity to stasis
            }
          } else if (size3(counter) > size2(counter)) {
            transition_type(counter) = 24; // Fecundity to growth
          } else if (size3(counter) < size2(counter)) {
            transition_type(counter) = 25; // Fecundity to shrinkage
          }
        } else if (size3(counter) == size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 12; // Full growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 17; // Growth to shrinkage
            } else {
              transition_type(counter) = 13; // Growth to stasis
            }
          } else if (repstatus2(counter) < repstatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 18; // Shrinkage to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 15; // Full shrinkage
            } else {
              transition_type(counter) = 16; // Shrinkage to stasis
            }
          } else {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 11; // Stasis to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 14; // Stasis to shrinkage
            } else {
              transition_type(counter) = 10; // Full stasis
            }
          }
        } else if (size3(counter) > size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus2(counter) < repstatus1(counter) ||
            entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else {
            transition_type(counter) = 11; // Stasis to growth
          }
        } else if (size3(counter) > size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 12; // Full growth
        } else if (size3(counter) == size2(counter) && size2(counter) > size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) ||
            entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus3(counter) < repstatus2(counter) ||
            entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else {
            transition_type(counter) = 13; // Growth to stasis
          }
        } else if (size3(counter) < size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else if (repstatus2(counter) < repstatus1(counter) ||
            entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 14; // Stasis to shrinkage
          }
        } else if (size3(counter) < size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 15; // Full shrinkage
        } else if (size3(counter) == size2(counter) && size2(counter) < size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) ||
            entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else if (repstatus3(counter) < repstatus2(counter) ||
            entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 16; // Shrinkage to stasis
          }
        } else if (size3(counter) < size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 17; // Growth to shrinkage
        } else if (size3(counter) > size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 18; // Shrinkage to growth
        }
        counter++;
      }
    }
  }
  
  StringVector names3(counter);
  StringVector names2(counter);
  StringVector names1(counter);
  for (int i = 0; i < counter; i++) {
    names3(i) = longnames3(i);
    names2(i) = longnames2(i);
    names1(i) = longnames1(i);
  }
  
  arma::uvec targetindices = find(hsindexl > -1);
  arma::ivec hsindex = hsindexl.elem(targetindices);
  arma::uvec t_type = transition_type.elem(targetindices);
  arma::vec size3c = size3.elem(targetindices);
  arma::vec size2c = size2.elem(targetindices);
  arma::vec size1c = size1.elem(targetindices);
  arma::uvec r_status3 = repstatus3.elem(targetindices);
  arma::uvec r_status2 = repstatus2.elem(targetindices);
  arma::uvec r_status1 = repstatus1.elem(targetindices);
  arma::uvec e_status3 = entrystatus3.elem(targetindices);
  arma::uvec e_status2 = entrystatus2.elem(targetindices);
  arma::uvec e_status1 = entrystatus1.elem(targetindices);
  
  DataFrame output = DataFrame::create(Named("index") = hsindex, _["transition"] = t_type,
    _["stage3"] = names3, _["size3"] = size3c, _["repstatus3"] = r_status3, _["entrystatus3"] = e_status3,
    _["stage2"] = names2, _["size2"] = size2c, _["repstatus2"] = r_status2, _["entrystatus2"] = e_status2,
    _["stage1"] = names1, _["size1"] = size1c, _["repstatus1"] = r_status1, _["entrystatus1"] = e_status1);
  
  return output;
}

//' Creates Size Index for Elasticity Summaries of ahMPMs
//' 
//' Function \code{bambi2()} creates an index of estimable elements in
//' ahistorical matrices, and details the kind of transition that it is.
//' 
//' @name bambi2
//' 
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' 
//' @return A data frame with the following elements:
//' \item{index}{Vector index of matrix element in C++ terms.}
//' \item{transition}{Category of transition.}
//' \item{stage3}{Stage in occasion \emph{t}+1.}
//' \item{size3}{Size in occasion \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in occasion \emph{t}+1.}
//' \item{entrystatus3}{Entry status in occasion \emph{t}+1.}
//' \item{stage2}{Stage in occasion \emph{t}.}
//' \item{size2}{Size in occasion \emph{t}.}
//' \item{repstatus2}{Reproductive status in occasion \emph{t}.}
//' \item{entrystatus2}{Entry status in occasion \emph{t}.}
//'
//' The kind of transitions conforms to the following code: \code{1}: stasis, 
//' \code{2}: growth, \code{3}: shrinkage, \code{4}: fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.bambi2)]]
DataFrame bambi2(DataFrame stages) {
  StringVector stagenames = as<StringVector>(stages["stage"]);
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  arma::vec sizes = as<arma::vec>(stages["original_size"]);
  arma::uvec repstatus = as<arma::uvec>(stages["repstatus"]);
  arma::uvec entrystage = as<arma::uvec>(stages["entrystage"]);
  int numstages = astages.n_elem;
  astages = astages - 1;
  int predictedsize = numstages * numstages;
  
  arma::ivec ahsindexl(predictedsize);
  arma::uvec transition_type(predictedsize);
  ahsindexl.fill(-1);
  transition_type.zeros();
  
  StringVector longstages3(predictedsize);
  StringVector longstages2(predictedsize);
  
  arma::vec size2(predictedsize);
  arma::vec size3(predictedsize);
  size2.fill(-1);
  size3.fill(-1);
  
  arma::uvec repstatus2(predictedsize);
  arma::uvec repstatus3(predictedsize);
  repstatus2.zeros();
  repstatus3.zeros();
  
  arma::uvec entrystatus2(predictedsize);
  arma::uvec entrystatus3(predictedsize);
  entrystatus2.zeros();
  entrystatus3.zeros();
  
  int counter = 0;
  
  for (int i1 = 0; i1 < numstages; i1++) {
    for (int i2 = 0; i2 < numstages; i2++) {
      ahsindexl(counter) = (i1 * numstages) + i2;
      
      int stage2 =  astages(i1);
      longstages2(counter) = stagenames(stage2);
      size2(counter) = sizes(stage2);
      repstatus2(counter) = repstatus(stage2);
      entrystatus2(counter) = entrystage(stage2);
      
      int stage3 = astages(i2);
      longstages3(counter) = stagenames(stage3);
      size3(counter) = sizes(stage3);
      repstatus3(counter) = repstatus(stage3);
      entrystatus3(counter) = entrystage(stage3);
      
      if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
        transition_type(counter) = 4; // Fecundity
      } else if (size3(counter) == size2(counter)) {
        if (repstatus3(counter) > repstatus2(counter) ||
          entrystatus3(counter) < entrystatus2(counter)) {
          transition_type(counter) = 2; // Growth
        } else if (repstatus3(counter) < repstatus2(counter) ||
          entrystatus3(counter) > entrystatus2(counter)) {
          transition_type(counter) = 3; // Shrinkage
        } else {
          transition_type(counter) = 1; // Stasis
        }
      } else if (size3(counter) > size2(counter)) {
        transition_type(counter) = 2; // Growth
      } else if (size3(counter) < size2(counter)) {
        transition_type(counter) = 3; // Shrinkage
      }
      
      counter++;
    }
  }
  
  arma::uvec targetindices = find(ahsindexl > -1);
  arma::ivec ahsindex = ahsindexl.elem(targetindices);
  
  DataFrame output = DataFrame::create(Named("index") = ahsindex,
    _["transition"] = transition_type, _["stage3"] = longstages3, _["size3"] = size3,
    _["repstatus3"] = repstatus3, _["entrystatus3"] = entrystatus3, _["stage2"] = longstages2,
    _["size2"] = size2, _["repstatus2"] = repstatus2, _["entrystatus2"] = entrystatus2);
  
  return output;
}

//' Creates Summary Data for Elasticity Matrix Inputs
//' 
//' Function \code{demolition3()} sums elasticity values from elasticity
//' matrices, and LTRE contributions from LTRE and sLTRE matrices, according to
//' the categories developed by functions \code{bambi2()} and \code{bambi3()}.
//' 
//' @name demolition3
//' 
//' @param e_amat A single elasticity, LTRE, or sLTRE matrix.
//' @param bambesque This is the output from \code{bambi2()} or \code{bambi3()}
//' corresponding to the current lefkoMat object. The format is a data frame
//' giving the indices and characteristics of all predicted potential non-zero
//' elements in the supplied matrix.
//' @param amat_ The A matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions.
//' @param fmat_ The F matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions.
//' 
//' @return A list with two data frames, one showing the summed elasticities for
//' the historical matrix supplied (if supplied), and the other showing the
//' ahistorical summary of the historical matrix or the summed elasticities of
//' a supplied ahistorical elasticity matrix. Also includes sums of only the
//' positive elements and only the negative elements, in all cases.
//' 
//' @section Notes:
//' If the original matrices are provided, then this function was made to split
//' co-occurring survival-fecundity elasticities according to the ratio of the
//' fecundity portion of the element to the survival portion of that element.
//' However, this transition splitting capability developed using the original
//' matrices does not currently work properly, and so it is better to use this
//' function without objects \code{amat_} and \code{fmat_}, forcing co-occurring
//' survival-fecundity transitions to be treated as fecundity only.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.demolition3)]]
List demolition3(arma::mat e_amat, DataFrame bambesque,
  Nullable<Rcpp::NumericVector> amat_ = R_NilValue,
  Nullable<Rcpp::NumericVector> fmat_ = R_NilValue) {
  arma::uvec eindices = as<arma::uvec>(bambesque["index"]);
  arma::uvec categories = as<arma::uvec>(bambesque["transition"]);
  
  int e_amatsize = e_amat.n_elem;
  int e_amatrows = e_amat.n_rows;
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  arma::mat amat(e_amatrows, e_amatrows);
  arma::mat fmat(e_amatrows, e_amatrows);
  
  if (maxelem > e_amatsize) {
    throw Rcpp::exception("Supplied info does not correspond to input matrices.",
      false);
  }
  
  if (amat_.isNotNull() && fmat_.isNotNull()) {
    amat = Rcpp::as<arma::mat>(amat_);
    fmat = Rcpp::as<arma::mat>(fmat_);
    
  } else {
    amat.ones();
    fmat.zeros();
    
    arma::uvec fec_trans = find(categories == 4);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 20);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 21);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 22);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 26);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
  }
  
  // This splits fecundity transitions if they include survival portions
  arma::mat corr_mat = amat;
  arma::uvec z_indices = find(corr_mat == 0);
  int z_indicesnem = z_indices.n_elem;
  for (int i = 0; i < z_indicesnem; i++) {
    corr_mat(z_indices(i)) = 1.0;
  }
  
  arma::mat fec_fraction = fmat / corr_mat;
  DataFrame histout;
  DataFrame ahistout;
  
  StringVector histcats {"Full stasis", "Stasis to growth", "Full growth",
    "Growth to stasis", "Stasis to shrinkage", "Full shrinkage",
    "Shrinkage to stasis", "Growth to shrinkage", "Shrinkage to growth",
    "Stasis to fecundity", "Growth to fecundity", "Shrinkage to fecundity",
    "Fecundity to stasis", "Fecundity to growth", "Fecundity to shrinkage",
    "Fecundity to fecundity"};
  arma::uvec histcatnums {10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23,
    24, 25, 26};
  arma::vec histsums(16);
  arma::vec histpos(16);
  arma::vec histneg(16);
  arma::vec hc_ahistsums(4);
  arma::vec hc_ahistpos(4);
  arma::vec hc_ahistneg(4);
  histsums.zeros();
  histpos.zeros();
  histneg.zeros();
  hc_ahistsums.zeros();
  hc_ahistpos.zeros();
  hc_ahistneg.zeros();
  
  StringVector ahistcats {"Stasis", "Growth", "Shrinkage", "Fecundity"};
  arma::uvec ahistcatnums {1, 2, 3, 4};
  arma::vec ahistsums(4);
  arma::vec ahistpos(4);
  arma::vec ahistneg(4);
  ahistsums.zeros();
  ahistpos.zeros();
  ahistneg.zeros();
  
  if (minindex > 9) {
    // Object minindex will only be above 9 if the supplied object is historical
    arma::vec size3 = as<arma::vec>(bambesque["size3"]);
    arma::vec size2 = as<arma::vec>(bambesque["size2"]);
    arma::vec size1 = as<arma::vec>(bambesque["size1"]);
    
    arma::uvec repstatus3 = as<arma::uvec>(bambesque["repstatus3"]);
    arma::uvec repstatus2 = as<arma::uvec>(bambesque["repstatus2"]);
    arma::uvec repstatus1 = as<arma::uvec>(bambesque["repstatus1"]);
    
    arma::uvec entrystatus3 = as<arma::uvec>(bambesque["entrystatus3"]);
    arma::uvec entrystatus2 = as<arma::uvec>(bambesque["entrystatus2"]);
    arma::uvec entrystatus1 = as<arma::uvec>(bambesque["entrystatus1"]);
    
    for (int i = 0; i < 16; i++) {
      arma::uvec currentguys = find(categories == histcatnums(i));
      int currentguysnem = currentguys.n_elem;
      
      if (histcatnums(i) == 20 || histcatnums(i) == 21 || histcatnums(i) == 22 ||
        histcatnums(i) == 26) { // Fecundity transitions
      
        // Splits transitions that are actually combos of fecundity and survival
        for (int j = 0; j < currentguysnem; j++) {
          int this_guy = eindices(currentguys(j));
          
          if (fec_fraction(this_guy) == 1) {
            hc_ahistsums(3) += (e_amat(this_guy));
            histsums(i) += (e_amat(this_guy));
            
            if (e_amat(this_guy) > 0) {
              hc_ahistpos(3) += (e_amat(this_guy));
              histpos(i) += (e_amat(this_guy));
            } else if (e_amat(this_guy) < 0) {
              hc_ahistneg(3) += (e_amat(this_guy));
              histneg(i) += (e_amat(this_guy));
            }
            
          } else {
            hc_ahistsums(3) += (e_amat(this_guy) * fec_fraction(this_guy));
            histsums(i) += (e_amat(this_guy) * fec_fraction(this_guy));
            
            if (e_amat(this_guy) > 0) {
              hc_ahistpos(3) += (e_amat(this_guy) * (fec_fraction(this_guy)));
              histpos(i) += (e_amat(this_guy) * (fec_fraction(this_guy)));
            } else if (e_amat(this_guy) < 0) {
              hc_ahistneg(3) += (e_amat(this_guy) * (fec_fraction(this_guy)));
              histneg(i) += (e_amat(this_guy) * (fec_fraction(this_guy)));
            }
            
            arma::uvec counter = find(eindices == this_guy);
            if (entrystatus2(counter(0)) == 1 && repstatus1(counter(0)) == 1) {
              if (size3(counter(0)) == size2(counter(0)) &&
                repstatus3(counter(0)) == repstatus2(counter(0)) &&
                  entrystatus3(counter(0)) == entrystatus2(counter(0))) {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (size3(counter(0)) > size2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (size3(counter(0)) < size2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
              
            } else if (size3(counter(0)) == size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
              if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
                
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
                
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
                
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                  
                } else {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
                
              } else {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
              
            } else if (size3(counter(0)) > size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
              if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
              
            } else if (size3(counter(0)) > size2(counter(0)) &&
              size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
              
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
              if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
              
            } else if (size3(counter(0)) == size2(counter(0)) &&
              size2(counter(0)) < size1(counter(0))) {
              if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
                
              } else {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(6) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(86) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(6) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            }
          }
        }
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 ||
        histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0);
        arma::uvec all_es_neg = find(all_es < 0);
        int all_es_pos_num = all_es_pos.n_elem;
        int all_es_neg_num = all_es_neg.n_elem;
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(2) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(2) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(2) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 ||
        histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0);
        arma::uvec all_es_neg = find(all_es < 0);
        int all_es_pos_num = all_es_pos.n_elem;
        int all_es_neg_num = all_es_neg.n_elem;
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(0) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(0) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(0) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 ||
        histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0);
        arma::uvec all_es_neg = find(all_es < 0);
        int all_es_pos_num = all_es_pos.n_elem;
        int all_es_neg_num = all_es_neg.n_elem;
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(1) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(1) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(1) += sum(all_es.elem(all_es_neg));
        }
      }
    }
    
    histout = DataFrame::create(Named("category") = histcats, _["elas"] = histsums,
      _["elas_pos"] = histpos, _["elas_neg"] = histneg);
    ahistout = DataFrame::create(Named("category") = ahistcats, _["elas"] = hc_ahistsums,
      _["elas_pos"] = hc_ahistpos, _["elas_neg"] = hc_ahistneg);
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      arma::vec all_es = e_amat.elem(eindices(currentguys));
      arma::uvec all_es_pos = find(all_es > 0);
      arma::uvec all_es_neg = find(all_es < 0);
      int all_es_pos_num = all_es_pos.n_elem;
      int all_es_neg_num = all_es_neg.n_elem;
      
      double getoutofdodge = sum(all_es);
      ahistsums(i) += getoutofdodge;
      if (all_es_pos_num > 0) {
        ahistpos(i) += sum(all_es.elem(all_es_pos));
      }
      if (all_es_neg_num > 0) {
        ahistneg(i) += sum(all_es.elem(all_es_neg));
      }
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats, _["elas"] = ahistsums,
      _["elas_pos"] = ahistpos, _["elas_neg"] = ahistneg);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}

//' Estimate LTRE of Any Population Matrix
//' 
//' \code{ltre3matrix()} returns the one-way fixed deterministic LTRE matrix of
//' a dense or sparse set of input matrices.
//' 
//' @name ltre3matrix
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from_ \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param mean A logical value indicating whether to use the element-wise mean
//' matrix as the reference.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a cube of LTRE contributions, with each slice
//' corresponding to each input matrix in Amats. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ltre3matrix)]]
arma::cube ltre3matrix(List Amats, Rcpp::IntegerVector refnum,
  Nullable<Rcpp::List> refmats_ = R_NilValue, bool mean = true,
  bool sparse = false) {
  int Amatnum = Amats.length();
  int matdim = as<arma::mat>(Amats(0)).n_rows;
  
  List eigenstuff;
  List mean_set;
  arma::cube senscube(matdim, matdim, Amatnum);
  senscube.zeros();
  
  if (!sparse) {
    // Dense matrix analysis
    arma::mat finalrefmat(matdim, matdim);
    finalrefmat.zeros();
    int refmatnum = refnum.length();
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::mat>(refmats(refnum(i))) / static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = as<arma::mat>(refmats(refnum(0)));
      }
      
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::mat>(Amats(refnum(i))) /
            static_cast<double>(Amatnum));
        }
      } else {
        finalrefmat = as<arma::mat>(Amats(refnum(0)));
      }
    }
    
    // Now we create the halfway matrices and run the sensitivities
    arma::mat halfmat(matdim, matdim);
    arma::mat diffmat(matdim, matdim);
    halfmat.zeros();
    diffmat.zeros();
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (as<arma::mat>(Amats(i)))) / static_cast<double>(2.0);
      diffmat = (as<arma::mat>(Amats(i))) - finalrefmat;
      
      senscube.slice(i) = diffmat % sens3matrix(halfmat, 0);
    }
    
    return senscube;
    
  } else {
    // Sparse matrix analysis
    arma::sp_mat finalrefmat(matdim, matdim);
    finalrefmat.zeros();
      
    int refmatnum = refnum.length();
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) / 
            static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = (arma::sp_mat(as<arma::mat>(refmats(refnum(0)))));
      }
      
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (arma::sp_mat(as<arma::mat>(Amats(refnum(i)))) / 
            static_cast<double>(Amatnum));
        }
        
      } else {
        finalrefmat = (arma::sp_mat(as<arma::mat>(Amats(refnum(0)))));
      }
    }
    
    // Now we create the halfway matrices and run the sensitivities
    arma::sp_mat halfmat(matdim, matdim);
    arma::sp_mat diffmat(matdim, matdim);
    halfmat.zeros();
    diffmat.zeros();
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (arma::sp_mat(as<arma::mat>(Amats(i))))) /
        static_cast<double>(2.0);
      diffmat = (arma::sp_mat(as<arma::mat>(Amats(i)))) - finalrefmat;
      
      senscube.slice(i) = arma::mat(diffmat % sens3sp_matrix(halfmat, diffmat));
    }
    
    return senscube;
  }
}

//' Estimate sLTRE of Any Population Matrix
//' 
//' \code{sltre3matrix()} returns the one-way stochastic LTRE matrix of
//' a dense or sparse set of input matrices.
//' 
//' @name sltre3matrix
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object).
//' @param loy A data frame showing the order of populations, patches, and
//' occasions of the matrices provided in object \code{Amats}.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param tweights_ Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' @param steps The number of occasions to project the stochastic simulation
//' forward, if performing an sLTRE. Defaults to \code{10000}. Note that the
//' total number of occasions projected equals this number plus the number given
//' in object \code{burnin}.
//' @param burnin The number of initial occasions to project the population
//' without calculating population metrics. Defaults to 3000.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a list of two lists of matrices. The first,
//' \code{cont_mean}, holds the sLTRE contributions of shifts in mean elements.
//' The second, \code{cont_sd}, holds the sLTRE contributions of shifts in
//' temporal standard deviations of matrix elements.
//' 
//' @section Notes:
//' This function uses the simulation approach developed in Davison et al.
//' (2010), which is a good approximation though not an analytical solution.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sltre3matrix)]]
Rcpp::List sltre3matrix(List Amats, DataFrame loy, Rcpp::IntegerVector refnum,
  Nullable<Rcpp::List> refmats_ = R_NilValue,
  Nullable<arma::vec> tweights_ = R_NilValue, int steps = 10000,
  int burnin = 3000, bool sparse = false) {
  int theclairvoyant = steps + burnin;
  arma::vec tweights;
  int Amatnum = Amats.length();
  int matdim = as<arma::mat>(Amats(0)).n_rows;
  int matlength = as<arma::mat>(Amats(0)).n_elem;
  int refmatnum = refnum.length();
  
  List eigenstuff;
  List mean_set;
  List poppatch_meanmat;
  List poppatch_sdmat;
  List ref_byyear;
  List elas_mean;
  List elas_sd;
  List diff_meanmat;
  List diff_sdmat;
  List cont_meanmat;
  List cont_sdmat;
  
  // This section creates the order vectors for pop-patch-year
  StringVector pops = as<StringVector>(loy["pop"]);
  arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
  StringVector patches = as<StringVector>(loy["patch"]);
  arma::uvec year2c = as<arma::uvec>(loy["year2c"]);
  arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
  arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
  arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  arma::uvec uniqueyears = unique(year2c);
  int numpoppatches = uniquepoppatches.n_elem;
  int numyears = uniqueyears.n_elem;
  
  arma::vec twinput_corr;
  arma::uvec theprophecy;
  
  // The main loop
  if (!sparse) {
    // Dense matrix analysis
    arma::mat mat_mean(matdim, matdim);
    arma::mat mat_sd(matdim, matdim);
    arma::mat ref_matmean(matdim, matdim);
    arma::mat ref_matsd(matdim, matdim);
    mat_mean.zeros();
    mat_sd.zeros();
    ref_matmean.zeros();
    ref_matsd.zeros();
    
    // First the pop/patch means and sds
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (as<arma::mat>(Amats(poppatch_chosen(j))) / 
          static_cast<double>(numpoppatch_chosen));
      }
      if (i == 0) {
        poppatch_meanmat = List::create(mat_mean);
      } else {
        poppatch_meanmat.push_back(mat_mean);
      }
      mat_mean.zeros();
    }
    
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      arma::vec mat_elems(numpoppatch_chosen);
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          mat_elems(k) = as<arma::mat>(Amats(poppatch_chosen(k)))(j);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(j) = stddev(mat_elems, 1);
        }
      }
      if (i == 0) {
        poppatch_sdmat = List::create(mat_sd);
      } else {
        poppatch_sdmat.push_back(mat_sd);
      }
      mat_sd.zeros();
    }
    
    // Here we create the reference matrix sets
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      mat_mean.zeros();
      mat_sd.zeros();
      ref_byyear = List::create(refmats);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (as<arma::mat>(refmats(refnum(i))) / static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum);
        mat_elems.zeros();
        
        for (int j = 0; j < refmatnum; j++) {
          mat_elems(j) = as<arma::mat>(refmats(refnum(j)))(i);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(i) = stddev(mat_elems, 1);
        }
      }
      ref_matmean = mat_mean;
      ref_matsd = mat_sd;
      
      mat_mean.zeros();
      mat_sd.zeros();
      
    } else {
      if (refmatnum == Amatnum) {
        
        // Reference by year
        for (int i = 0; i < numyears; i++) {
          arma::uvec year_chosen = find(year2c == uniqueyears(i));
          int numyear_chosen = year_chosen.n_elem;
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (as<arma::mat>(Amats(year_chosen(j))) /
              static_cast<double>(numyear_chosen));
          }
          if (i == 0) {
            ref_byyear = List::create(mat_mean);
          } else {
            ref_byyear.push_back(mat_mean);
          }
          mat_mean.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < numyears; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) / static_cast<double>(numyears));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < numyears; j++) {
            mat_elems(j) = as<arma::mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        mat_mean.zeros();
        mat_sd.zeros();
        
      } else {
        // Reference by year
        for (int i = 0; i < refmatnum; i++) {
          if (i == 0) {
            ref_byyear = List::create(as<arma::mat>(Amats(i)));
          } else {
            ref_byyear.push_back(as<arma::mat>(Amats(i)));
          }
          mat_mean.zeros();
          mat_sd.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) /
            static_cast<double>(refmatnum));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < refmatnum; j++) {
            mat_elems(j) = as<arma::mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        mat_mean.zeros();
        mat_sd.zeros();
      }
    }
    
    // Time weights
    if (tweights_.isNotNull()) {
      tweights = as<arma::vec>(tweights_);
      if (tweights.n_elem != numyears) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions in the reference matrix set.", false);
      }
    } else {
      tweights.resize(numyears);
      tweights.ones();
    }
    
    // Here we set up the vector of chosen occasions, sampled from all possible occasions
    tweights = tweights / sum(tweights);
    theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears, theclairvoyant, true, tweights);
    
    arma::vec startvec(matdim);
    startvec.ones();
    startvec = startvec / matdim; // This is the start vector for w and v calculations
    
    // Stochastic elasticities
    
    // The next section creates stable stage and rep value vectors arranged in
    // matrix format. The first two are general for whatever has been input,
    // whether historical or ahistorical, while the next two are specifically
    // for ahistorical versions of historical inputs
    arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::rowvec Rvecmat(theclairvoyant);
    wprojection.zeros();
    vprojection.zeros();
    Rvecmat.zeros();
    
    // Here we run the control loop to create the w and v values we need from the reference annual matrices
    arma::mat crazy_prophet = proj3(startvec, ref_byyear, theprophecy, 1, 0, 0);
    wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1),
      theclairvoyant);
    vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0,
      ((startvec.n_elem * 3) - 1), theclairvoyant);
    Rvecmat = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3),
      theclairvoyant); // Rvec
    
    arma::mat sensmat(matdim, matdim);
    arma::mat elasmean1(matdim, matdim);
    arma::mat elassd1(matdim, matdim);
    sensmat.zeros();
    elasmean1.zeros();
    elassd1.zeros();
    
    // All references should go to elasmat
    for (int j = burnin; j < theclairvoyant; j++) { // This is the main time loop for the sensitivity matrices, 
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wtplus1 = wprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // Numerator of key matrix equation
      arma::mat currentsens_den = (Rvecmat(j) * vtplus1.as_row() * wtplus1); // Denominator of key equation
      double cd_double = currentsens_den(0,0);
      sensmat = currentsens_num / (cd_double);
      
      elasmean1 = elasmean1 + ((sensmat % ref_matmean) /
        static_cast<double>(theclairvoyant - burnin));
      elassd1 = elassd1 + ((sensmat % ((as<arma::mat>(ref_byyear(theprophecy(j)))) - ref_matmean)) / 
        static_cast<double>(theclairvoyant - burnin));
    }
    elas_mean = List::create(elasmean1);
    elas_sd = List::create(elassd1);
    
    // Now the difference matrix estimation
    arma::mat diffmean1(matdim, matdim);
    arma::mat diffsd1(matdim, matdim);
    diffmean1.zeros();
    diffsd1.zeros();
    
    // This creates indices of 0 elements for use in the next control loop
    diffmean1 = log(ref_matmean);
    diffsd1 = log(ref_matsd);
    diffmean1.elem(find_nonfinite(diffmean1)).zeros();
    diffsd1.elem(find_nonfinite(diffsd1)).zeros();
    ref_matmean = diffmean1;
    ref_matsd = diffsd1;
    
    for (int i = 0; i < numpoppatches; i++) {
      diffmean1.zeros();
      diffsd1.zeros();
      
      diffmean1 = log(as<arma::mat>(poppatch_meanmat(i)));
      diffsd1 = log(as<arma::mat>(poppatch_sdmat(i)));
      diffmean1.elem(find_nonfinite(diffmean1)).zeros();
      diffsd1.elem(find_nonfinite(diffsd1)).zeros();
      
      poppatch_meanmat(i) = diffmean1;
      poppatch_sdmat(i) = diffsd1;
      
      diffmean1.zeros();
      diffsd1.zeros();
      diffmean1 = as<arma::mat>(poppatch_meanmat(i)) - ref_matmean;
      diffsd1 = as<arma::mat>(poppatch_sdmat(i)) - ref_matsd;
      
      // In Davison's original code, all elements equal to 0 in the reference
      // matrices must also be 0s in the difference matrices
      
      // Now the contributions
      diffmean1 = diffmean1 % as<arma::mat>(elas_mean(0));
      diffsd1 = diffsd1 % as<arma::mat>(elas_sd(0));
      
      if (i == 0) {
        cont_meanmat = List::create(diffmean1);
        cont_sdmat = List::create(diffsd1);
      } else {
        cont_meanmat.push_back(diffmean1);
        cont_sdmat.push_back(diffsd1);
      }
    }
  } else {
    // Sparse matrix analysis
    arma::sp_mat mat_mean(matdim, matdim);
    arma::sp_mat mat_sd(matdim, matdim);
    arma::sp_mat ref_matmean(matdim, matdim);
    arma::sp_mat ref_matsd(matdim, matdim);
    mat_mean.zeros();
    mat_sd.zeros();
    ref_matmean.zeros();
    ref_matsd.zeros();
  
    // First the pop/patch means and sds
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(poppatch_chosen(j)))) / 
          static_cast<double>(numpoppatch_chosen));
      }
      if (i == 0) {
        poppatch_meanmat = List::create(mat_mean);
      } else {
        poppatch_meanmat.push_back(mat_mean);
      }
      mat_mean.zeros();
    }
    
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      arma::vec mat_elems(numpoppatch_chosen);
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          mat_elems(k) = as<arma::mat>(Amats(poppatch_chosen(k)))(j);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(j) = stddev(mat_elems, 1);
        }
      }
      if (i == 0) {
        poppatch_sdmat = List::create(mat_sd);
      } else {
        poppatch_sdmat.push_back(mat_sd);
      }
      mat_sd.zeros();
    }
    
    // Here we create the reference matrix sets
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      mat_mean.zeros();
      mat_sd.zeros();
      ref_byyear = List::create(refmats);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) /
          static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum);
        mat_elems.zeros();
        
        for (int j = 0; j < refmatnum; j++) {
          mat_elems(j) = arma::sp_mat(as<arma::mat>(refmats(refnum(j))))(i);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(i) = stddev(mat_elems, 1);
        }
      }
      
      ref_matmean = mat_mean;
      ref_matsd = mat_sd;
      
    } else {
      if (refmatnum == Amatnum) {
        // Reference by year
        for (int i = 0; i < numyears; i++) {
          arma::uvec year_chosen = find(year2c == uniqueyears(i));
          int numyear_chosen = year_chosen.n_elem;
          
        for (int j = 0; j < numyear_chosen; j++) {
          mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(year_chosen(j)))) / 
            static_cast<double>(numyear_chosen));
        }
        if (i == 0) {
          ref_byyear = List::create(mat_mean);
        } else {
          ref_byyear.push_back(mat_mean);
        }
        mat_mean.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < numyears; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) / static_cast<double>(numyears));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < numyears; j++) {
            mat_elems(j) = as<arma::sp_mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        mat_mean.zeros();
        mat_sd.zeros();
        
      } else {
        // Reference by year
        for (int i = 0; i < refmatnum; i++) {
          if (i == 0) {
            ref_byyear = List::create(arma::sp_mat(as<arma::mat>(Amats(i))));
          } else {
            ref_byyear.push_back(arma::sp_mat(as<arma::mat>(Amats(i))));
          }
          mat_mean.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) /
            static_cast<double>(refmatnum));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < refmatnum; j++) {
            mat_elems(j) = as<arma::sp_mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        mat_mean.zeros();
        mat_sd.zeros();
      }
    }
    
    // Time weights
    if (tweights_.isNotNull()) {
      tweights = as<arma::vec>(tweights_);
      if (tweights.n_elem != numyears) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions in the reference matrix set.", false);
      }
      
    } else {
      tweights.resize(numyears);
      tweights.ones();
    }
    
    // Here we set up the vector of chosen occasions, sampled from all possible occasions
    tweights = tweights / sum(tweights);
    theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears, theclairvoyant, true, tweights);
    
    arma::vec startvec(matdim);
    startvec.ones();
    startvec = startvec / matdim; // This is the start vector for w and v calculations
    
    // Stochastic elasticities
    
    // The next section creates stable stage and rep value vectors arranged in
    // matrix format. The first two are general for whatever has been input,
    // whether historical or ahistorical, while the next two are specifically
    // for ahistorical versions of historical inputs
    arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::rowvec Rvecmat(theclairvoyant);
    wprojection.zeros();
    vprojection.zeros();
    Rvecmat.zeros();
    
    // Here we run the control loop to create the w and v values we need from the reference annual matrices
    arma::mat crazy_prophet = proj3sp(startvec, ref_byyear, theprophecy, 1, 0, 0);
    wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1),
      theclairvoyant);
    vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0,
      ((startvec.n_elem * 3) - 1), theclairvoyant);
    Rvecmat = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3),
      theclairvoyant); // Rvec
    
    arma::sp_mat sensmat(matdim, matdim);
    arma::sp_mat elasmean1(matdim, matdim);
    arma::sp_mat elassd1(matdim, matdim);
    sensmat.zeros();
    elasmean1.zeros();
    elassd1.zeros();
    
    // All references should go to elasmat
    for (int j = burnin; j < theclairvoyant; j++) { // Main loop for sensitivity matrices
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wtplus1 = wprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // Numerator of key matrix equation
      arma::mat currentsens_den = (Rvecmat(j) * vtplus1.as_row() * wtplus1); // Denominator of key equation
      double cd_double = currentsens_den(0,0);
      sensmat = arma::sp_mat(currentsens_num / (cd_double));
      
      elasmean1 = elasmean1 + ((sensmat % ref_matmean) /
        (static_cast<double>(theclairvoyant - burnin)));
      elassd1 = elassd1 + ((sensmat % ((as<arma::sp_mat>(ref_byyear(theprophecy(j)))) - ref_matmean)) / 
        (static_cast<double>(theclairvoyant - burnin)));
    }
    elas_mean = List::create(elasmean1);
    elas_sd = List::create(elassd1);
    
    // Now the difference and contribution matrix estimation
    ref_matmean = spmat_log(ref_matmean);
    ref_matsd = spmat_log(ref_matsd);
    arma::sp_mat diffmean1(matdim, matdim);
    arma::sp_mat diffsd1(matdim, matdim);
    
    for (int i = 0; i < numpoppatches; i++) {
      diffmean1.zeros();
      diffsd1.zeros();
      
      poppatch_meanmat(i) = spmat_log(as<arma::sp_mat>(poppatch_meanmat(i)));
      poppatch_sdmat(i) = spmat_log(as<arma::sp_mat>(poppatch_sdmat(i)));
      
      diffmean1 = as<arma::sp_mat>(poppatch_meanmat(i)) - ref_matmean;
      diffsd1 = as<arma::sp_mat>(poppatch_sdmat(i)) - ref_matsd;
      diffmean1 = diffmean1 % as<arma::sp_mat>(elas_mean(0));
      diffsd1 = diffsd1 % as<arma::sp_mat>(elas_sd(0));
      
      if (i == 0) {
        cont_meanmat = List::create(arma::mat(diffmean1));
        cont_sdmat = List::create(arma::mat(diffsd1));
      } else {
        cont_meanmat.push_back(arma::mat(diffmean1));
        cont_sdmat.push_back(arma::mat(diffsd1));
      }
    }
  }
  return List::create(Named("cont_mean") = cont_meanmat,
    _["cont_sd"] = cont_sdmat);
}