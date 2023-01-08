#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;

//' Create Stageframe for Population Matrix Projection Analysis
//' 
//' Function \code{sf_create()} returns a data frame describing each ahistorical
//' life history stage in the life history model. This data frame can be used as 
//' input into MPM creation functions including \code{\link{flefko3}()}, 
//' \code{\link{flefko2}()}, \code{\link{aflefko2}()}, \code{\link{rlefko3}()},
//' and \code{\link{rlefko2}()}, in which it determines how each stage is
//' treated during matrix estimation.
//' 
//' @name sf_create
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
//' Vectors used to create a stageframe may not mix \code{NA} values with
//' non-\code{NA} values.
//' 
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
//' If importing an IPM rather than building one with \code{lefko3}: Using the
//' \code{vrm_input} approach to building function-based MPMs with provided
//' linear model slope coefficients requires careful attention to the
//' stageframe. Although no hfv data frame needs to be entered in this instance,
//' stages for which vital rates are to be estimated via linear models
//' parameterized with coefficients provided via function
//' \code{\link{vrm_import}()} should be marked as occurring within the dataset,
//' while stages for which the provided coefficients should not be used should
//' be marked as not occurring within the dataset.
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
  }
  
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
    Rf_warningcall(R_NilValue, "High ipmbin numbers may lead to dramatic decreases in statistical power and overparameterized matrices.\n");
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

//' Calculate Actual Stage, Age, Stage-Pair, or Age-Stage Distributions
//' 
//' Function \code{actualstage3()} shows the frequencies and proportions of
//' each stage, stage pair, age-stage, or age in each year.
//' 
//' @name actualstage3
//' 
//' @param data A demographic dataset in hfv format.
//' @param check_stage A logical value indicating whether to assess frequencies
//' and proportions of stages. Defaults to \code{TRUE}.
//' @param check_age A logical value indicating whether to assess frequencies and
//' proportions of ages. Defaults to \code{FALSE}.
//' @param historical A logical value indicating whether the stage structure
//' should be ahistorical (\code{FALSE}) or historical (\code{TRUE}). Defaults to
//' \code{FALSE}.
//' @param year2 A string value indicating the name of the variable coding for
//' monitoring occasion at time \emph{t}. Defaults to \code{"year2"}.
//' @param indices A vector of three strings, indicating the stage indices for
//' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
//' Defaults to \code{c("stage3index", "stage2index", "stage1index")}.
//' @param stagecol A vector of three strings, indicating the stage name columns
//' for times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
//' Defaults to \code{stagecol = c("stage3", "stage2", "stage1")}.
//' @param agecol A single string indicating the age of individuals in time
//' \emph{t}. Defaults to \code{"obsage"}.
//' @param remove_stage A string vector indicating the names of stages to remove
//' from consideration. Defaults to \code{"NotAlive"}.
//' @param t1_allow A string vector indicating which stages to be removed should
//' be allowed in the stage at time \emph{t}-1 portion of historical stage
//' pairs, if \code{historical = TRUE}. Defaults to \code{"NotAlive"}. Can also
//' be set to \code{"none"}.
//' 
//' @return A data frame with the following variables:
//' \item{rowid}{A string identifier term, equal to the monitoring occasion in
//' time \emph{t} and the stage index.}
//' \item{stageindex}{The stageframe index of the stage. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage}{The name of each stage, or \code{NA}. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage2}{The name of the stage in time \emph{t}. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage1}{The name of the stage in time \emph{t}-1, or \code{NA}. Only
//' output if \code{check_stage = TRUE}.}
//' \item{age}{The age at time \emph{t}. Only output if \code{check_age = TRUE}.}
//' \item{year2}{Monitoring occasion in time \emph{t}.}
//' \item{frequency}{The number of individuals in the respective stage and time.}
//' \item{actual_prop}{The proportion of individuals alive in time \emph{t} in
//' the respective stage.}
//' 
//' @section Notes:
//' This function produces frequencies and proportions of stages in hfv formatted
//' data using stage index variables rather than stage name variables, and so
//' requires the former. The latter is only required if the user wants to know
//' the associated stage names.
//' 
//' Frequencies and proportions will be calculated for all times, including the
//' last time, which is generally found in the \code{stage3} columns of the last
//' \code{year2} entry in object \code{data}. The default is to treat the
//' \code{year2} entry for that time as \code{max(year2) + 1}.
//' 
//' If \code{check_stage = TRUE} and \code{check_age = FALSE}, then this function
//' will assess frequencies and proportions of stages or historical stage-pairs.
//' If both \code{check_stage = TRUE} and \code{check_age = TRUE}, then this
//' function will assess frequencies and proportions of age-stages. If
//' \code{check_stage = FALSE} and \code{check_age = TRUE}, then the frequencies
//' and proportions of ages only will be assessed.
//' 
//' Note that no stageframe is required for this function to operate. Stage
//' names and their order are inferred directly from the object \code{data}.
//' 
//' @examples
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 3, 6, 11, 19.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1.5, 1.5, 3.5, 5)
//' comments <- c("Dormant seed", "1st yr protocorm", "2nd yr protocorm",
//'   "3rd yr protocorm", "Seedling", "Dormant adult",
//'   "Extra small adult (1 shoot)", "Small adult (2-4 shoots)",
//'   "Medium adult (5-7 shoots)", "Large adult (8-14 shoots)",
//'   "Extra large adult (>14 shoots)")
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset, 
//'   binhalfwidth = binvec, comments = comments)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE, age_offset = 4)
//' 
//' all_stage_props_ah <- actualstage3(cypraw_v1)
//' all_stage_props_h <- actualstage3(cypraw_v1, historical = TRUE)
//' all_stage_props_h_NANotAllow <- actualstage3(cypraw_v1, historical = TRUE,
//'   t1_allow = "none")
//' all_stage_props_as <- actualstage3(cypraw_v1, check_age = TRUE)
//' all_age_props <- actualstage3(cypraw_v1, check_stage = FALSE,
//'   check_age = TRUE)
//' 
//' @export actualstage3
// [[Rcpp::export(actualstage3)]]
List actualstage3(RObject data, bool check_stage = true, bool check_age = false,
  bool historical = false, Nullable<RObject> year2 = R_NilValue,
  Nullable<RObject> indices = R_NilValue, Nullable<RObject> stagecol = R_NilValue,
  Nullable<RObject> agecol= R_NilValue, Nullable<RObject> remove_stage = R_NilValue,
  Nullable<RObject> t1_allow = R_NilValue) {
  
  DataFrame data_;
  
  String year2_;
  String agecol_;
  StringVector remove_stage_;
  StringVector t1_allow_;
  
  int year2_num_;
  int agecol_num_;
  IntegerVector remove_stage_num_;
  IntegerVector remove_stage_num_t1a;
  IntegerVector t1_allow_num_;
  
  bool year2_num_yn = false;
  bool agecol_num_yn = false;
  bool remove_stage_num_yn = false;
  bool t1_allow_num_yn = false;
  
  StringVector indices_;
  StringVector stagecol_;
  
  IntegerVector indices_num_;
  IntegerVector stagecol_num_;
  
  bool indices_num_yn = false;
  bool stagecol_num_yn = false;
  
  int stage3index_;
  int stage2index_;
  int stage1index_;
  int stage3_;
  int stage2_;
  int stage1_;
  
  bool stages3_supplied = false;
  bool stages2_supplied = false;
  bool stages1_supplied = false;
  bool indices3_supplied = false;
  bool indices2_supplied = false;
  bool indices1_supplied = false;
  bool ages_supplied = false;
  bool years_supplied = false;
  
  if (!check_age && !check_stage) {
    throw Rcpp::exception("Options check_age and check_stage cannot both be FALSE.", false);
  }
  
  if (year2.isNotNull()) {
    if (is<StringVector>(year2)) {
      StringVector year2_long = as<StringVector>(year2);
      year2_ = year2_long(0);
      
    } else if (is<IntegerVector>(year2)) {
      IntegerVector year2_long = as<IntegerVector>(year2);
      year2_num_ = year2_long(0);
      year2_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object year2 must be either a variable name or a column number from object data.", 
        false);
    }
    
  } else year2_ = "year2";
  
  if (agecol.isNotNull()) {
    if (is<StringVector>(agecol)) {
      StringVector agecol_long = as<StringVector>(agecol);
      agecol_ = agecol_long(0);
      
    } else if (is<IntegerVector>(agecol)) {
      IntegerVector agecol_long = as<IntegerVector>(agecol);
      agecol_num_ = agecol_long(0);
      agecol_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object agecol must be either a variable name or a column number from object data.", 
        false);
    }
    
  } else agecol_ = "obsage";
  
  if (remove_stage.isNotNull()) {
    if (is<StringVector>(remove_stage)) {
      remove_stage_ = as<StringVector>(remove_stage);
      
      if (remove_stage_.length() == 1) {
        if (remove_stage_(0) == "none") {
          remove_stage_ = {};
          remove_stage_num_ = {};
          
        }
      }
      
    } else if (is<IntegerVector>(remove_stage)) {
      remove_stage_num_ = as<IntegerVector>(remove_stage);
      remove_stage_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object remove_stage must be either a string of variable names or a vector of column number from object data.", 
        false);
    }
    
  } else remove_stage_ = {"NotAlive"};
  
  if (t1_allow.isNotNull()) {
    if (is<StringVector>(t1_allow)) {
      t1_allow_ = as<StringVector>(t1_allow);
      
      if (t1_allow_.length() == 1) {
        if (t1_allow_(0) == "none") {
          t1_allow_ = {};
          t1_allow_num_ = {};
          
        }
      }
      
    } else if (is<IntegerVector>(t1_allow)) {
      t1_allow_num_ = as<IntegerVector>(t1_allow);
      t1_allow_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object t1_allow must be either a string of variable names or a vector of column number from object data.", false);
    }
  } else t1_allow_ = {"NotAlive"};
  
  if (indices.isNotNull()) {
    if (is<StringVector>(indices)) {
      indices_ = as<StringVector>(indices);
      
    } else if (is<IntegerVector>(indices)) {
      indices_num_ = as<IntegerVector>(indices);
      indices_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object indices must be a vector of variable names or column numbers from object data.", 
        false);
    }
    
  } else indices_ = {"stage3index", "stage2index", "stage1index"};
  
  if (stagecol.isNotNull()) {
    if (is<StringVector>(stagecol)) {
      stagecol_ = as<StringVector>(stagecol);
      
    } else if (is<IntegerVector>(stagecol)) {
      stagecol_num_ = as<IntegerVector>(stagecol);
      stagecol_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object stagecol must be a vector of variable names or column numbers from object data.", 
        false);
    }
    
  } else stagecol_ = {"stage3", "stage2", "stage1"};
  
  if (is<DataFrame>(data)) {
    data_ = as<DataFrame>(data);
    
  } else {
    throw Rcpp::exception("Object data must be a hfv data frame.", false);
  }
  
  // Check data frame for the correct variables, and pull out variables as vectors
  CharacterVector data_vars = data_.attr("names");
  int data_vars_num = data_vars.length();
  int data_rows = data_.nrows();
  int indices_length = indices_.length();
  int stagecol_length = stagecol_.length();
  
  if (historical) {
    if (check_age) {
      throw Rcpp::exception("Package lefko3 does not currently support historical age-by-stage analyses.", false);
    }
    
    if (indices_length < 3) {
      throw Rcpp::exception("Object indices must contain stage index variables for times t+1, t, and t-1 if historical = TRUE.", false);
    }
    
    if (stagecol_length < 3) {
      throw Rcpp::exception("Object stagecol must contain stage classifications for times t+1, t, and t-1 if historical = TRUE.", false);
    }
  } else {
    if (indices_length < 2) {
      throw Rcpp::exception("Object indices must contain stage index variables for times t+1 and t if historical = FALSE.", false);
    }
    
    if (stagecol_length < 2) {
      throw Rcpp::exception("Object stagecol must contain stage classifications for times t+1 and t if historical = FALSE.", false);
    }
  }
  
  // Loop to find column numbers of variables given as strings
  for (int i = 0; i < data_vars_num; i++) {
    if (!year2_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), year2_)) {
        year2_num_ = i;
        years_supplied = true;
      }
    }
    
    if (!agecol_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), agecol_)) {
        agecol_num_ = i;
        ages_supplied = true;
      }
    }
    
    if (!indices_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(0)))) {
        stage3index_ = i;
        indices3_supplied = true;
      }
      
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(1)))) {
        stage2index_ = i;
        indices2_supplied = true;
      }
      
      if (historical) {
        if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(2)))) {
          stage1index_ = i;
          indices1_supplied = true;
        }
      }
    }
    
    if (!stagecol_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(0)))) {
        stage3_ = i;
        stages3_supplied = true;
      }
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(1)))) {
        stage2_ = i;
        stages2_supplied = true;
      }
      
      if (historical) {
        if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(2)))) {
          stage1_ = i;
          stages1_supplied = true;
        }
      }
    }
  }
  
  // Check if variable column numbers actually fit within data frame provided
  if (year2_num_yn) {
    if (year2_num_ < 0 || year2_num_ >= data_vars_num) {
      throw Rcpp::exception("Variable year2 not defined within data frame as input.", false);
      
    } else {
      years_supplied = true;
    }
  }
  
  if (agecol_num_yn) {
    if (agecol_num_ < 0 || agecol_num_ >= data_vars_num) {
      throw Rcpp::exception("Variable agecol not defined within data frame as input.", false);
      
    } else {
      ages_supplied = true;
    }
  }
  
  if (indices_num_yn) {
    if (indices_num_(0) >= 0 && indices_num_(0) < data_vars_num) {
      stage3index_ = indices_num_(0);
      indices3_supplied = true;
    }
    
    if (indices_num_(1) >= 0 && indices_num_(1) < data_vars_num) {
      stage2index_ = indices_num_(1);
      indices2_supplied = true;
    }
    if (historical) {
      if (indices_num_(2) >= 0 && indices_num_(2) < data_vars_num) {
        stage1index_ = indices_num_(2);
        indices1_supplied = true;
      }
    }
  }
  
  if (stagecol_num_yn) {
    if (stagecol_num_(0) >= 0 && stagecol_num_(0) < data_vars_num) {
      stage3_ = stagecol_num_(0);
      stages3_supplied = true;
    }
    
    if (stagecol_num_(1) >= 0 && stagecol_num_(1) < data_vars_num) {
      stage2_ = stagecol_num_(1);
      stages2_supplied = true;
    }
    
    if (historical) {
      if (stagecol_num_(2) >= 0 && stagecol_num_(2) < data_vars_num) {
        stage1_ = stagecol_num_(2);
        stages1_supplied = true;
      }
    }
  }
  
  // New data frame variable assignments
  arma::ivec data_year2;
  if (years_supplied) data_year2 = as<arma::ivec>(data_[year2_num_]);
  
  if (data_year2.n_elem != data_rows) {
    throw Rcpp::exception("Object year2 does not contain a valid variable.", false);
  }
  
  IntegerVector data_agecol;
  
  if (check_age) {
    if (ages_supplied) {
      data_agecol = data_[agecol_];
      
      if (data_agecol.length() != data_rows) {
        throw Rcpp::exception("Object agecol does not contain a valid variable for age at time t.", false);
      }
    } else {
      throw Rcpp::exception("Object agecol must be provided if check_age = TRUE.", false);
    }
  }
  
  CharacterVector data_stage3;
  CharacterVector data_stage2;
  CharacterVector data_stage1;
  
  arma::ivec data_stage3index;
  arma::ivec data_stage2index;
  arma::ivec data_stage1index;
  
  if (check_stage) {
    if (stages3_supplied) {
      data_stage3 = as<CharacterVector>(data_[stage3_]);
      
      // Rcout << "\n Length of data_stage3: " << data_stage3.length() << "\n";
      // Rcout << "First entry in data_stage3: " << data_stage3(0) << "\n";
      
      if (data_stage3.length() != data_rows) {
        throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t+1.", false);
      }
    }
    if (stages2_supplied) {
      data_stage2 = as<CharacterVector>(data_[stage2_]);
      
      // Rcout << "\n Length of data_stage2: " << data_stage2.length() << "\n";
      // Rcout << "First entry in data_stage2: " << data_stage2(0) << "\n";
      
      if (data_stage2.length() != data_rows) {
        throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t.", false);
      }
    }
    
    if (indices3_supplied) {
      data_stage3index = as<arma::ivec>(data_[stage3index_]);
      
      // Rcout << "\n Length of data_stage3index: " << data_stage3index.n_elem << "\n";
      // Rcout << "First entry in data_stage3index: " << data_stage3index(0) << "\n";
      
      if (data_stage3index.n_elem != data_rows) {
        throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t+1.", false);
      }
    }
    if (indices2_supplied) {
      data_stage2index = as<arma::ivec>(data_[stage2index_]);
      
      // Rcout << "First entry in data_stage2index: " << data_stage2index(0) << "\n";
      // Rcout << "\n Length of data_stage2index: " << data_stage2index.n_elem << "\n";
      
      if (data_stage2index.n_elem != data_rows) {
        throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t.", false);
      }
    }
    
    if (!stages3_supplied && !indices3_supplied) {
      throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.", false);
    }
    if (!stages2_supplied && !indices2_supplied) {
      throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.", false);
    }
    
    if (historical) {
      if (stages1_supplied) {
        data_stage1 = data_[stage1_];
        
        if (data_stage1.length() != data_rows) {
          throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t-1.", false);
        }
      }
      if (indices3_supplied) {
        data_stage1index = as<arma::ivec>(data_[stage1index_]);
        
        if (data_stage1index.n_elem != data_rows) {
          throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t-1.", false);
        }
      }
      
      if (!stages1_supplied && !indices1_supplied) {
        throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.", false);
      }
    }
  }
  
  // Now we develop vectors of all values
  arma::ivec all_year2 = unique(data_year2);
  int year2_num = all_year2.n_elem;
  int years_num = year2_num + 1;
  
  arma::ivec all_years (years_num);
  for (int i = 0; i < year2_num; i++) {
    all_years(i) = all_year2(i);
  }
  all_years(year2_num) = all_years(year2_num - 1) + 1;
  
  CharacterVector all_stage3;
  CharacterVector all_stage2;
  CharacterVector all_stage1;
  CharacterVector all_stages;
  CharacterVector all_stages_t1a;
  
  arma::ivec all_stage3index;
  arma::ivec all_stage2index;
  arma::ivec all_stage1index;
  arma::ivec all_stageindices;
  arma::ivec all_stageindices_t1a;
  
  int num_stages = 0;
  int num_stages_t1a = 0;
  
  if (check_stage) {
    if (indices3_supplied && indices2_supplied) {
      all_stage3index = unique(data_stage3index);
      all_stage2index = unique(data_stage2index);
      arma::ivec all_stage32indices = unique(join_cols(all_stage3index, all_stage2index));
      
      if (historical) {
        all_stage1index = unique(data_stage1index);
        all_stageindices = unique(join_cols(all_stage32indices, all_stage1index));
      } else {
        all_stageindices = all_stage32indices;
      }
      
      int s_length = all_stageindices.n_elem;
      // Rcout << "\n Length of all_stageindices: " << s_length << "\n";
      
      if (stages3_supplied && stages2_supplied) {
        CharacterVector temp_stages (s_length);
        for (int i = 0; i < s_length; i++) {
          arma::uvec found_stages3 = find(data_stage3index == all_stageindices(i));
          
          if (found_stages3.n_elem > 0) {
            temp_stages(i) = data_stage3(found_stages3(0));
          } else {
            arma::uvec found_stages2 = find(data_stage2index == all_stageindices(i));
            
            if (found_stages2.n_elem > 0) {
              temp_stages(i) = data_stage2(found_stages2(0));
              
            } else if (historical && stages1_supplied) {
              arma::uvec found_stages1 = find(data_stage1index == all_stageindices(i));
              
              if (found_stages1.n_elem > 0) {
                temp_stages(i) = data_stage1(found_stages1(0));
              }
            }
          }
        }
        all_stages = temp_stages;
        num_stages = all_stages.length();
        num_stages_t1a = all_stages.length();
      } else {
        num_stages = all_stageindices.n_elem;
        num_stages_t1a = all_stageindices.n_elem;
      }
    } else if (stages3_supplied && stages2_supplied) {
      all_stage3 = sort_unique(data_stage3);
      all_stage2 = sort_unique(data_stage2);
      
      CharacterVector all_stage32 = union_(all_stage3, all_stage2);
      
      if (historical) {
        all_stage1 = sort_unique(data_stage1);
        
        all_stages = union_(all_stage32, all_stage1);
        
      } else {
        all_stages = all_stage32;
      }
      
      num_stages = all_stages.length();
      num_stages_t1a = all_stages.length();
      IntegerVector temp_stageindices = seq(1, num_stages);
      all_stageindices = as<arma::ivec>(temp_stageindices);
    }
    
    if (remove_stage_.length() != 0 && !remove_stage_num_yn) {
      IntegerVector transfer_remove_stage_num (remove_stage_.length());
      IntegerVector transfer_remove_stage_num_t1a (remove_stage_.length());
      int internal_counter = 0;
      int internal_counter_t1a = 0;
      
      for (int i = 0; i < num_stages; i++) {
        for (int j = 0; j < remove_stage_.length(); j++) {
          if (all_stages(i) == remove_stage_(j)) {
            transfer_remove_stage_num(internal_counter) = i;
            remove_stage_num_yn = true;
            internal_counter++;
            
            if (historical) {
              if (t1_allow_.length() > 0) {
                
                bool found_stage = false;
                for (int k = 0; k < t1_allow_.length(); k++) {
                  if (all_stages(i) == t1_allow_(k)) {
                    found_stage = true;
                    // Rcout << "Found stage (1) " << all_stages(i) << "\n";
                    
                    transfer_remove_stage_num_t1a = shrink(transfer_remove_stage_num_t1a);
                  }
                }
                
                if (!found_stage) {
                  transfer_remove_stage_num_t1a(internal_counter_t1a) = i;
                  internal_counter_t1a++;
                }
              } else if (t1_allow_num_yn) {
                bool found_stage = false;
                for (int k = 0; k < t1_allow_num_.length(); k++) {
                  if (all_stageindices(i) == t1_allow_num_(k)) {
                    found_stage = true;
                    // Rcout << "Found stage (2) " << all_stages(i) << "\n";
                    
                    transfer_remove_stage_num_t1a = shrink(transfer_remove_stage_num_t1a);
                  }
                }
                
                if (!found_stage) {
                  transfer_remove_stage_num_t1a(internal_counter_t1a) = i;
                  internal_counter_t1a++;
                }
              }
            }
          }
        }
      }
      remove_stage_num_ = transfer_remove_stage_num;
      if (historical) {
        remove_stage_num_t1a = clone(transfer_remove_stage_num_t1a);
      } else {
        remove_stage_num_t1a = clone(transfer_remove_stage_num);
      }
      
      
      
      // Error checks
      // Rcout << "all_stages\n";
      // for (int i=0; i < all_stages.length(); i++) {
      //   Rcout << "stage " << i << " " << all_stages(i) << "\n";
      // }
      // 
      // Rcout << "remove_stage_num_\n";
      // for (int i=0; i < remove_stage_num_.length(); i++) {
      //   Rcout << "remove_stage_num_ " << i << " " << all_stages(remove_stage_num_(i)) << "\n";
      // }
      // 
      // Rcout << "remove_stage_num_\n";
      // for (int i=0; i < remove_stage_num_t1a.length(); i++) {
      //   Rcout << "remove_stage_num_t1a " << i << " " << all_stages(remove_stage_num_t1a(i)) << "\n";
      // }
      
      if (!remove_stage_num_yn) {
        Rf_warningcall(R_NilValue,
          "Stage(s) provided in option remove_stage could not be found and so will be ignored.");
      }
      
      IntegerVector unique_rsn = unique(remove_stage_num_);
      IntegerVector unique_rsn_t1a = unique(remove_stage_num_t1a);
      
      if (unique_rsn.length() < remove_stage_num_.length() || 
          unique_rsn_t1a.length() < remove_stage_num_t1a.length()) {
        Rf_warningcall(R_NilValue, "Some stages supplied in either option remove_stage or option t1_allow are duplicates and will be ignored.");
        
        remove_stage_num_ = clone(unique_rsn);
        unique_rsn_t1a = clone(unique_rsn_t1a);
      }
      
      // Redetermination of the number of stages to be used
      IntegerVector sorted_remst = remove_stage_num_.sort();
      IntegerVector rev_sorted_remst = rev(sorted_remst);
      
      CharacterVector rep_all_stages = clone(all_stages);
      CharacterVector orig_all_stages = clone(all_stages);
      IntegerVector rep_all_stageindices = as<IntegerVector>(wrap(all_stageindices));
      IntegerVector orig_rep_all_stageindices = clone(rep_all_stageindices);
      
      for (int i = 0; i < rev_sorted_remst.length(); i++) {
        rep_all_stages.erase(rev_sorted_remst(i));
        rep_all_stageindices.erase(rev_sorted_remst(i));
      }
      all_stages = rep_all_stages;
      all_stageindices = as<arma::ivec>(rep_all_stageindices);
      num_stages = all_stages.length();
      
      if (historical) {
        IntegerVector sorted_remst_t1a = remove_stage_num_t1a.sort();
        IntegerVector rev_sorted_remst_t1a = rev(sorted_remst_t1a);
        
        CharacterVector rep_all_stages_t1a = clone(orig_all_stages);
        IntegerVector rep_all_stageindices_t1a = clone(orig_rep_all_stageindices);
        
        for (int i = 0; i < rev_sorted_remst_t1a.length(); i++) {
          rep_all_stages_t1a.erase(rev_sorted_remst_t1a(i));
          rep_all_stageindices_t1a.erase(rev_sorted_remst_t1a(i));
        }
        all_stages_t1a = rep_all_stages_t1a;
        all_stageindices_t1a = as<arma::ivec>(rep_all_stageindices_t1a);
        num_stages_t1a = all_stages_t1a.length();
      }
    } else if (historical) {
      all_stages_t1a = all_stages;
      all_stageindices_t1a = all_stageindices;
      num_stages_t1a = all_stages.length();
    }
  }
  
  IntegerVector all_ages;
  int num_ages = 0;
  
  if (check_age) {
    IntegerVector all_age2 = sort_unique(data_agecol);
    int num_age2 = all_age2.length();
    num_ages = num_age2 + 1;
    
    IntegerVector all_ages_ (num_age2 + 1);
    
    for (int i = 0; i < num_age2; i++) {
      all_ages_(i) = all_age2(i);
    }
    all_ages_(num_age2) = all_ages_(num_age2 - 1) + 1;
    
    all_ages = all_ages_;
  }
  
  // New data frame variables
  int new_df_rows = 0;
  
  if (check_stage && !check_age) {
    if (!historical) {
      new_df_rows = years_num * num_stages;
      
    } else {
      new_df_rows = years_num * num_stages * num_stages_t1a;
    }
    
  } else if (!check_stage && check_age) {
    new_df_rows = years_num * num_ages;
    
  } else if (check_stage && check_age) {
    new_df_rows = years_num * num_stages * num_ages;
  }
  
  IntegerVector new_row_names (new_df_rows);
  CharacterVector new_row_id (new_df_rows);
  IntegerVector new_stageindex (new_df_rows);
  CharacterVector new_stage (new_df_rows);
  CharacterVector new_stage2 (new_df_rows);
  CharacterVector new_stage1 (new_df_rows);
  IntegerVector new_age (new_df_rows);
  IntegerVector new_year (new_df_rows);
  IntegerVector new_frequency (new_df_rows);
  NumericVector new_actualprop (new_df_rows);
  
  // Rcout << "\n new_df_rows: " << new_df_rows << "\n";
  
  if (check_stage && !check_age) {
    if (!historical) {
      IntegerVector first_years = rep(all_years(0), num_stages);
      CharacterVector first_stages2 = clone(all_stages);
      IntegerVector first_indices = as<IntegerVector>(wrap(all_stageindices));
      
      CharacterVector dud_stages (num_stages);
      for (int i = 0; i < num_stages; i++) {
        dud_stages(i) = "";
      }
      
      CharacterVector first_stages1 (num_stages);
      for (int i = 0; i < num_stages; i++) {
        first_stages1(i) = "";
      }
      
      for (int i = 1; i < years_num; i++) {
        IntegerVector next_years = rep(all_years(i), num_stages);
        first_years = concat_int(first_years, next_years);
        
        first_stages2 = concat_str(first_stages2, all_stages);
        first_stages1 = concat_str(first_stages1, dud_stages);
        
        first_indices = concat_int(first_indices, as<IntegerVector>(wrap(all_stageindices)));
      }
      new_year = first_years;
      
      new_stage2 = first_stages2;
      new_stage1 = first_stages1;
      new_stageindex = first_indices;
      
    } else {
      IntegerVector first_years = rep(all_years(0), num_stages * num_stages_t1a);
      CharacterVector first_stages2 = clone(all_stages);
      
      CharacterVector first_stages1_short (num_stages_t1a);
      CharacterVector first_stages1;
      
      for (int i = 0; i < num_stages_t1a; i++) {
        for (int j = 0; j < num_stages_t1a; j++) {
          first_stages1_short(j) = all_stages_t1a(i);
        }
        
        if (i > 0) {
          first_stages1 = concat_str(first_stages1, first_stages1_short);
          first_stages2 = concat_str(first_stages2, all_stages);
        } else {
          first_stages1 = clone(first_stages1_short);
        }
      }
      
      CharacterVector new_first_stages1 = clone(first_stages1);
      CharacterVector new_first_stages2 = clone(first_stages2);
      
      for (int i = 1; i < years_num; i++) {
        IntegerVector next_years = rep(all_years(i), num_stages * num_stages_t1a);
        first_years = concat_int(first_years, next_years);
        
        new_first_stages1 = concat_str(new_first_stages1, first_stages1);
        new_first_stages2 = concat_str(new_first_stages2, first_stages2);
      }
      
      new_year = first_years;
      
      new_stage2 = new_first_stages2;
      new_stage1 = new_first_stages1;
    }
    
  } else if (!check_stage && check_age) {
  
    IntegerVector first_years = rep(all_years(0), num_ages);
    IntegerVector first_ages = clone(all_ages);
    
    for (int i = 1; i < years_num; i++) {
      IntegerVector next_years = rep(all_years(i), num_ages);
      first_years = concat_int(first_years, next_years);
      
      first_ages = concat_int(first_ages, all_ages);
    }
    
    new_year = first_years;
    new_age = first_ages;
    
  } else if (check_stage && check_age) {
    
    CharacterVector first_stages2 = clone(all_stages);
    IntegerVector first_indices = as<IntegerVector>(wrap(all_stageindices));
    
    CharacterVector dud_stages (num_stages);
    for (int i = 0; i < num_stages; i++) {
      dud_stages(i) = "";
    }
    
    CharacterVector first_stages1 (num_stages);
    for (int i = 0; i < num_stages; i++) {
      first_stages1(i) = "";
    }
    
    IntegerVector first_ages = rep(all_ages(0), num_stages);
    for (int i = 1; i < num_ages; i++) {
      IntegerVector next_ages = rep(all_ages(i), num_stages);
      first_ages = concat_int(first_ages, next_ages);
      
      first_stages2 = concat_str(first_stages2, all_stages);
      first_stages1 = concat_str(first_stages1, dud_stages);
    }
    
    IntegerVector first_years = rep(all_years(0), num_stages * num_ages);
    IntegerVector cfirst_ages = clone(first_ages);
    CharacterVector cfirst_stages2 = clone(first_stages2);
    CharacterVector cfirst_stages1 = clone(first_stages1);
    
    for (int i = 1; i < years_num; i++) {
      IntegerVector next_years = rep(all_years(i), num_stages * num_ages);
      first_years = concat_int(first_years, next_years);
      
      first_ages = concat_int(first_ages, cfirst_ages);
      
      first_stages2 = concat_str(first_stages2, cfirst_stages2);
      first_stages1 = concat_str(first_stages1, cfirst_stages1);
      
      first_indices = concat_int(first_indices, as<IntegerVector>(wrap(all_stageindices)));
    }
    
    new_year = first_years;
    new_age = first_ages;
    new_stage2 = first_stages2;
    new_stage1 = first_stages1;
    new_stageindex = first_indices;
    
  }
  
  /*
  Rcout << "Pre-loop\n";
  Rcout << "\n Length of new_row_names: " << new_row_names.length() << "\n";
  Rcout << " First entry: " << new_row_names(0) << "\n";
  
  Rcout << "\n Length of new_row_id: " << new_row_id.length() << "\n";
  Rcout << " First entry: " << new_row_id(0) << "\n";
  
  Rcout << "\n Length of new_age: " << new_age.length() << "\n";
  Rcout << " First entry: " << new_age(0) << "\n";
  
  Rcout << "\n Length of new_stage: " << new_stage.length() << "\n";
  Rcout << " First entry: " << new_stage(0) << "\n";
  
  Rcout << "\n Length of new_stage2: " << new_stage2.length() << "\n";
  Rcout << " First entry: " << new_stage2(0) << "\n";
  
  Rcout << "\n Length of new_stage1: " << new_stage1.length() << "\n";
  Rcout << " First entry: " << new_stage1(0) << "\n";
  
  Rcout << "V";
  */
  
  // Main analysis loop
  IntegerVector pop_by_years (years_num);
  
  for (int i = 0; i < new_df_rows; i++) {
    new_row_names(i) = i + 1;
    
    new_row_id(i) = new_year(i);
    new_row_id(i) += " ";
    
    // Here we set the stage and row id designations
    if (check_stage && !check_age) {
      if (!historical) {
        new_stage(i) = new_stage2(i);
        
      } else {
        new_stage(i) = new_stage2(i);
        new_stage(i) += " ";
        new_stage(i) += new_stage1(i);
        
      }
      new_row_id(i) += new_stage(i);
      
    } else if (!check_stage && check_age) {
      new_row_id(i) += new_age(i);
      
    } else if (check_stage && check_age) {
      new_stage(i) = new_stage2(i);
      
      new_row_id(i) += new_stage(i);
      new_row_id(i) += " ";
      new_row_id(i) += new_age(i);
      
    }
    
    // Now we'll find the frequencies of age-stage-year combos
    arma::uvec year_guys = find(data_year2 == new_year(i));
    int year_guys_length = year_guys.n_elem;
    
    if (check_stage && !check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (stringcompare_hard(as<std::string>(data_stage2(year_guys(j))), as<std::string>(new_stage2(i)))) {
            if (!historical) {
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            } else {
            
              if (stringcompare_hard(as<std::string>(data_stage1(year_guys(j))), as<std::string>(new_stage1(i)))) {
                arma::uvec year_found = find(all_years == new_year(i));
                pop_by_years(year_found(0))++;
                
                new_frequency(i)++;
              }
            }
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2 == new_year(i) - 1);
        int year_guys_length3 = year_guys3.n_elem;
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (stringcompare_hard(as<std::string>(data_stage3(year_guys3(j))), as<std::string>(new_stage2(i)))) {
            if (!historical) {
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            } else {
              if (stringcompare_hard(as<std::string>(data_stage2(year_guys3(j))), as<std::string>(new_stage1(i)))) {
                arma::uvec year_found = find(all_years == new_year(i));
                pop_by_years(year_found(0))++;
                
                new_frequency(i)++;
              }
            }
          }
        }
      }
      
      
      
    } else if (check_stage && check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (stringcompare_hard(as<std::string>(data_stage2(year_guys(j))), as<std::string>(new_stage2(i)))) {
            if (data_agecol(year_guys(j)) == new_age(i)) {
              
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            }
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2 == new_year(i) - 1);
        int year_guys_length3 = year_guys3.n_elem;
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (stringcompare_hard(as<std::string>(data_stage3(year_guys3(j))), as<std::string>(new_stage2(i)))) {
            if (data_agecol(year_guys3(j)) == new_age(i) - 1) {
              
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            }
          }
        }
      }
    } else if (!check_stage && check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (data_agecol(year_guys(j)) == new_age(i)) {
            
            arma::uvec year_found = find(all_years == new_year(i));
            pop_by_years(year_found(0))++;
            
            new_frequency(i)++;
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2 == new_year(i) - 1);
        int year_guys_length3 = year_guys3.n_elem;
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (data_agecol(year_guys3(j)) == new_age(i) - 1) {
            
            arma::uvec year_found = find(all_years == new_year(i));
            pop_by_years(year_found(0))++;
            
            new_frequency(i)++;
          }
        }
      }
    }
  }
  
  /*
  Rcout << "Post-loop\n";
  Rcout << "\n Length of new_row_id: " << new_row_id.length() << "\n";
  Rcout << " First entry: " << new_row_id(0) << "\n";
  
  Rcout << "\n Length of new_stage: " << new_stage.length() << "\n";
  Rcout << " First entry: " << new_stage(0) << "\n";
  
  Rcout << "\n Length of new_age: " << new_age.length() << "\n";
  Rcout << " First entry: " << new_age(0) << "\n";
  
  Rcout << "\n Length of new_year: " << new_year.length() << "\n";
  Rcout << " First entry: " << new_year(0) << "\n";
  
  Rcout << "\n Length of new_frequency: " << new_frequency.length() << "\n";
  Rcout << " First entry: " << new_frequency(0) << "\n";
  
  Rcout << "\n Length of new_actualprop: " << new_actualprop.length() << "\n";
  Rcout << " First entry: " << new_actualprop(0) << "\n";
  */
  
  for (int i = 0; i < new_df_rows; i++) {
    arma::uvec year_found = find(all_years == new_year(i));
    
    if (pop_by_years(year_found(0)) > 0) {
      new_actualprop(i) = new_frequency(i) / static_cast<double>(pop_by_years(year_found(0)));
    } else {
      new_actualprop(i) = 0;
    }
  }
  
  
  // Structure the output data frame
  List output;
  
  if (check_stage && !check_age) {
    List raw_output (8);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_stageindex;
    raw_output(2) = new_stage;
    raw_output(3) = new_stage2;
    raw_output(4) = new_stage1;
    raw_output(5) = new_year;
    raw_output(6) = new_frequency;
    raw_output(7) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "stageindex", "stage", "stage2",
      "stage1", "year2", "Freq", "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
    
  } else if (!check_stage && check_age) {
    List raw_output (5);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_age;
    raw_output(2) = new_year;
    raw_output(3) = new_frequency;
    raw_output(4) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "age", "year2", "Freq",
      "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
    
  } else if (check_stage && check_age) {
    List raw_output (9);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_stageindex;
    raw_output(2) = new_stage;
    raw_output(3) = new_stage2;
    raw_output(4) = new_stage1;
    raw_output(5) = new_age;
    raw_output(6) = new_year;
    raw_output(7) = new_frequency;
    raw_output(8) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "stageindex", "stage", "stage2",
      "stage1", "age", "year2", "Freq", "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
  
  }
  
  output.attr("row.names") = new_row_names;
  output.attr("class") = "data.frame";
  return output;
}

//' Check and Reorganize Density Input Table Into Usable Format
//' 
//' Function \code{density_reassess()} takes a density input table as supplied
//' by the \code{\link{density_input}()} function, and checks and rearranges it
//' into a single, complete density input table.
//' 
//' @name density_reassess
//' 
//' @param stageframe The correct stageframe, already modified by
//' \code{\link{sf_reassess}()}.
//' @param dens_inp The density input data frame as is toward the end of
//' \code{\link{density_input}()}.
//' @param agestages The agestages element from the used \code{lefkoMat} object.
//' Only needed if an age-by-stage MPM will be used.
//' @param historical A logical value denoting whether MPM is historical.
//' Defaults to \code{FALSE}.
//' @param agebystage A logical value denoting whether MPM is age-by-stage.
//' Defaults to \code{FALSE}.
//' 
//' @return A corrected density input deta frame, usable in density-dependent
//' MPM creation.
//' 
//' @keywords internal
//' @noRd
Rcpp::DataFrame density_reassess(DataFrame stageframe, DataFrame dens_inp,
  Nullable<DataFrame> agestages, bool historical = false,
  bool agebystage = false) {
  
  StringVector stagevec = as<StringVector>(stageframe["stage"]);
  arma::ivec stageidvec = as<arma::ivec>(stageframe["stage_id"]);
  NumericVector minagevec = as<NumericVector>(stageframe["min_age"]);
  NumericVector maxagevec = as<NumericVector>(stageframe["max_age"]);
  arma::uvec repvec = as<arma::uvec>(stageframe["repstatus"]);
  arma::uvec obsvec = as<arma::uvec>(stageframe["obsstatus"]);
  arma::uvec propvec = as<arma::uvec>(stageframe["propstatus"]);
  arma::uvec immvec = as<arma::uvec>(stageframe["immstatus"]);
  arma::uvec matvec = as<arma::uvec>(stageframe["matstatus"]);
  arma::uvec indvec = as<arma::uvec>(stageframe["indataset"]);
  arma::ivec groupvec = as<arma::ivec>(stageframe["group"]);
  
  int no_stages = indvec.n_elem;
  arma::ivec alive(no_stages, fill::ones);

  // Identify all groups in the stageframe
  arma::ivec all_groups = unique(groupvec);
  int no_groups = all_groups.n_elem;
  StringVector group_text(no_groups);
  
  for (int i = 0; i < no_groups; i++) {
    group_text(i) = "group";
    group_text(i) += std::to_string(all_groups(i));
  }
  
  StringVector stage3_di = dens_inp["stage3"];
  StringVector stage2_di = dens_inp["stage2"];
  StringVector stage1_di = dens_inp["stage1"];
  IntegerVector age2_di = dens_inp["age2"];
  IntegerVector style_di = dens_inp["style"];
  IntegerVector time_delay_di = dens_inp["time_delay"];
  NumericVector alpha_di = dens_inp["alpha"];
  NumericVector beta_di = dens_inp["beta"];
  IntegerVector type_di = dens_inp["type"];
  IntegerVector type_t12_di = dens_inp["type_t12"];
  int di_rows = stage3_di.length();
  
  StringVector unique_stages = unique(stagevec);
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
  
  DataFrame agestages_;
  arma::ivec agestages_stageid;
  StringVector agestages_stage;
  arma::ivec agestages_age;
  int agestages_rows = 0;
  int age_min = 0;
  
  if (agestages.isNotNull()) {
    if (is<NumericVector>(agestages)) {
      if (agebystage) {
        Rf_warningcall(R_NilValue, "Function density_input() requires an agestages object if input MPM is age-by-stage.");
        agebystage = false;
      }
      
    } else if (is<DataFrame>(agestages)) {
      agestages_ = as<DataFrame>(agestages);
      
      if (agestages_.length() == 1) {
        agebystage = false;
        
      } else {
        agestages_stage = as<StringVector>(agestages_["stage"]);
        agestages_stageid = as<arma::ivec>(agestages_["stage_id"]);
        agestages_age = as<arma::ivec>(agestages_["age"]);
        agestages_rows = agestages_stage.length();
        age_min = agestages_age.min();
        
        agebystage = true;
      }
    } else {
      throw Rcpp::exception("Object input as agestages is not recognized.", false);
    }
  }
  
  if (historical && agebystage) {
    throw Rcpp::exception("MPMs cannot be both historical and age-by-stage.", false);
  }
  
  // Check for good entries density input data frame
  for (int i = 0; i < stage3_di.length(); i++) {
    int s3di_count {0};
    int s2di_count {0};
    int s1di_count {0};
    
    for (int j = 0; j < all_possible_stage_terms.length(); j++) {
      std::string s3used = as<std::string>(stage3_di(i));
      std::string s2used = as<std::string>(stage2_di(i));
      std::string s1used = as<std::string>(stage1_di(i));
      
      for (int k = 0; k < s3used.size(); k++) {
        s3used[k] = tolower(s3used[k]);
      }
      for (int k = 0; k < s3used.size(); k++) {
        s2used[k] = tolower(s2used[k]);
      }
      for (int k = 0; k < s3used.size(); k++) {
        s1used[k] = tolower(s1used[k]);
      }
      
      if (stringcompare_hard(s3used, as<std::string>(all_possible_stage_terms(j)))) s3di_count++;
      if (stringcompare_hard(as<std::string>(stage3_di(i)), as<std::string>(all_possible_stage_terms(j)))) s3di_count++;
      
      if (stringcompare_hard(s2used, as<std::string>(all_possible_stage_terms(j)))) s2di_count++;
      if (stringcompare_hard(as<std::string>(stage2_di(i)), as<std::string>(all_possible_stage_terms(j)))) s2di_count++;
      
      if (stringcompare_hard(s3used, "notalive")) {
        throw Rcpp::exception("Stage NotAlive is not allowed.", false);
      }
      if (stringcompare_hard(s2used, "notalive")) {
        throw Rcpp::exception("Stage NotAlive is not allowed.", false);
      }
      
      if (historical) {
        if (stringcompare_hard(s1used, as<std::string>(all_possible_stage_terms(j)))) s1di_count++;
        if (stringcompare_hard(as<std::string>(stage1_di(i)), as<std::string>(all_possible_stage_terms(j)))) s1di_count++;
        
        if (stringcompare_hard(s1used, "notalive")) {
          throw Rcpp::exception("Stage NotAlive is not allowed.", false);
        }
      } 
    }
    
    if (s3di_count == 0) {
      throw Rcpp::exception("Stage names in density input frame (stage3) must match stageframe",
        false);
    }
    if (s2di_count == 0) {
      throw Rcpp::exception("Stage names in density input frame (stage2) must match stageframe",
        false);
    }
    if (historical) {
      if (s1di_count == 0) {
        throw Rcpp::exception("Stage names in density input frame (stage1) must match stageframe",
          false);
      }
    }
  }
    
  IntegerVector s1_calls (di_rows, 1);
  IntegerVector s2_calls (di_rows, 1);
  IntegerVector s3_calls (di_rows, 1);
  IntegerVector s3_planned (di_rows, 1);
  IntegerVector s2_planned (di_rows, 1);
  IntegerVector s1_planned (di_rows, 1);
  
  IntegerVector s123_calls (di_rows, 1);
  
  List age3_calls (di_rows);
  List age2_calls (di_rows);
  
  List stageid3_calls (di_rows);
  List stageid2_calls (di_rows);
  
  arma::uvec prop_stages = find(propvec);
  arma::uvec prop0_stages = find(propvec == 0);
  arma::uvec imm_stages = find(immvec);
  arma::uvec alive_stages = find(alive);
  arma::uvec mat_stages = find(matvec);
  arma::uvec rep_stages = find(repvec);
  arma::uvec rep0_stages = find(repvec == 0);
  arma::uvec mat_rep0_stages = intersect(mat_stages, rep0_stages);
  arma::uvec obs_stages = find(obsvec);
  arma::uvec obs0_stages = find(obsvec == 0);
  arma::uvec all_stages = find(alive); // 7 "all"
  int no_current_group {0};
  
  // Now we build the expanded and edited density input frame
  for (int i = 0; i < di_rows; i++) {
    
    std::string s3used = as<std::string>(stage3_di(i));
    std::string s2used = as<std::string>(stage2_di(i));
    std::string s1used = as<std::string>(stage1_di(i));
    
    for (int j = 0; j < s3used.size(); j++) {
      s3used[j] = tolower(s3used[j]);
    }
    for (int j = 0; j < s2used.size(); j++) {
      s2used[j] = tolower(s2used[j]);
    }
    for (int j = 0; j < s1used.size(); j++) {
      s1used[j] = tolower(s1used[j]);
    }
    
    // Time t+1
    if (stringcompare_hard(s3used, "prop")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = prop_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "npr")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
         s3_calls(i) = prop0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "immat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = imm_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "mat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
         s3_calls(i) = mat_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "rep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = rep_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "nrep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = mat_rep0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "obs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = obs_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "nobs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = obs0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "all")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = all_stages.n_elem;
      }
      
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage3_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = current_group.n_elem;
          
          if (agebystage) {
            int found_stages = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (type_di(i) == 1) {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (IntegerVector::is_na(age2_di(i))) {
                      found_stages++;
                    } else if (agestages_age(k) == (age2_di(i) + 1)) {
                      found_stages++;
                    }
                  }
                } else {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (agestages_age(k) == age_min) found_stages++;
                  }
                }
              }
            }
            
            s3_calls(i) = found_stages;
            IntegerVector age3_vec(found_stages);
            IntegerVector stageid3_vec(found_stages);
            int a3_counter = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (type_di(i) == 1) {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (IntegerVector::is_na(age2_di(i))) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    } else if (agestages_age(k) == (age2_di(i) + 1)) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    }
                  }
                } else {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (agestages_age(k) == age_min) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    }
                  }
                }
              }
            }
            age3_calls(i) = age3_vec;
            stageid3_calls(i) = stageid3_vec;
            
          } else {
            s3_calls(i) = no_current_group;
          }
        }
      }
    }
    if (s3_calls(i) == 0) s3_calls(i) = 1;
    
    // Time t
    if (stringcompare_hard(s2used, "prop")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = prop_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "npr")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = prop0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "immat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (imm_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (imm_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = imm_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "mat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = mat_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "rep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (rep_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (rep_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = rep_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "nrep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = mat_rep0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "obs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = obs_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "nobs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = obs0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "all")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (all_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (all_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = all_stages.n_elem;
      }
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage2_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = current_group.n_elem;
          
          if (agebystage) {
            int found_stages = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (current_group(j) == agestages_stageid(k) - 1) {
                  if (IntegerVector::is_na(age2_di(i))) {
                    found_stages++;
                  } else if (agestages_age(k) == age2_di(i)) {
                    found_stages++;
                  }
                }
              }
            }
            
            s2_calls(i) = found_stages;
            IntegerVector age2_vec(found_stages);
            IntegerVector stageid2_vec(found_stages);
            int a2_counter = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (current_group(j) == agestages_stageid(k) - 1) {
                  if (IntegerVector::is_na(age2_di(i))) {
                    age2_vec(a2_counter) = agestages_age(k);
                    stageid2_vec(a2_counter) = agestages_stageid(k);
                    a2_counter++;
                  } else if (agestages_age(k) == age2_di(i)) {
                    age2_vec(a2_counter) = agestages_age(k);
                    stageid2_vec(a2_counter) = agestages_stageid(k);
                    a2_counter++;
                  }
                }
              }
            }
            age2_calls(i) = age2_vec;
            stageid2_calls(i) = stageid2_vec;
            
          } else {
            s2_calls(i) = no_current_group;
          }
        }
      }
    }
    if (s2_calls(i) == 0) s2_calls(i) = 1;
    
    // Time t-1
    if (stringcompare_hard(s1used, "prop")) {
      s1_calls(i) = prop_stages.n_elem;
    } else if (stringcompare_hard(s1used, "npr")) {
      s1_calls(i) = prop0_stages.n_elem;
    } else if (stringcompare_hard(s1used, "immat")) {
      s1_calls(i) = imm_stages.n_elem;
    } else if (stringcompare_hard(s1used, "mat")) {
      s1_calls(i) = mat_stages.n_elem;
    } else if (stringcompare_hard(s1used, "rep")) {
      s1_calls(i) = rep_stages.n_elem;
    } else if (stringcompare_hard(s1used, "nrep")) {
      s1_calls(i) = mat_rep0_stages.n_elem;
    } else if (stringcompare_hard(s1used, "obs")) {
      s1_calls(i) = obs_stages.n_elem;
    } else if (stringcompare_hard(s1used, "nobs")) {
      s1_calls(i) = obs0_stages.n_elem;
    } else if (stringcompare_hard(s1used, "all")) {
      s1_calls(i) = all_stages.n_elem;
    } else if (StringVector::is_na(stage1_di(i))) {
      s1_calls(i) = 1;
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage1_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = current_group.n_elem;
          
          s1_calls(i) = no_current_group;
        }
      }
    }
    if (s1_calls(i) == 0) s1_calls(i) = 1;
    
    s123_calls(i) = s3_calls(i) * s2_calls(i) * s1_calls(i);
  }
  
  NumericVector basepoints(di_rows, 0.0);
  for (int i = 0; i < (di_rows - 1); i++) {
    basepoints(i+1) = basepoints(i) + s123_calls(i);
  }
  
  // New output data frame set-up
  int newdi_rows = sum(s123_calls);
  
  StringVector stage3_newdi(newdi_rows);
  StringVector stage2_newdi(newdi_rows);
  StringVector stage1_newdi(newdi_rows);
  IntegerVector age2_newdi(newdi_rows);
  IntegerVector style_newdi(newdi_rows);
  NumericVector alpha_newdi(newdi_rows);
  NumericVector beta_newdi(newdi_rows);
  IntegerVector time_delay_newdi(newdi_rows);
  IntegerVector type_newdi(newdi_rows);
  IntegerVector type_t12_newdi(newdi_rows);
  
  int overall_counter {0};
  // int group_check {0};
  
  int group_baseline3 {0};
  int group_baseline2 {0};
  int group_baseline1 {0};
  
  int group_ratchet3 {0};
  int group_ratchet2 {0};
  int group_ratchet1 {0};
  
  int prevl3 {0};
  int prevl2 {0};
  int prevl1 {0};
  
  for (int i = 0; i < di_rows; i++) {
    overall_counter = 0;
    
    int age3_counter = 0;
    int age2_counter = 0;
    
    std::string s3used = as<std::string>(stage3_di(i));
    std::string s2used = as<std::string>(stage2_di(i));
    std::string s1used = as<std::string>(stage1_di(i));
    
    for (int j = 0; j < s3used.size(); j++) {
      s3used[j] = tolower(s3used[j]);
    }
    for (int j = 0; j < s2used.size(); j++) {
      s2used[j] = tolower(s2used[j]);
    }
    for (int j = 0; j < s1used.size(); j++) {
      s1used[j] = tolower(s1used[j]);
    }
    
    for (int j = 0; j < s1_calls(i); j++) {
      for (int k = 0; k < s2_calls(i); k++) {
        for (int l = 0; l < s3_calls(i); l++) {
          
          // Time t+1
          if (stringcompare_hard(s3used, "prop")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "npr")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "immat")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "mat")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "rep")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "nrep")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "obs")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "nobs")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "all")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(l));
            }
            
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage3_di(i) == group_text(m)) {
                if (agebystage) {
                  group_check = 1;
                  
                  IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
                  IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
                  int a3v_length = age3_vec.length();
                  
                  stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
                  
                  age3_counter++;
                  if (age3_counter == a3v_length) age3_counter = 0;
                
                } else {
                  if (l == 0) group_ratchet3 = 0;
                  if (l != prevl3 && l != 0) group_ratchet3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(groupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet3 > (current_group_length - 1)) {
                    group_ratchet3 = 0;
                  }
                  
                  if (group_ratchet3 == 0) {
                    group_baseline3 = l;
                  }
                  
                  stage3_newdi(basepoints(i) + overall_counter) = 
                    stagevec(current_group(l - group_baseline3));
                  
                  prevl3 = l;
                  
                }
              }
            }
            if (group_check == 0) {
              stage3_newdi(basepoints(i) + overall_counter) = stage3_di(i);
            }
            
            group_check = 0;
          }
          
          // Set up of most core variables in output data frame
          age2_newdi(basepoints(i) + overall_counter) = age2_di(i);
          style_newdi(basepoints(i) + overall_counter) = style_di(i);
          alpha_newdi(basepoints(i) + overall_counter) = alpha_di(i);
          beta_newdi(basepoints(i) + overall_counter) = beta_di(i);
          time_delay_newdi(basepoints(i) + overall_counter) = time_delay_di(i);
          type_newdi(basepoints(i) + overall_counter) = type_di(i);
          type_t12_newdi(basepoints(i) + overall_counter) = type_t12_di(i);
          
          // Time t
          if (stringcompare_hard(s2used, "prop")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(k));
            }
          } else if (stringcompare_hard(s2used, "npr")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "immat")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(k));
            }
          } else if (stringcompare_hard(s2used, "mat")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(k));
            }
          } else if (stringcompare_hard(s2used, "rep")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(k));
            }
          } else if (stringcompare_hard(s2used, "nrep")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "obs")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(k));
            }
          } else if (stringcompare_hard(s2used, "nobs")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "all")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(k));
            }
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage2_di(i) == group_text(m)) {
                if (agebystage) {
                  group_check = 1;
                  
                  IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
                  IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
                  int a2v_length = age2_vec.length();
                  
                  stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
                  age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
                  
                  age2_counter++;
                  if (age2_counter == a2v_length) age2_counter = 0;
                
                } else {
                  if (k == 0) group_ratchet2 = 0;
                  if (k != prevl2 && k != 0) group_ratchet2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(groupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet2 > (current_group_length - 1)) {
                    group_ratchet2 = 0;
                  }
                  
                  if (group_ratchet2 == 0) {
                    group_baseline2 = k;
                  }
                  
                  stage2_newdi(basepoints(i) + overall_counter) =
                    stagevec(current_group(k - group_baseline2));
                  
                  prevl2 = k;
                }
              }
            }
             
            if (group_check == 0) {
             stage2_newdi(basepoints(i) + overall_counter) = stage2_di(i);
            }
            
            group_check = 0;
          }
          
          // Time t-1
          if (stringcompare_hard(s1used, "prop")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(j));
          } else if (stringcompare_hard(s1used, "npr")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(j));
          } else if (stringcompare_hard(s1used, "immat")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(j));
          } else if (stringcompare_hard(s1used, "mat")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(j));
          } else if (stringcompare_hard(s1used, "rep")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(j));
          } else if (stringcompare_hard(s1used, "nrep")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(j));
          } else if (stringcompare_hard(s1used, "obs")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(j));
          } else if (stringcompare_hard(s1used, "nobs")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(j));
          } else if (stringcompare_hard(s1used, "all")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(j));
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage1_di(i) == group_text(m)) {
                if (j == 0) group_ratchet1 = 0;
                if (j != prevl1 && j != 0) group_ratchet1 += 1;
                
                group_check = 1;
                arma::uvec current_group = find(groupvec == m);
                int current_group_length = current_group.n_elem;
                if (group_ratchet1 > (current_group_length - 1)) {
                  group_ratchet1 = 0;
                }
                
                if (group_ratchet1 == 0) {
                  group_baseline1 = j;
                }
                
                stage1_newdi(basepoints(i) + overall_counter) =
                  stagevec(current_group(j - group_baseline1));
                
                prevl1 = j;
              }
            }
            
            if (group_check == 0) {
              stage1_newdi(basepoints(i) + overall_counter) = stage1_di(i);
            }
            
            group_check = 0;
          }
          
          overall_counter++;
        }
      }
    }
  }
  
  // Output final set-up
  Rcpp::List new_di(10);
  
  new_di(0) = stage3_newdi;
  new_di(1) = stage2_newdi;
  new_di(2) = stage1_newdi;
  new_di(3) = age2_newdi;
  new_di(4) = style_newdi;
  new_di(5) = alpha_newdi;
  new_di(6) = beta_newdi;
  new_di(7) = time_delay_newdi;
  new_di(8) = type_newdi;
  new_di(9) = type_t12_newdi;
  
  CharacterVector namevec = {"stage3", "stage2", "stage1", "age2", "style",
    "alpha", "beta", "time_delay", "type", "type_t12"};
  CharacterVector newclass = {"data.frame"};
  new_di.attr("names") = namevec;
  new_di.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newdi_rows);
  new_di.attr("class") = newclass;
  
  return new_di;
}

//' Create a Data Frame of Density Dependence Relationships in Matrix Elements
//' 
//' Function \code{density_input()} provides all necessary data to incorporate
//' density dependence into a \code{lefkoMat} object, a list of matrices, or a
//' single matrix. Four forms of density dependence are allowed, including the
//' Ricker function, the Beverton-Holt function, the Usher function, and the
//' logistic function. In each case, density must have an effect with at least a
//' one time-step delay (see Notes). The resulting data frame provides a guide
//' for other \code{lefko3} functions to modify matrix elements by density.
//'
//' @name density_input
//' 
//' @param mpm The \code{lefkoMat} object that will be subject to density
//' dependent projection.
//' @param stage3 A vector showing the name or number of the stage in occasion
//' \emph{t}+1 in the transitions to be affected by density. Abbreviations for
//' groups of stages are also usable (see Notes).
//' @param stage2 A vector showing the name or number of the stage in occasion
//' \emph{t} in the transition to be affected by density. Abbreviations for
//' groups of stages are also usable (see Notes).
//' @param stage1 A vector showing the name or number of the stage in occasion
//' \emph{t}-1 in the transition to be affected by density. Only needed if a
//' historical MPM is used. Abbreviations for groups of stages are also usable
//' (see Notes).
//' @param age2 A vector showing the age of the stage in occasion \emph{t} in the
//' transition to be affected by density. Only needed if an age-by-stage MPM is
//' used.
//' @param style A vector coding for the style of density dependence on each
//' transition subject to density dependence. Options include \code{1},
//' \code{ricker}, \code{ric}, or \code{r} for the Ricker function; \code{2},
//' \code{beverton}, \code{bev}, and \code{b} for the Beverton-Holt function;
//' \code{3}, \code{usher}, \code{ush}, and \code{u} for the Usher function; and
//' \code{4}, \code{logistic}, \code{log}, and \code{l} for the logistic
//' function. If only a single code is provided, then all noted transitions are
//' assumed to be subject to this style of density dependence. Defaults to
//' \code{ricker}.
//' @param time_delay An integer vector indicating the number of occasions back
//' on which density dependence operates. Defaults to \code{1}, and may not equal
//' any integer less than 1. If a single number is input, then all noted
//' transitions are assumed to be subject to this time delay.  Defaults to
//' \code{1}.
//' @param alpha A vector indicating the numeric values to use as the
//' alpha term in the two parameter Ricker, Beverton-Holt, or Usher function, or
//' the value of the carrying capacity \emph{K} to use in the logistic equation
//' (see \code{Notes} section for more on this term). If a single number is
//' provided, then all noted transitions are assumed to be subject to this value
//' of alpha. Defaults to \code{1}.
//' @param beta A vector indicating the numeric values to use as the beta term in
//' the two parameter Ricker, Beverton-Holt, or Usher function. Used to indicate
//' whether to use \emph{K} as a hard limit in the logistic equation (see section
//' \code{Notes} below). If a single number is provided, then all noted
//' transitions are assumed to be subject to this value of \code{beta}. Defaults
//' to \code{1}.
//' @param type A vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition; or 
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition. If empty or not provided, then defaults to \code{1} for survival
//' transition.
//' @param type_t12 An optional vector denoting the kind of transition between
//' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
//' deVries format is desired. This should be entered as \code{1}, \code{S}, or
//' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
//' fecundity transitions. Defaults to \code{1} for survival transition, with
//' impacts only on the construction of deVries-format hMPMs.
//' 
//' @return A data frame of class \code{lefkoDens}. This object can be used as
//' input in function \code{\link{projection3}()}.
//' 
//' Variables in this object include the following:
//' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
//' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
//' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced,
//' if applicable.}
//' \item{age2}{Age at occasion \emph{t} in the transition to be replaced, if
//' applicable.}
//' \item{style}{Style of density dependence, coded as 1, 2, 3, or 4 for the
//' Ricker, Beverton-Holt, Usher, or logistic function, respectively.}
//' \item{time_delay}{The time delay on density dependence, in time steps.}
//' \item{alpha}{The value of alpha in the Ricker, Beverton-Holt, or Usher
//' function, or the value of carrying capacity, \emph{K}, in the logistic
//' function.}
//' \item{beta}{The value of beta in the Ricker, Beverton-Holt, or Usher
//' function.}
//' \item{type}{Designates whether the transition from occasion \emph{t} to
//' occasion \emph{t}+1 is a survival transition probability (1), or a fecundity
//' rate (2).}
//' \item{type_t12}{Designates whether the transition from occasion \emph{t}-1 to
//' occasion \emph{t} is a survival transition probability (1), a fecundity rate
//' (2).}
//' 
//' @section Notes:
//' This function provides inputs when density dependence is operationalized
//' directly on matrix elements. It can be used in both \code{projection3()} and
//' \code{f_projection3()}. Users wishing to modify vital rate functions by
//' density dependence functions for use in function-based projections with
//' function \code{f_projection3()} should use function \code{density_vr()} to
//' provide the correct inputs.
//' 
//' The parameters \code{alpha} and \code{beta} are applied according to the
//' two-parameter Ricker function, the two-parameter Beverton-Holt function, the
//' two-parameter Usher function, or the one-parameter logistic function.
//' Although the default is that a 1 time step delay is assumed, greater time
//' delays can be set through the \code{time_delay} option.
//' 
//' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
//' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
//' are to be used, \code{nrep} if all mature but non-reproductive stages are to
//' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
//' immature stages are to be used, \code{prop} if all propagule stages are to be
//' used, \code{npr} if all non-propagule stages are to be used, \code{obs} if
//' all observable stages are to be used, \code{nobs} if all unobservable stages
//' are to be used, and leave empty or use \code{all} if all stages in stageframe
//' are to be used.
//' 
//' When using the logistic function, it is possible that the time delay used in
//' density dependent simulations will cause matrix elements to become negative.
//' To prevent this behavior, set the associated \code{beta} term to \code{1.0}.
//' Doing so will set \code{K} as the hard limit in the logistic equation,
//' essentially setting a minimum limit at \code{0} for all matrix elements
//' modified.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{projection3}()}
//' 
//' @examples
//' \donttest{
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
//' 
//' e3d <- density_input(ehrlen3mean, stage3 = c("Sd", "Sdl"),
//'   stage2 = c("rep", "rep"), stage1 = c("all", "all"), style = 1,
//'   time_delay = 1, alpha = 1, beta = 0, type = c(2, 2), type_t12 = c(1, 1))
//' 
//' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE, substoch = 2,
//'   density = e3d)
//' }
//' 
//' @export density_input
// [[Rcpp::export(density_input)]]
DataFrame density_input(List mpm, RObject stage3, RObject stage2,
  Nullable<RObject> stage1 = R_NilValue, Nullable<RObject> age2 = R_NilValue,
  Nullable<RObject> style = R_NilValue, Nullable<RObject> time_delay = R_NilValue,
  Nullable<RObject> alpha = R_NilValue, Nullable<RObject> beta = R_NilValue,
  Nullable<RObject> type = R_NilValue, Nullable<RObject> type_t12 = R_NilValue) {
  
  bool historical = false;
  bool agebystage = false;
  
  // Check quality of mpm input
  StringVector mpm_class_vec;
  if (mpm.hasAttribute("class")) {
    mpm_class_vec = mpm.attr("class");
  } else mpm_class_vec = {"list"};
  std::string mpm_class = as<std::string>(mpm_class_vec(0));
  
  CharacterVector mpm_elems = mpm.names();
  int mpm_name_check = 0;
  for (int i = 0; i < mpm_elems.length(); i++) {
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "ahstages")) mpm_name_check++;
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "hstages")) mpm_name_check++;
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "agestages")) mpm_name_check++;
  }
  
  if (!(mpm_name_check == 3 && stringcompare_hard(mpm_class, "lefkoMat"))) {
    throw Rcpp::exception("This function requires a lefkoMat object as input.", false);
  }
  
  DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
  DataFrame hstages = as<DataFrame>(mpm["hstages"]);
  DataFrame agestages = as<DataFrame>(mpm["agestages"]);
  
  if (hstages.length() > 1) {
    historical = true;
  }
  if (agestages.length() > 1) {
    agebystage = true;
  } 
  if (stageframe.length() == 1) {
    throw Rcpp::exception("Input lefkoMat object does not appear to have a stageframe, which should be element ahstages.",
      false);
  }
  
  CharacterVector ahstages_elems = stageframe.names();
  int ahs_name_check = 0;
  for (int i = 0; i < ahstages_elems.length(); i++) {
    if (stringcompare_hard(as<std::string>(ahstages_elems(i)), "stage_id")) ahs_name_check++;
    if (stringcompare_hard(as<std::string>(ahstages_elems(i)), "stage")) ahs_name_check++;
  }
  
  if (ahs_name_check != 2) {
    throw Rcpp::exception("Stageframe appears to be modified. Please make sure that element ahstages in the lefkoMat object includes both a stage column holding stage names and a stage_id column holding unique, stage identifying integers.", 
      false);
  }
  
  IntegerVector stage_id_sf = as<IntegerVector>(stageframe["stage_id"]);
  StringVector stage_sf = as<StringVector>(stageframe["stage"]);
  int no_stages = stage_sf.length();
  
  // Density input vector standardization
  StringVector stage3_names;
  StringVector stage2_names;
  StringVector stage1_names;
  
  if (is<StringVector>(stage3)) {
    StringVector stage3_names_ = as<StringVector>(stage3);
    stage3_names = stage3_names_;
    
  } else if (is<IntegerVector>(stage3)) {
    arma::ivec stage3_ids = as<arma::ivec>(stage3);
    int stage3_entries = stage3_ids.n_elem;
    
    arma::uvec bad_lows = find(stage3_ids < 1);
    arma::uvec bad_highs = find(stage3_ids > no_stages);
    
    if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
      throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
    }
    
    StringVector new_stage3 (stage3_entries);
    for (int i = 0; i < stage3_entries; i++) {
      new_stage3(i) = stage_sf(stage3_ids(i) - 1);
    }
    
    stage3_names = new_stage3;
    
  } else {
    throw Rcpp::exception("Option stage3 is not a recognized input type.", false);
  }
  
  if (is<StringVector>(stage2)) {
    StringVector stage2_names_ = as<StringVector>(stage2);
    stage2_names = stage2_names_;
    
  } else if (is<IntegerVector>(stage2)) {
    arma::ivec stage2_ids = as<arma::ivec>(stage2);
    int stage2_entries = stage2_ids.n_elem;
    
    arma::uvec bad_lows = find(stage2_ids < 1);
    arma::uvec bad_highs = find(stage2_ids > no_stages);
    
    if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
      throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
    }
    
    StringVector new_stage2 (stage2_entries);
    for (int i = 0; i < stage2_entries; i++) {
      new_stage2(i) = stage_sf(stage2_ids(i) - 1);
    }
    
    stage2_names = new_stage2;
    
  } else {
    throw Rcpp::exception("Option stage2 is not a recognized input type.", false);
  }
  
  if (stage1.isNotNull()) {
    if (is<IntegerVector>(stage1)) {
      arma::ivec stage1_ids = as<arma::ivec>(stage1);
      int stage1_entries = stage1_ids.n_elem;
      
      arma::uvec bad_lows = find(stage1_ids < 1);
      arma::uvec bad_highs = find(stage1_ids > no_stages);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
      }
      
      StringVector new_stage1 (stage1_entries);
      for (int i = 0; i < stage1_entries; i++) {
        new_stage1(i) = stage_sf(stage1_ids(i) - 1);
      }
      
      stage1_names = new_stage1;
      
    } else if (is<StringVector>(stage1)) {
      StringVector stage1_names_ = as<StringVector>(stage1);
      stage1_names = stage1_names_;
      
    } else {
      throw Rcpp::exception("Option stage1 is not a recognized input type.", false);
    }
    
    if (!historical) {
      throw Rcpp::exception("Do not include option stage1 for ahistorical MPMs.", false);
    }
    
  } else if (historical) {
    throw Rcpp::exception("Option stage1 is needed for historical MPMs.", false);
    
  } else {
    StringVector stage1_names_ = {NA_STRING};
    stage1_names = stage1_names_;
  }
  
  if (stage3_names.length() != stage2_names.length()) {
    throw Rcpp::exception("All transitions to modify require information at least for stage2 and stage3. These inputs must also be of equal length.", 
      false);
  }
  
  if (historical) {
    if (stage3_names.length() != stage1_names.length()) {
      throw Rcpp::exception("All historical transitions to modify require information for stage1, stage2, and stage3. These inputs must also be of equal length.", 
        false);
    }
  }
  
  IntegerVector age2_vec;
  IntegerVector style_vec;
  IntegerVector time_delay_vec;
  NumericVector alpha_vec;
  NumericVector beta_vec;
  IntegerVector type_vec;
  IntegerVector type_t12_vec;
  
  if (age2.isNotNull()) {
    if (!agebystage) {
      throw Rcpp::exception("Do not use the age2 option unless the MPM is age-by-stage.", false);
    }
    
    if (is<IntegerVector>(age2)) {
      IntegerVector age2_vec_ = as<IntegerVector>(age2);
      age2_vec = age2_vec_;
      
      arma::ivec age2_arma = as<arma::ivec>(age2_vec_);
      arma::uvec bad_lows = find(age2_arma < 0);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Negative ages are not allowed.", false);
      }
    } else {
      throw Rcpp::exception("Only integers are allowed in option age2.", false);
    }
  } else {
    IntegerVector age2_vec_ = {NA_INTEGER};
    age2_vec = age2_vec_;
  }
  
  if (style.isNotNull()) {
    StringVector ricker_style = {"1", "ricker", "ricke", "rick", "ric", "ri", "r"};
    StringVector beverton_style = {"2", "beverton", "beverto", "bevert", "bever",
      "beve", "bev", "be", "b", "holt", "hol", "ho", "h"};
    StringVector usher_style = {"3", "usher", "ushe", "ush", "us", "u"};
    StringVector logistic_style = {"4", "logistic", "logisti", "logist", "logis",
      "logi", "log", "lo", "l"};
    
    if (is<IntegerVector>(style)) {
      IntegerVector style_ = as<IntegerVector>(style);
      style_vec = style_;
      
      arma::ivec style_arma = as<arma::ivec>(style_);
      
      arma::uvec bad_lows = find(style_arma < 1);
      arma::uvec bad_highs = find(style_arma > 4);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Invalid density dependence style entered.", false);
      }
      
    } else if (is<NumericVector>(style)) {
      IntegerVector style_ = as<IntegerVector>(style);
      style_vec = style_;
      
      arma::ivec style_arma = as<arma::ivec>(style_);
      
      arma::uvec bad_lows = find(style_arma < 1);
      arma::uvec bad_highs = find(style_arma > 4);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Invalid density dependence style entered.", false);
      }
      
    } else if (is<StringVector>(style)) {
      StringVector style_stringvec = as<StringVector>(style);
      int style_elems = style_stringvec.length();
      
      IntegerVector style_vec_ (style_elems, 0);
      
      for (int i = 0; i < style_elems; i++) {
        std::string ssv = as<std::string>(style_stringvec(i));
        for (int j = 0; j < ssv.size(); j++) {
          ssv[j] = tolower(ssv[j]);
        }
        
        for (int j = 0; j < ricker_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(ricker_style(j)))) {
            style_vec_(i) = 1;
          }
        }
        
        for (int j = 0; j < beverton_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(beverton_style(j)))) {
            style_vec_(i) = 2;
          }
        }
        
        for (int j = 0; j < usher_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(usher_style(j)))) {
            style_vec_(i) = 3;
          }
        }
        
        for (int j = 0; j < logistic_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(logistic_style(j)))) {
            style_vec_(i) = 1;
          }
        }
      }
      
      if (min(style_vec_) == 0) {
        throw Rcpp::exception("Some density dependence styles were not recognized.", false);
      }
      
      style_vec = style_vec_;
      
    } else {
      throw Rcpp::exception("Invalid entry for option style.", false);
      
    }
  } else {
    IntegerVector style_vec_ = {1};
    style_vec = style_vec_;
  }
  
  if (time_delay.isNotNull()) {
    if (is<IntegerVector>(time_delay)) {
      IntegerVector time_delay_vec_ = as<IntegerVector>(time_delay);
      time_delay_vec = time_delay_vec_;
      
      arma::ivec time_delay_arma = as<arma::ivec>(time_delay_vec_);
      arma::uvec bad_lows = find(time_delay_arma < 1);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Time delays less than 1 time step are not allowed.", false);
      }
    } else if (is<NumericVector>(time_delay)) {
      IntegerVector time_delay_vec_ = as<IntegerVector>(time_delay);
      time_delay_vec = time_delay_vec_;
      
      arma::ivec time_delay_arma = as<arma::ivec>(time_delay_vec_);
      arma::uvec bad_lows = find(time_delay_arma < 1);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Time delays less than 1 time step are not allowed.", false);
      }
    } else {
      throw Rcpp::exception("Only integers are allowed in option time_delay.", false);
    }
  } else {
    IntegerVector time_delay_vec_ = {1};
    time_delay_vec = time_delay_vec_;
  }
  
  if (alpha.isNotNull()) {
    if (is<NumericVector>(alpha)) {
      NumericVector alpha_vec_ = as<NumericVector>(alpha);
      alpha_vec = alpha_vec_;
      
    } else {
      throw Rcpp::exception("Only floating point decimals are allowed in option alpha.", false);
    }
    
  } else {
    NumericVector alpha_vec_ = {1};
    alpha_vec = alpha_vec_;
  }
  
  if (beta.isNotNull()) {
    if (is<NumericVector>(beta)) {
      NumericVector beta_vec_ = as<NumericVector>(beta);
      beta_vec = beta_vec_;
      
    } else {
      throw Rcpp::exception("Only floating point decimals are allowed in option beta.", false);
    }
    
  } else {
    NumericVector beta_vec_ = {1};
    beta_vec = beta_vec_;
  }
  
  if (type.isNotNull()) {
    if (is<IntegerVector>(type)) {
      IntegerVector type_vec_ = as<IntegerVector>(type);
      type_vec = type_vec_;
      
      arma::ivec type_arma = as<arma::ivec>(type_vec_);
      arma::uvec bad_lows = find(type_arma < 1);
      arma::uvec bad_highs = find(type_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Transition types may only be type 1 or 2.", false);
      }
    } else if (is<StringVector>(type)) {
      StringVector type_stringvec = as<StringVector>(type);
      int type_elems = type_stringvec.length();
      
      IntegerVector type_vec_ (type_elems, 0);
      
      for (int i = 0; i < type_elems; i++) {
        if (stringcompare_simple(as<std::string>(type_stringvec(i)), "s", true)) {
          type_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "1", true)) {
          type_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "f", true)) {
          type_vec_(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "2", true)) {
          type_vec_(i) = 2;
        } else {
          throw Rcpp::exception("Invalid entry in option type.", false);
        }
      }
      type_vec = type_vec_;
      
    } else if (is<NumericVector>(type)) {
      IntegerVector type_vec_ = as<IntegerVector>(type);
      type_vec = type_vec_;
      
      arma::ivec type_arma = as<arma::ivec>(type_vec_);
      arma::uvec bad_lows = find(type_arma < 1);
      arma::uvec bad_highs = find(type_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Transition types may only be type 1 or 2.", false);
      }
      
    } else {
      throw Rcpp::exception("Only integers are allowed in option type.", false);
    }
    
  } else {
    IntegerVector type_vec_ = {1};
    type_vec = type_vec_;
  }
  
  if (type_t12.isNotNull()) {
    if (is<IntegerVector>(type_t12)) {
      IntegerVector type_t12_vec_ = as<IntegerVector>(type_t12);
      type_t12_vec = type_t12_vec_;
      
      arma::ivec type_t12_arma = as<arma::ivec>(type_t12_vec_);
      arma::uvec bad_lows = find(type_t12_arma < 1);
      arma::uvec bad_highs = find(type_t12_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Historical transition types may only be type 1 or 2.", false);
      }
    } else if (is<NumericVector>(type_t12)) {
      IntegerVector type_t12_vec_ = as<IntegerVector>(type_t12);
      type_t12_vec = type_t12_vec_;
      
      arma::ivec type_t12_arma = as<arma::ivec>(type_t12_vec_);
      arma::uvec bad_lows = find(type_t12_arma < 1);
      arma::uvec bad_highs = find(type_t12_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Historical transition types may only be type 1 or 2.", false);
      }
      
    } else if (is<StringVector>(type_t12)) {
      StringVector type_t12_stringvec = as<StringVector>(type_t12);
      int type_t12_elems = type_t12_stringvec.length();
      
      IntegerVector type_t12_vec_ (type_t12_elems, 0);
      
      for (int i = 0; i < type_t12_elems; i++) {
        if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "s", true)) {
          type_t12_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "1", true)) {
          type_t12_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "f", true)) {
          type_t12_vec_(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "2", true)) {
          type_t12_vec_(i) = 2;
        } else {
          throw Rcpp::exception("Invalid entry in option type_t12.", false);
        }
      }
      type_t12_vec = type_t12_vec_;
      
    } else {
      throw Rcpp::exception("Only integers are allowed in option type_t12.", false);
    }
    
  } else {
    IntegerVector type_t12_vec_ = {1};
    type_t12_vec = type_t12_vec_;
  }
  
  
  StringVector new_stage3_names;
  StringVector new_stage2_names;
  StringVector new_stage1_names;
  IntegerVector new_age2_vec;
  IntegerVector new_style_vec;
  IntegerVector new_time_delay_vec;
  NumericVector new_alpha_vec;
  NumericVector new_beta_vec;
  IntegerVector new_type_vec;
  IntegerVector new_type_t12_vec;
  
  IntegerVector vec_lengths = {static_cast<int>(stage3_names.length()),
    static_cast<int>(stage2_names.length()), static_cast<int>(stage1_names.length()),
    static_cast<int>(age2_vec.length()), static_cast<int>(style_vec.length()),
    static_cast<int>(time_delay_vec.length()), static_cast<int>(alpha_vec.length()),
    static_cast<int>(beta_vec.length()), static_cast<int>(type_vec.length()),
    static_cast<int>(type_t12_vec.length())};
  int vec_max = max(vec_lengths);
  
  if (stage3_names.length() != vec_max) {
    if (stage3_names.length() == 1) {
      StringVector ns3n (vec_max, stage3_names(0));
      new_stage3_names = ns3n;
      
    } else {
      throw Rcpp::exception("Vector stage3 is not the correct length.", false);
    }
    
  } else {
    new_stage3_names = stage3_names;
  }
  
  if (stage2_names.length() != vec_max) {
    if (stage2_names.length() == 1) {
      StringVector ns2n (vec_max, stage2_names(0));
      new_stage2_names = ns2n;
      
    } else {
      throw Rcpp::exception("Vector stage2 is not the correct length.", false);
    }
    
  } else {
    new_stage2_names = stage2_names;
  }
  
  if (stage1_names.length() != vec_max) {
    if (stage1_names.length() == 1) {
      StringVector ns1n (vec_max, stage1_names(0));
      new_stage1_names = ns1n;
      
    } else {
      throw Rcpp::exception("Vector stage1 is not the correct length.", false);
    }
    
  } else {
    new_stage1_names = stage1_names;
  }
  
  if (age2_vec.length() != vec_max) {
    if (age2_vec.length() == 1) {
      IntegerVector na2v (vec_max, age2_vec(0));
      new_age2_vec = na2v;
      
    } else {
      throw Rcpp::exception("Vector age2 is not the correct length.", false);
    }
    
  } else {
    new_age2_vec = age2_vec;
  }
  
  if (style_vec.length() != vec_max) {
    if (style_vec.length() == 1) {
      IntegerVector ns2v (vec_max, style_vec(0));
      new_style_vec = ns2v;
      
    } else {
      throw Rcpp::exception("Vector style is not the correct length.", false);
    }
    
  } else {
    new_style_vec = style_vec;
  }
  
  if (time_delay_vec.length() != vec_max) {
    if (time_delay_vec.length() == 1) {
      IntegerVector ntd2v (vec_max, time_delay_vec(0));
      new_time_delay_vec = ntd2v;
      
    } else {
      throw Rcpp::exception("Vector time_delay is not the correct length.", false);
    }
    
  } else {
    new_time_delay_vec = time_delay_vec;
  }
  
  if (type_vec.length() != vec_max) {
    if (type_vec.length() == 1) {
      IntegerVector nt2v (vec_max, type_vec(0));
      new_type_vec = nt2v;
      
    } else {
      throw Rcpp::exception("Vector type is not the correct length.", false);
    }
    
  } else {
    new_type_vec = type_vec;
  }
  
  if (type_t12_vec.length() != vec_max) {
    if (type_t12_vec.length() == 1) {
      IntegerVector ntt122v (vec_max, type_t12_vec(0));
      new_type_t12_vec = ntt122v;
      
    } else {
      throw Rcpp::exception("Vector type_t12 is not the correct length.", false);
    }
    
  } else {
    new_type_t12_vec = type_t12_vec;
  }
  
  if (alpha_vec.length() != vec_max) {
    if (alpha_vec.length() == 1) {
      NumericVector nav (vec_max, alpha_vec(0));
      new_alpha_vec = nav;
      
    } else {
      throw Rcpp::exception("Vector alpha is not the correct length.", false);
    }
    
  } else {
    new_alpha_vec = alpha_vec;
  }
  
  if (beta_vec.length() != vec_max) {
    if (beta_vec.length() == 1) {
      NumericVector nbv (vec_max, beta_vec(0));
      new_beta_vec = nbv;
      
    } else {
      throw Rcpp::exception("Vector beta is not the correct length.", false);
    }
    
  } else {
    new_beta_vec = beta_vec;
  }
  
  DataFrame output_tab = DataFrame::create(_["stage3"] = new_stage3_names,
    _["stage2"] = new_stage2_names, _["stage1"] = new_stage1_names,
    _["age2"] = new_age2_vec, _["style"] = new_style_vec,
    _["time_delay"] = new_time_delay_vec, _["alpha"] = new_alpha_vec,
    _["beta"] = new_beta_vec, _["type"] = new_type_vec,
    _["type_t12"] = new_type_t12_vec);
  
  DataFrame output = density_reassess(stageframe, output_tab, agestages,
    historical, agebystage);
  
  StringVector needed_classes {"data.frame", "lefkoDens"};
  output.attr("class") = needed_classes;

  StringVector stage3_final = as<StringVector>(output["stage3"]);
  StringVector stage2_final = as<StringVector>(output["stage2"]);
  StringVector stage1_final = as<StringVector>(output["stage1"]);
  StringVector age2_final = as<StringVector>(output["age2"]);
  
  int final_rows = stage3_final.length();
  
  StringVector check_elems(final_rows);
  
  for (int i = 0; i < final_rows; i++) {
    check_elems(i) = stage3_final(i);
    check_elems(i) += " ";
    check_elems(i) += stage2_final(i);
    check_elems(i) += " ";
    if (!StringVector::is_na(stage1_final(i))) {
      check_elems(i) += stage1_final(i);
      check_elems(i) += " ";
    }
    if (!StringVector::is_na(age2_final(i))) check_elems(i) += age2_final(i);
  }
  
  StringVector unique_elems = unique(check_elems);
  if (unique_elems.length() != final_rows) {
    Rf_warningcall(R_NilValue, "Some transitions appear to be listed multiple times. This may cause errors in analysis.");
  }
  
  return output;
}

