#ifndef LEFKOUTILS_mat_stuff_H
#define LEFKOUTILS_mat_stuff_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

namespace LefkoMats {

  //' Create General Element Index for Any lefkoMat Matrix
  //' 
  //' This function creates a general element index by taking a list of matrices
  //' and finding all non-zero values, which it then return C++ indices of
  //' within an arma::uvec object.
  //' 
  //' @name general_index
  //' 
  //' @param mats A list of matrices.
  //' 
  //' @return An arma::uvec object listing indices of non-zero values in order.
  //' This object is essentially a union of all non-zero indices across all
  //' matrices in the list.
  //' 
  //' @keywords internal
  //' @noRd
  inline arma::uvec general_index (Rcpp::List mats) {
    int mat_length = mats.length();
    
    arma::uvec torture_chamber;
    
    for (int i = 0; i < mat_length; i++) {
      arma::uvec iron_maiden;
      if (is<S4>(mats(i))) {
        iron_maiden = find(as<arma::sp_mat>(mats(i)));
      } else {
        iron_maiden = find(as<arma::mat>(mats(i)));
      }
      
      if (i == 0) {
        torture_chamber = iron_maiden;
      } else {
        torture_chamber = unique(join_cols(torture_chamber, iron_maiden));
      }
    }
    
    return torture_chamber;
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
  inline Rcpp::List decomp3(arma::mat Amat) {
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
  inline Rcpp::List decomp3sp(arma::mat Amat) {
    arma::sp_mat spAmat(Amat);
    arma::sp_mat t_spAmat = spAmat.t();
    arma::cx_vec Aeigval;
    arma::cx_vec Aeigvall;
    arma::cx_mat Aeigvecl;
    arma::cx_mat Aeigvecr;
    
    eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
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
  inline Rcpp::List decomp3sp_inp(arma::sp_mat spAmat) {
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
  inline arma::mat ovreplace(arma::vec allst321, arma::vec idx321old,
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
  
  //' Creates Base Skeleton Stageframe
  //' 
  //' Function \code{sf_core()} creates a skeleton stageframe.
  //' 
  //' @name sf_core
  //' 
  //' @param num_stages An integer denoting the number of stages to include.
  //' @param reassessed A logical value indicating whether to create a
  //' processed stageframe, or an original stageframe.
  //' @param small A logical value indicating whether to create a small
  //' stageframe composed of only five variables, or a full stageframe.
  //' 
  //' @return A data frame with empty values for all variables except
  //' \code{"stage"}.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::DataFrame sf_core(int num_stages, bool reassessed = false,
    bool small = false) {
    
    if (num_stages < 1) {
      throw Rcpp::exception("Stageframe cannot be made for fewer than 1 stage.", false);
    }
    
    List output;
    
    StringVector stage (num_stages);
    NumericVector original_size (num_stages);
    IntegerVector entrystage (num_stages);
    IntegerVector repstatus (num_stages);
    NumericVector zero_double (num_stages);
    IntegerVector zero_int (num_stages);
    IntegerVector one_int (num_stages, 1);
    
    for (int i = 0; i < num_stages; i++) {
      String new_stage = "stage ";
      new_stage += String((i + 1));
      
      stage(i) = new_stage;
      
      original_size(i) = static_cast<double>(i + 1);
    }
    
    if (small) {
      IntegerVector stage_id (num_stages);
      for (int i = 0; i < num_stages; i++) {
        stage_id(i) = i+1;
      }
      
      List output_new (5);
      output_new(0) = stage_id;
      output_new(1) = stage;
      output_new(2) = original_size;
      output_new(3) = repstatus;
      output_new(4) = entrystage;
      
      output = output_new;
      Rcpp::CharacterVector varnames = {"stage_id", "stage", "original_size", 
        "repstatus", "entrystage"};
      output.attr("names") = varnames;
      
    } else {
      StringVector comments (num_stages);
      NumericVector bin_half (num_stages, 0.5);
      NumericVector bin_full (num_stages, 1.0);
      
      for (int i = 0; i < num_stages; i++) {
        String new_comment = "stage ";
        new_comment += String((i + 1));
        new_comment += " comment";
        
        comments(i) = new_comment;
      }
      
      if (reassessed) {
        IntegerVector stage_id (num_stages);
        for (int i = 0; i < num_stages; i++) {
          stage_id(i) = i+1;
        }
        
        List output_new (33);
        output_new(0) = stage_id;
        output_new(1) = stage;
        output_new(2) = original_size;
        output_new(3) = clone(zero_double);
        output_new(4) = clone(zero_double);
        output_new(5) = clone(zero_int);
        output_new(6) = clone(zero_int);
        output_new(7) = repstatus;
        output_new(8) = clone(one_int);
        output_new(9) = clone(zero_int);
        
        output_new(10) = clone(zero_int);
        output_new(11) = clone(zero_int);
        output_new(12) = entrystage;
        output_new(13) = clone(one_int);
        output_new(14) = bin_half;
        output_new(15) = clone(zero_double);
        output_new(16) = clone(zero_double);
        output_new(17) = clone(original_size);
        output_new(18) = bin_full;
        output_new(19) = clone(zero_double);
        
        output_new(20) = clone(zero_double);
        output_new(21) = clone(zero_double);
        output_new(22) = clone(zero_double);
        output_new(23) = clone(zero_double);
        output_new(24) = clone(zero_double);
        output_new(25) = clone(zero_double);
        output_new(26) = clone(zero_double);
        output_new(27) = clone(zero_double);
        output_new(28) = clone(zero_double);
        output_new(29) = clone(zero_int);
        
        output_new(30) = comments;
        output_new(31) = clone(one_int);
        output_new(32) = clone(zero_int);
        
        output = output_new;
        Rcpp::CharacterVector varnames = {"stage_id", "stage", "original_size",
          "original_size_b", "original_size_c", "min_age", "max_age",
          "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus",
          "entrystage", "indataset", "binhalfwidth_raw", "sizebin_min",
          "sizebin_max", "sizebin_center", "sizebin_width", "binhalfwidthb_raw",
          "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
          "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max",
          "sizebinc_center", "sizebinc_width", "group", "comments", "alive",
          "almostborn"};
        output.attr("names") = varnames;
        
      } else {
        List output_new (29);
        output_new(0) = stage;
        output_new(1) = original_size;
        output_new(2) = clone(zero_double);
        output_new(3) = clone(zero_double);
        output_new(4) = clone(zero_int);
        output_new(5) = clone(zero_int);
        output_new(6) = repstatus;
        output_new(7) = clone(one_int);
        output_new(8) = clone(zero_int);
        output_new(9) = clone(zero_int);
        
        output_new(10) = clone(zero_int);
        output_new(11) = clone(one_int);
        output_new(12) = bin_half;
        output_new(13) = clone(zero_double);
        output_new(14) = clone(zero_double);
        output_new(15) = clone(original_size);
        output_new(16) = bin_full;
        output_new(17) = clone(zero_double);
        output_new(18) = clone(zero_double);
        output_new(19) = clone(zero_double);
        
        output_new(20) = clone(zero_double);
        output_new(21) = clone(zero_double);
        output_new(22) = clone(zero_double);
        output_new(23) = clone(zero_double);
        output_new(24) = clone(zero_double);
        output_new(25) = clone(zero_double);
        output_new(26) = clone(zero_double);
        output_new(27) = clone(zero_int);
        output_new(28) = comments;
        
        output = output_new;
        Rcpp::CharacterVector varnames = {"stage", "size", "size_b", "size_c",
          "min_age", "max_age", "repstatus", "obsstatus", "propstatus",
          "immstatus", "matstatus", "indataset", "binhalfwidth_raw",
          "sizebin_min", "sizebin_max", "sizebin_center", "sizebin_width",
          "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
          "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw",
          "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
          "group", "comments"};
        output.attr("names") = varnames;
      }
    }
    
    CharacterVector newclasses = {"data.frame", "stageframe"};
    output.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, num_stages);
    output.attr("class") = newclasses;
    
    return output;
  }
  
  //' Base Skeleton Data Frame for Paramnames Objects
  //' 
  //' Creates a simple skeleton \code{paramnames} object that can be entered as
  //' input in functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
  //' \code{\link{aflefko2}()}.
  //' 
  //' @name paramnames_skeleton
  //' 
  //' @param name_terms A logical value indicating whether to start each variable
  //' name as \code{none} if \code{FALSE}, or as the default \code{modelparams}
  //' name if \code{TRUE}. Defaults to \code{FALSE}.
  //' 
  //' @return A three column data frame, of which the first describes the
  //' parameters in reasonably plain English, the second gives the name of the
  //' parameter within the MPM generating functions, and the third is to be
  //' edited with the names of the variables as they appear in the models.
  //' 
  //' @section Notes:
  //' The third column in the resulting object should be edited with the names only
  //' of those variables actually used in vital rate modeling. This
  //' \code{paramnames} object should apply to all models used in a single MPM
  //' building exercise. So, for example, if the models used include random terms,
  //' then they should all have the same random terms. Fixed terms can vary,
  //' however.
  //' 
  //' @keywords internal
  //' @noRd
  inline DataFrame paramnames_skeleton(bool name_terms = false) {
    
    CharacterVector parameter_names = {"time t", "individual", "patch",
      "alive in time t+1", "observed in time t+1", "sizea in time t+1",
      "sizeb in time t+1", "sizec in time t+1", "reproductive status in time t+1",
      "fecundity in time t+1", "fecundity in time t", "sizea in time t",
      "sizea in time t-1", "sizeb in time t", "sizeb in time t-1",
      "sizec in time t", "sizec in time t-1", "reproductive status in time t",
      "reproductive status in time t-1", "maturity status in time t+1",
      "maturity status in time t", "age in time t", "density in time t",
      "individual covariate a in time t", "individual covariate a in time t-1",
      "individual covariate b in time t", "individual covariate b in time t-1",
      "individual covariate c in time t", "individual covariate c in time t-1",
      "stage group in time t", "stage group in time t-1"};
    
    CharacterVector mainparams = {"year2", "individ", "patch", "surv3", "obs3",
      "size3", "sizeb3", "sizec3", "repst3", "fec3", "fec2", "size2", "size1",
      "sizeb2", "sizeb1", "sizec2", "sizec1", "repst2", "repst1", "matst3",
      "matst2", "age", "density", "indcova2", "indcova1", "indcovb2", "indcovb1",
      "indcovc2", "indcovc1", "group2", "group1"};
    
    CharacterVector modelparams = {"none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none", "none", "none", "none", "none",
      "none", "none", "none", "none", "none", "none", "none", "none"};
    
    if (name_terms) modelparams = mainparams;
    
    DataFrame output = DataFrame::create(_["parameter_names"] = parameter_names,
      _["mainparams"] = mainparams, _["modelparams"] = modelparams);
      
    return output;
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
  //' @param mat_input A logical value indicating whether the input matrix class
  //' is a standard \code{NumericMatrix}.
  //' @param sparse_switch A number denoting whether to force simple matrices to
  //' be projected as sparse (\code{1}) or not (\code{0}).
  //' 
  //' @return A list using the structure of a lefkoMat object.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List turbogeodiesel(DataFrame& loy, List Umats, List Fmats,
    DataFrame hstages, DataFrame agestages, DataFrame stages, bool patchmats,
    bool popmats, bool mat_input, int sparse_switch) {
    
    StringVector pops = as<StringVector>(loy["pop"]);
    arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
    StringVector patches = as<StringVector>(loy["patch"]);
    arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
    arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
    arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
    arma::uvec uniquepops = unique(pop_num);
    arma::uvec uniquepoppatches = unique(poppatchc);
    int loydim = pops.length();
    int numofpops = uniquepops.n_elem;
    int numofpatches = uniquepoppatches.n_elem;
    
    if (numofpatches == 1) popmats = false;
    
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
    
    // Needed matrix number and order determination
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
    
    if (patchmats && !popmats) {
      totalmatrices = toestimate.n_elem;
    } else if (!patchmats && popmats) {
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
    
    // Predicts which elements will be targeted for arithmetic mean estimation
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
            
          } else if (static_cast<int>(hstage2nin(i2)) == numstages || 
              static_cast<int>(hstage3in(i1)) == numstages) {
            hsindexl(counter) = (i1 * numhstages) + i2;
            counter++;
          }
        }
      }
    }
    
    arma::uvec hsgood = find(hsindexl);
    arma::uvec hsindex = hsindexl.elem(hsgood);
    arma::uvec zerovec(1, fill::zeros);
    arma::uvec allindices = join_cols(zerovec, hsindex);
    
    // Build U & F matrices of element-wise arithmetic means
    // Each column holds predicted non-zero elements of each mean matrix
    // Each matrix is a column vector within the overall matrix
    int core_elem = counter;
    
    arma::mat umatvec(core_elem, totalmatrices, fill::zeros);
    arma::mat fmatvec(core_elem, totalmatrices, fill::zeros);
    
    int patchchoice {0};
    int popchoice {0};
    
    pop_num = pop_num - 1;
    poppatchc = poppatchc - 1;
    
    for (int i = 0; i < loydim; i ++) {
      if (patchmats) {
        patchchoice = poppatchc(i);
        
        if (mat_input) {
          arma::mat Umats_invaded = as<arma::mat>(Umats(i));
          arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
          
          umatvec.col(patchchoice) = umatvec.col(patchchoice) +
            (Umats_invaded.elem(allindices) / yearsinpatch(i));
          fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
            (Umats_invaded.elem(allindices) / yearsinpatch(i));
        } else {
          arma::sp_mat Umats_invaded = as<arma::sp_mat>(Umats(i));
          arma::sp_mat Fmats_invaded = as<arma::sp_mat>(Fmats(i));
          
          arma::vec Umats_invaded_allindices (allindices.n_elem, fill::zeros);
          arma::vec Fmats_invaded_allindices (allindices.n_elem, fill::zeros);
          
          for (int j = 0; j < allindices.n_elem; j++) {
            Umats_invaded_allindices(j) = Umats_invaded(allindices(j));
            Fmats_invaded_allindices(j) = Fmats_invaded(allindices(j));
          }
        
          umatvec.col(patchchoice) = umatvec.col(patchchoice) +
            (Umats_invaded_allindices / yearsinpatch(i));
          fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
            (Umats_invaded_allindices / yearsinpatch(i));
        }
        
        if (popmats) {
          if (patchmats) {
            popchoice = numofpatches + pop_num(i);
          } else {
            popchoice = pop_num(i);
          }
          
          if (mat_input) {
            arma::mat Umats_invaded = as<arma::mat>(Umats(i));
            arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
            
            umatvec.col(popchoice) = umatvec.col(popchoice) +
              (Umats_invaded.elem(allindices) / (yearsinpatch(i) * patchesinpop(i)));
            fmatvec.col(popchoice) = fmatvec.col(popchoice) +
              (Umats_invaded.elem(allindices) / (yearsinpatch(i) * patchesinpop(i)));
          } else {
            arma::sp_mat Umats_invaded = as<arma::sp_mat>(Umats(i));
            arma::sp_mat Fmats_invaded = as<arma::sp_mat>(Fmats(i));
            
            arma::vec Umats_invaded_allindices (allindices.n_elem, fill::zeros);
            arma::vec Fmats_invaded_allindices (allindices.n_elem, fill::zeros);
            
            for (int j = 0; j < allindices.n_elem; j++) {
              Umats_invaded_allindices(j) = Umats_invaded(allindices(j));
              Fmats_invaded_allindices(j) = Fmats_invaded(allindices(j));
            }
          
            umatvec.col(popchoice) = umatvec.col(popchoice) +
              (Umats_invaded_allindices / (yearsinpatch(i) * patchesinpop(i)));
            fmatvec.col(popchoice) = fmatvec.col(popchoice) +
              (Umats_invaded_allindices / (yearsinpatch(i) * patchesinpop(i)));
          }
        }
      }
    }
    
    // New labels data frame
    int cheatsheetlength {1};
    
    if (patchmats && !popmats) {
      cheatsheetlength = numofpatches;
      
    } else if (!patchmats && popmats) {
      cheatsheetlength = numofpops;
      
    } else if (popmats && patchmats) {
      cheatsheetlength = numofpops + numofpatches;
    }
    
    StringVector poporder_redone(cheatsheetlength);
    StringVector patchorder_redone(cheatsheetlength);
    
    int new_labels_counter {0};
    
    if (patchmats) {
      for (int i = 0; i < numofpatches; i++) {
        poporder_redone(i) = poporderlong(i);
        patchorder_redone(i) = patchorderlong(i);
        
        new_labels_counter++;
      }
    }
    
    if (popmats) {
      for (int i = 0; i < numofpops; i++) {
        poporder_redone(new_labels_counter) = uniquepops_str(i);
        patchorder_redone(new_labels_counter) = "0";
        
        new_labels_counter++;
      }
    }
    
    DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone,
      _["patch"] = patchorder_redone);
    
    // Create the main list objects holding matrices
    List U (totalmatrices);
    List F (totalmatrices);
    List A (totalmatrices);
    int totalutrans {0};
    int totalftrans {0};
    
    if (mat_input && sparse_switch == 0) {
      for (int i = 0; i < totalmatrices; i++) {
        arma::mat umat_base(numhstages, numhstages, fill::zeros);
        arma::mat fmat_base(numhstages, numhstages, fill::zeros);
        
        umat_base.elem(allindices) = umatvec.col(i);
        fmat_base.elem(allindices) = fmatvec.col(i);
        arma::mat amat_base = umat_base + fmat_base;
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = utrans.n_elem;
      totalftrans = ftrans.n_elem;
      
    } else {
      for (int i = 0; i < totalmatrices; i++) {
        arma::sp_mat umat_base(numhstages, numhstages);
        arma::sp_mat fmat_base(numhstages, numhstages);
        
        for (int j = 0; j < allindices.n_elem; j++) {
          umat_base(allindices(j)) = umatvec(j, i);
          fmat_base(allindices(j)) = fmatvec(j, i);
        }
        arma::sp_mat amat_base = umat_base + fmat_base;
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = utrans.n_elem;
      totalftrans = ftrans.n_elem;
    }
    
    NumericVector matrixqc(3);
    matrixqc(0) = totalutrans; // summed number of non-zero u transitions
    matrixqc(1) = totalftrans; // summed number of non-zero f transitions
    matrixqc(2) = totalmatrices;
    
    // Final output
    List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F,
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
  //' @param mat_input A logical value indicating whether the input matrix class
  //' is a standard \code{NumericMatrix}.
  //' @param sparse_switch A number denoting whether to force simple matrices to
  //' be projected as sparse (\code{1}) or not (\code{0}).
  //' 
  //' @return A list using the structure of a LefkoMat object.
  //' 
  //' @keywords internal
  //' @noRd
  inline Rcpp::List geodiesel(DataFrame& loy, List Umats, List Fmats,
    DataFrame agestages, DataFrame stages, bool patchmats, bool popmats,
    bool mat_input, int sparse_switch) {
    
    StringVector pops = as<StringVector>(loy["pop"]);
    arma::uvec pop_num = as<arma::uvec>(loy["popc"]);
    StringVector patches = as<StringVector>(loy["patch"]);
    arma::uvec poppatchc = as<arma::uvec>(loy["poppatchc"]);
    arma::uvec patchesinpop = as<arma::uvec>(loy["patchesinpop"]);
    arma::uvec yearsinpatch = as<arma::uvec>(loy["yearsinpatch"]);
    arma::uvec uniquepops = unique(pop_num);
    arma::uvec uniquepoppatches = unique(poppatchc);
    int loydim = pops.length();
    int numofpops = uniquepops.n_elem;
    int numofpatches = uniquepoppatches.n_elem;
    
    if (numofpatches == 1) popmats = false;
    
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
    
    // Needed matrix number and order determination
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
    
    if (patchmats && !popmats) {
      totalmatrices = toestimate.n_elem;
    } else if (!patchmats && popmats) {
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
    
    // Predicts which elements will be targeted for arithmetic mean estimation
    arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
    int initialstages = astages.n_elem;
    
    // Test for the presence of ages, and determine matrix dimensions
    int colsused {0};
    
    if (mat_input && sparse_switch == 0) {
      arma::mat initUmat = as<arma::mat>(Umats(0));
      colsused = initUmat.n_cols;
    } else if (mat_input) {
      arma::sp_mat initUmat(as<arma::mat>(Umats(0)));
      colsused = initUmat.n_cols;
    } else {
      arma::sp_mat initUmat = as<arma::sp_mat>(Umats(0));
      colsused = initUmat.n_cols;
    }
    
    int agemultiplier = colsused / initialstages;
    
    int numstages = astages.n_elem * agemultiplier;
    
    // Build U & F matrices of element-wise arithmetic means
    // Each column holds predicted non-zero elements of each mean matrix
    // Each matrix is a column vector within the overall matrix
    int core_elem = numstages * numstages;
    
    arma::mat umatvec;
    arma::mat fmatvec;
    arma::mat amatvec;
    arma::sp_mat umatvec_sp;
    arma::sp_mat fmatvec_sp;
    arma::sp_mat amatvec_sp;
    
    if (mat_input && sparse_switch == 0) {
      arma::mat umatvec_(core_elem, totalmatrices);
      arma::mat fmatvec_(core_elem, totalmatrices);
      umatvec_.zeros();
      fmatvec_.zeros();
      
      int patchchoice {0};
      int popchoice {0};
      
      pop_num = pop_num - 1;
      poppatchc = poppatchc - 1;
      
      for (int i = 0; i < loydim; i ++) {
        if (patchmats) {
          patchchoice = poppatchc(i);
          
          arma::mat Umats_invaded = as<arma::mat>(Umats(i));
          arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
          
          umatvec_.col(patchchoice) = umatvec_.col(patchchoice) +
            (arma::vectorise(Umats_invaded) / yearsinpatch(i));
          fmatvec_.col(patchchoice) = fmatvec_.col(patchchoice) +
            (arma::vectorise(Fmats_invaded) / yearsinpatch(i));
        }
        
        if (popmats) {
          if (patchmats) {
            popchoice = numofpatches + pop_num(i);
          } else {
            popchoice = pop_num(i);
          }
          
          arma::mat Umats_invaded = as<arma::mat>(Umats(i));
          arma::mat Fmats_invaded = as<arma::mat>(Fmats(i));
          
          umatvec_.col(popchoice) = umatvec_.col(popchoice) +
            (arma::vectorise(Umats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
          fmatvec_.col(popchoice) = fmatvec_.col(popchoice) +
            (arma::vectorise(Fmats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
        }
      }
      umatvec = umatvec_;
      fmatvec = fmatvec_;
      amatvec = umatvec + fmatvec;
      
    } else {
      arma::sp_mat umatvec_sp_(core_elem, totalmatrices);
      arma::sp_mat fmatvec_sp_(core_elem, totalmatrices);
      
      int patchchoice {0};
      int popchoice {0};
      
      pop_num = pop_num - 1;
      poppatchc = poppatchc - 1;
      
      for (int i = 0; i < loydim; i ++) {
        if (patchmats) {
          patchchoice = poppatchc(i);
          
          arma::sp_mat Umats_invaded;
          arma::sp_mat Fmats_invaded;
          
          if (mat_input) { 
            arma::sp_mat Umats_invaded_(as<arma::mat>(Umats(i)));
            arma::sp_mat Fmats_invaded_(as<arma::mat>(Fmats(i)));
            Umats_invaded = Umats_invaded_;
            Fmats_invaded = Fmats_invaded_;
          } else { 
            Umats_invaded = as<arma::sp_mat>(Umats(i));
            Fmats_invaded = as<arma::sp_mat>(Fmats(i));
          }
          
          umatvec_sp_.col(patchchoice) = umatvec_sp_.col(patchchoice) +
            (arma::vectorise(Umats_invaded) / yearsinpatch(i));
          fmatvec_sp_.col(patchchoice) = fmatvec_sp_.col(patchchoice) +
            (arma::vectorise(Fmats_invaded) / yearsinpatch(i));
        }
        
        if (popmats) {
          if (patchmats) {
            popchoice = numofpatches + pop_num(i);
          } else {
            popchoice = pop_num(i);
          }
          
          arma::sp_mat Umats_invaded;
          arma::sp_mat Fmats_invaded;
          
          if (mat_input) {
            arma::sp_mat Umats_invaded_(as<arma::mat>(Umats(i)));
            arma::sp_mat Fmats_invaded_(as<arma::mat>(Fmats(i)));
            Umats_invaded = Umats_invaded_;
            Fmats_invaded = Fmats_invaded_;
          } else {
            Umats_invaded = as<arma::sp_mat>(Umats(i));
            Fmats_invaded = as<arma::sp_mat>(Fmats(i));
          }
          
          umatvec_sp_.col(popchoice) = umatvec_sp_.col(popchoice) +
            (arma::vectorise(Umats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
          fmatvec_sp_.col(popchoice) = fmatvec_sp_.col(popchoice) +
            (arma::vectorise(Fmats_invaded) / (yearsinpatch(i) * patchesinpop(i)));
        }
      }
      
      umatvec_sp = umatvec_sp_;
      fmatvec_sp = fmatvec_sp_;
      amatvec_sp = umatvec_sp + fmatvec_sp;
    }
    
    // New labels data frame
    int cheatsheetlength {1};
    
    if (patchmats && !popmats) {
      cheatsheetlength = numofpatches;
      
    } else if (!patchmats && popmats) {
      cheatsheetlength = numofpops;
      
    } else if (popmats && patchmats) {
      cheatsheetlength = numofpops + numofpatches;
    }
    StringVector poporder_redone(cheatsheetlength);
    StringVector patchorder_redone(cheatsheetlength);
    
    int new_labels_counter {0};
    
    if (patchmats) {
      for (int i = 0; i < numofpatches; i++) {
        poporder_redone(i) = poporderlong(i);
        patchorder_redone(i) = patchorderlong(i);
        
        new_labels_counter++;
      }
    }
    
    if (popmats) {
      for (int i = 0; i < numofpops; i++) {
        poporder_redone(new_labels_counter) = uniquepops_str(i);
        patchorder_redone(new_labels_counter) = "0";
        
        new_labels_counter++;
      }
    }
    
    DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, 
      _["patch"] = patchorder_redone);
    
    // Create the main list objects to hold the matrices
    List U (totalmatrices);
    List F (totalmatrices);
    List A (totalmatrices);
    int totalutrans {0};
    int totalftrans {0};
    
    if (mat_input && sparse_switch == 0) {
      for (int i = 0; i < totalmatrices; i++) {
        arma::mat umat_base = umatvec.col(i);
        arma::mat fmat_base = fmatvec.col(i);
        arma::mat amat_base = amatvec.col(i);
        
        umat_base.reshape(numstages, numstages);
        fmat_base.reshape(numstages, numstages);
        amat_base.reshape(numstages, numstages);
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = utrans.n_elem;
      totalftrans = ftrans.n_elem;
      
    } else {
      for (int i = 0; i < totalmatrices; i++) {
        arma::sp_mat umat_base = umatvec_sp.col(i);
        arma::sp_mat fmat_base = fmatvec_sp.col(i);
        arma::sp_mat amat_base = amatvec_sp.col(i);
        
        umat_base.reshape(numstages, numstages);
        fmat_base.reshape(numstages, numstages);
        amat_base.reshape(numstages, numstages);
        
        U(i) = umat_base;
        F(i) = fmat_base;
        A(i) = amat_base;
      }
      
      // Matrix QC output
      arma::uvec utrans = find(umatvec);
      arma::uvec ftrans = find(fmatvec);
      totalutrans = utrans.n_elem;
      totalftrans = ftrans.n_elem;
    }
    
    NumericVector matrixqc(3);
    matrixqc(0) = totalutrans; // summed number of U transitions
    matrixqc(1) = totalftrans; // summed number of F transitions
    matrixqc(2) = totalmatrices;
    
    // Final output
    List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F,
      _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
    
    return output;
  }
}
#endif