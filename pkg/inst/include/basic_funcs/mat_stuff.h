#ifndef LEFKOUTILS_mat_stuff_H
#define LEFKOUTILS_mat_stuff_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

namespace LefkoMats {
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
  inline arma::vec flagrantcrap(const arma::mat& Xmat,
    const arma::uvec& allindices) {
    
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
  inline arma::vec moreflagrantcrap(const arma::mat& Xmat) {
    
    arma::vec newcol = arma::vectorise(Xmat);
    
    return newcol;
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
  inline Rcpp::List turbogeodiesel(DataFrame& loy, List Umats, List Fmats,
    DataFrame hstages, DataFrame agestages, DataFrame stages, bool patchmats,
    bool popmats) {
    
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
    arma::uvec zerovec(1);
    zerovec.zeros();
    arma::uvec allindices = join_cols(zerovec, hsindex);
    
    // Build U and F matrices of element-wise arithmetic means
    // Each column corresponds to the predicted non-zero elements of each mean
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
      if (patchmats) {
        patchchoice = poppatchc(i);
        
        umatvec.col(patchchoice) = umatvec.col(patchchoice) +
          (flagrantcrap(as<arma::mat>(Umats[i]), allindices) / yearsinpatch(i));
        fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
          (flagrantcrap(as<arma::mat>(Fmats[i]), allindices) / yearsinpatch(i));
      }
      
      if (popmats) {
        if (patchmats) {
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
  inline Rcpp::List geodiesel(DataFrame& loy, List Umats, List Fmats,
    DataFrame agestages, DataFrame stages, bool patchmats, bool popmats) {
    
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
    arma::mat initUmat = Umats(0);
    int colsused = initUmat.n_cols;
    int agemultiplier = colsused / initialstages;
    
    int numstages = astages.n_elem * agemultiplier;
    
    // Build U and F matrices of element-wise arithmetic means
    // Each column corresponds to the predicted non-zero elements of each mean
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
      if (patchmats) {
        patchchoice = poppatchc(i);
        
        umatvec.col(patchchoice) = umatvec.col(patchchoice) +
          (moreflagrantcrap(as<arma::mat>(Umats[i])) / yearsinpatch(i));
        fmatvec.col(patchchoice) = fmatvec.col(patchchoice) +
          (moreflagrantcrap(as<arma::mat>(Fmats[i])) / yearsinpatch(i));
      }
      
      if (popmats) {
        if (patchmats) {
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
      _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
    
    return output;
  }
  
}

#endif