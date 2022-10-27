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

}

#endif