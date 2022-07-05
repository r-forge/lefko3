#ifndef LEFKOUTILS_H
#define LEFKOUTILS_H

using namespace Rcpp;
using namespace arma;

//' Append NumericVector to the End of Another NumericVector
//' 
//' This function appends one NumericVector fully to another.
//' 
//' @name concat_dbl
//' 
//' @param A Any NumericVector.
//' @param B Any other NumericVector.
//' 
//' @return Returns a new NumericVector with elements of vector A followed by
//' elements of vector B.
//'
//' @keywords internal
//' @noRd
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
//' @name concat_int
//' 
//' @param A Any IntegerVector.
//' @param B Any other IntegerVector.
//' 
//' @return Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @keywords internal
//' @noRd
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
//' @name concat_str
//' 
//' @param A Any StringVector.
//' @param B Any other StringVector.
//' 
//' @return Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @keywords internal
//' @noRd
StringVector concat_str(StringVector x, StringVector y) {
  
  std::vector<std::string> xconv = as<std::vector<std::string> >(x);
  std::vector<std::string> yconv = as<std::vector<std::string> >(y);
  std::vector<std::string> z(x.size() + y.size());
  
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  
  StringVector zconv(z.begin(), z.end());
  
  return(zconv);
}

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
  
  if (str1_length >= str2_length && str2_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (!lower) {
        if (str1[i] != str2[rem_check]) {
          rem_check = 0;
        } else {
          rem_check += 1;
          if (rem_check >= str2_length) break;
        }
      } else {
        if (tolower(str1[i]) != tolower(str2[rem_check])) {
          rem_check = 0;
        } else {
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
//' variables in the zero-inflation model. Here, given as \code{NULL}.}
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
//' variables in the zero-inflation model. Not used in \code{lm}/\code{glm}/
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
//' variables in the zero-inflation model. Not used in \code{vglm} objects.}
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
//' variables in the zero-inflation model. Not used in \code{zeroinfl} objects.}
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
//' variables in the zero-inflation model. Not used in lme4.}
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
//' variables in the zero-inflation model.}
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
//' @param object An S3 vital rate model other than objects of class
//' \code{vrm_input}. Currently, this should be output from functions
//' \code{lm()}, \code{glm()}, \code{glm.nb()}, \code{zeroinfl()}, and
//' \code{glmmTMB()}.
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
//' variables in the zero-inflation model.}
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
//' variables in the zero-inflation model.}
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
  
  List output;
  
  if (stringcompare_hard(model_class, "vglm")) {
    output = vglm_extractor(object);
    
  } else if (stringcompare_hard(model_class, "lmerMod") || 
      stringcompare_hard(model_class, "glmerMod")) {
    output = lme4_extractor(object);
    
  } else {
    throw Rcpp::exception("Model type unrecognized.", false);
  }
  
  return output;
}

//' Extract Key Components of vrm_input Objects
//' 
//' This function extracts the components of a \code{vrm_input} object for
//' functions \code{jerzeibalowski()} and \code{motherbalowski()}.
//' 
//' @name vrm_extractor
//' 
//' @param object A \code{vrm_input} object.
//' 
//' @return A list with the following elements:
//' \item{class}{The exact class of \code{object}. Will be \code{vrm_input}.}
//' \item{family}{The response distribution.}
//' \item{dist}{An integer representing the response distribution. \code{0} = 
//' poisson, \code{1} = negbin, \code{2} = gaussian, \code{3} = gamma, \code{4}
//' = binomial, and \code{5} = constant.}
//' \item{zero_inflated}{A logical value indicating whether the distribution is
//' zero-inflated. Not used in \code{lm}/\code{glm}/\code{negbin} objects. Will
//' equal \code{FALSE} if \code{dist = 5}.}
//' \item{zero_truncated}{A logical value indicating whether the distribution is
//' zero-truncated. Not used in \code{lm}/\code{glm}/\code{negbin} objects. Will
//' equal \code{FALSE} if \code{dist = 5}.}
//' \item{all_vars}{A vector holding the names of each variable used by
//' \code{object}. Will equal only the intercept if \code{dist = 5}.}
//' \item{fixed_vars}{A string vector holding the names of the fixed variables.
//' Will equal only the intercept if \code{dist = 5}.}
//' \item{fixed_slopes}{A numeric vector holding the slope coefficients of the
//' fixed variables, in the same order as \code{fixed_vars}. Will equal only the
//' intercept if \code{dist = 5}.}
//' \item{fixed_zi_vars}{A string vector holding the names of the zero-inflated
//' fixed variables. Not used in \code{lm}/\code{glm}/\code{negbin} objects.
//' Will equal \code{NULL} is \code{dist = 5}.}
//' \item{fixed_zi_slopes}{A numeric vector holding the slope coefficients of
//' the zero-inflated fixed variables, in the same order as
//' \code{fixed_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
//' objects. Will equal \code{NULL} is \code{dist = 5}.}
//' \item{random_zi_vars}{A string vector holding the names of the random
//' variables in the zero-inflation model. Not used in \code{lm}/\code{glm}/
//' \code{negbin} objects. Will equal \code{NULL} is \code{dist = 5}.}
//' \item{random_zi_slopes}{A numeric vector holding the slope coefficients of
//' the random variables in the zero-inflation model, in the same order as
//' \code{random_zi_vars}. Not used in \code{lm}/\code{glm}/\code{negbin}
//' objects. Will equal \code{NULL} is \code{dist = 5}.}
//' \item{sigma}{The residual standard deviation of the model. Defaults to
//' \code{1.0}.}
//' \item{theta}{The estimated theta, if the response is negative binomial.
//' Otherwise, will equal \code{1.0}.}
//' 
//' @keywords internal
//' @noRd
List vrm_extractor(List object) {
  std::string object_class = "vrm_input";
  
  String distrib = as<String>(object["dist"]);
  double sigma_theta = object["sigma_theta"];
  
  String resp_family = "gaussian";
  int dist = 2;
  double theta {1.0};
  double sigma = sigma_theta;
  bool zero_inflation = false;
  bool zero_truncation = false;
  
  if (stringcompare_simple(distrib, "poisson")) {
    resp_family = "poisson";
    dist = 0;
  } else if (stringcompare_simple(distrib, "negbin")) {
    resp_family = "negbin";
    dist = 1;
    theta = sigma_theta;
  } else if (stringcompare_simple(distrib, "gamma")) {
    resp_family = "gamma";
    dist = 3;
  } else if (stringcompare_simple(distrib, "binom")) {
    resp_family = "binomial";
    dist = 4;
  } else if (stringcompare_simple(distrib, "cons")) {
    resp_family = "constant";
    dist = 5;
  }
  
  if (stringcompare_simple(distrib, "trunc")) {
    if (dist == 0 || dist == 1) zero_truncation = true; 
  }
  
  NumericVector coefs = as<NumericVector>(object["fixed_slopes"]);
  NumericVector zi_coefs = as<NumericVector>(object["fixed_zi"]);
  
  CharacterVector all_vars = as<CharacterVector>(object("effects_names"));
  
  NumericVector year_slopes = as<NumericVector>(object["year_slopes"]);
  NumericVector patch_slopes = as<NumericVector>(object["patch_slopes"]);
  NumericVector indcova2_slopes = as<NumericVector>(object["indcova2_slopes"]);
  NumericVector indcova1_slopes = as<NumericVector>(object["indcova1_slopes"]);
  NumericVector indcovb2_slopes = as<NumericVector>(object["indcovb2_slopes"]);
  NumericVector indcovb1_slopes = as<NumericVector>(object["indcovb1_slopes"]);
  NumericVector indcovc2_slopes = as<NumericVector>(object["indcovc2_slopes"]);
  NumericVector indcovc1_slopes = as<NumericVector>(object["indcovc1_slopes"]);
  
  NumericVector year_zi = as<NumericVector>(object["year_zi"]);
  NumericVector patch_zi = as<NumericVector>(object["patch_zi"]);
  NumericVector indcova2_zi = as<NumericVector>(object["indcova2_zi"]);
  NumericVector indcova1_zi = as<NumericVector>(object["indcova1_zi"]);
  NumericVector indcovb2_zi = as<NumericVector>(object["indcovb2_zi"]);
  NumericVector indcovb1_zi = as<NumericVector>(object["indcovb1_zi"]);
  NumericVector indcovc2_zi = as<NumericVector>(object["indcovc2_zi"]);
  NumericVector indcovc1_zi = as<NumericVector>(object["indcovc1_zi"]);
  
  double zi_year_sum = sum(year_zi);
  double zi_patch_sum = sum(patch_zi);
  double zi_indcova2_sum = sum(indcova2_zi);
  double zi_indcova1_sum = sum(indcova1_zi);
  double zi_indcovb2_sum = sum(indcovb2_zi);
  double zi_indcovb1_sum = sum(indcovb1_zi);
  double zi_indcovc2_sum = sum(indcovc2_zi);
  double zi_indcovc1_sum = sum(indcovc1_zi);
  
  double zi_sum = zi_year_sum + zi_patch_sum + zi_indcova2_sum + zi_indcova1_sum +
    zi_indcovb2_sum + zi_indcovb1_sum + zi_indcovc2_sum + zi_indcovc1_sum;
  
  if (zi_sum > 0.0) {
    zero_inflation = true;
  }
  
  if (zero_truncation && zero_inflation) {
    throw Rcpp::exception("Models cannot be both zero-inflated and zero-truncated.",
      false);
  }
  
  CharacterVector year_names = as<CharacterVector>(object["year_names"]);
  CharacterVector patch_names = as<CharacterVector>(object["patch_names"]);
  CharacterVector indcova_names = as<CharacterVector>(object["indcova_names"]);
  CharacterVector indcovb_names = as<CharacterVector>(object["indcovb_names"]);
  CharacterVector indcovc_names = as<CharacterVector>(object["indcovc_names"]);
  
  year_slopes.attr("names") = year_names;
  patch_slopes.attr("names") = patch_names;
  indcova2_slopes.attr("names") = indcova_names;
  indcova1_slopes.attr("names") = indcova_names;
  indcovb2_slopes.attr("names") = indcovb_names;
  indcovb1_slopes.attr("names") = indcovb_names;
  indcovc2_slopes.attr("names") = indcovc_names;
  indcovc1_slopes.attr("names") = indcovc_names;
  
  CharacterVector randomvar_names = {"year2", "patch", "indcova2", "indcova1",
    "indcovb2", "indcovb1", "indcovc2", "indcovc1"};
  
  List random_slopes (8);
  random_slopes(0) = year_slopes;
  random_slopes(1) = patch_slopes;
  random_slopes(2) = indcova2_slopes;
  random_slopes(3) = indcova1_slopes;
  random_slopes(4) = indcovb2_slopes;
  random_slopes(5) = indcovb1_slopes;
  random_slopes(6) = indcovc2_slopes;
  random_slopes(7) = indcovc1_slopes;
  
  List random_zi (8);
  random_zi(0) = year_zi;
  random_zi(1) = patch_zi;
  random_zi(2) = indcova2_zi;
  random_zi(3) = indcova1_zi;
  random_zi(4) = indcovb2_zi;
  random_zi(5) = indcovb1_zi;
  random_zi(6) = indcovc2_zi;
  random_zi(7) = indcovc1_zi;
  
  List random_names (8);
  random_names(0) = year_names;
  random_names(1) = patch_names;
  random_names(2) = indcova_names;
  random_names(3) = indcova_names;
  random_names(4) = indcovb_names;
  random_names(5) = indcovb_names;
  random_names(6) = indcovc_names;
  random_names(7) = indcovc_names;
  
  random_slopes.attr("names") = randomvar_names;
  random_zi.attr("names") = randomvar_names;
  random_names.attr("names") = randomvar_names;
  
  NumericVector group2_slopes = as<NumericVector>(object["group2_slopes"]);
  NumericVector group1_slopes = as<NumericVector>(object["group1_slopes"]);
  NumericVector group2_zi = as<NumericVector>(object["group2_zi"]);
  NumericVector group1_zi = as<NumericVector>(object["group1_zi"]);
  CharacterVector group2_names = clone(as<CharacterVector>(object["group_names"]));
  CharacterVector group1_names = clone(as<CharacterVector>(object["group_names"]));
  CharacterVector group2_zi_names = clone(as<CharacterVector>(object["group_names"]));
  CharacterVector group1_zi_names = clone(as<CharacterVector>(object["group_names"]));
  
  CharacterVector zi_vars = clone(all_vars);
  
  int no_group2_slopes = group2_slopes.length();
  int no_group1_slopes = group1_slopes.length();
  
  if (no_group2_slopes > 0) {
    for (int i = 0; i < no_group2_slopes; i++) {
      group2_names[i] = "group2" + group2_names[i];
    }
    all_vars = concat_str(all_vars, group2_names);
    coefs = concat_dbl(coefs, group2_slopes);
  }
  
  if (no_group1_slopes > 0) {
    for (int i = 0; i < no_group1_slopes; i++) {
      group1_names[i] = "group1" + group1_names[i];
    }
    all_vars = concat_str(all_vars, group1_names);
    coefs = concat_dbl(coefs, group1_slopes);
  }
  
  int no_group2_zi = group2_zi.length();
  int no_group1_zi = group1_zi.length();
  
  if (no_group2_zi > 0) {
    for (int i = 0; i < no_group2_zi; i++) {
      group1_zi_names[i] = "group2" + group2_zi_names[i];
    }
    zi_vars = concat_str(zi_vars, group2_zi_names);
    zi_coefs = concat_dbl(zi_coefs, group2_zi);
  }
  
  if (no_group1_zi > 0) {
    for (int i = 0; i < no_group1_zi; i++) {
      group1_zi_names[i] = "group1" + group1_zi_names[i];
    }
    zi_vars = concat_str(zi_vars, group1_zi_names);
    zi_coefs = concat_dbl(zi_coefs, group1_zi);
  }
  
  Rcpp::List output;
  if (dist == 5) {
    CharacterVector simple_vars = {"Intercept"};
    NumericVector simple_coefs (1);
    simple_coefs(0) = coefs(0);
    
    output = List::create(_["class"] = object_class, _["family"] = resp_family,
      _["dist"] = dist, _["zero_inflated"] = false, _["zero_truncated"] = false,
      _["all_vars"] = simple_vars, _["fixed_vars"] = simple_vars,
      _["fixed_slopes"] = simple_coefs, _["fixed_zi_vars"] = R_NilValue,
      _["fixed_zi_slopes"] = R_NilValue, _["random_vars"] = R_NilValue,
      _["random_slopes"] = R_NilValue, _["random_zi_vars"] = R_NilValue,
      _["random_zi_slopes"] = R_NilValue, _["sigma"] = sigma, _["theta"] = theta);
    
  } else  {
    output = List::create(_["class"] = object_class, _["family"] = resp_family,
      _["dist"] = dist, _["zero_inflated"] = zero_inflation,
      _["zero_truncated"] = zero_truncation, _["all_vars"] = all_vars,
      _["fixed_vars"] = all_vars, _["fixed_slopes"] = coefs,
      _["fixed_zi_vars"] = zi_vars, _["fixed_zi_slopes"] = zi_coefs,
      _["random_vars"] = random_names, _["random_slopes"] = random_slopes,
      _["random_zi_vars"] = random_names, _["random_zi_slopes"] = random_zi,
      _["sigma"] = sigma, _["theta"] = theta);
  }
  
  return output;
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

//' Extract Coefficients From Linear Vital Rate Models
//' 
//' Function \code{modelextract()} extracts coefficient values from linear
//' models estimated through various linear modeling functions in R, to estimate
//' vital rates in \code{lefko3}. Used to supply coefficients to
//' \code{\link{flefko3}()}, \code{\link{flefko2}()}, \code{\link{fleslie()}},
//' and \code{\link{aflefko2}()}.
//' 
//' @name modelextract
//' 
//' @param object A linear model estimated through one of the methods used in
//' function \code{\link{modelsearch}()}, or a \code{vrm_input} object.
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
//' @param nodata A logical value used to determine whether to use
//' \code{vrm_input} methods. Defaults to \code{FALSE}, in which case models
//' were developed with function \code{modelsearch()}.
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
List modelextract(RObject object, DataFrame paramnames, NumericVector mainyears,
  CharacterVector mainpatches, RObject maingroups, RObject mainindcova,
  RObject mainindcovb, RObject mainindcovc, bool nodata = false) {
  
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
    
  } else if (is<List>(object) && !nodata) {
    core_components = S3_extractor(as<List>(object));
    
  } else if (nodata) {
    core_components = vrm_extractor(as<List>(object));
    
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
  
  int no_fixed_zi_slopes = fixed_zi_slopes.length();
  
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
  
  // Individual covariates
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
  
  // Rcout << "no_fixed_vars: " << no_fixed_vars << "\n";
  // Rcout << "no_fixed_zi_vars: " << no_fixed_zi_vars << "\n";
  // Rcout << "no_fixed_zi_slopes: " << no_fixed_zi_slopes << "\n";
  
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
  
  for (int i = 0; i < no_fixed_zi_slopes; i++) {
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
    
    if (no_fixed_zi_slopes > 0) {
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
    int no_random_vars = random_slopes.length();
    
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
        // Rcout << "\n random indcova2 sorter reached \n";
        CharacterVector ran_inda2_names = random_vars[i];
        NumericVector ran_inda2_slopes = random_slopes[i];
        int no_ran_inda2_slopes = ran_inda2_names.length();
      
        for (int j = 0; j < no_ran_inda2_slopes; j++) {
          for (int k = 0; k < no_indcova_names; k++) {
            // Rcout << "ran_inda2_names(j) = " << ran_inda2_names(j) << "\n";
            // Rcout << "indcova_names(k) = " << indcova_names(k) << "\n";
            if (stringcompare_hard(as<std::string>(ran_inda2_names(j)), as<std::string>(indcova_names(k)))) {
              // Rcout << "\n random indcova2 sorting... \n";
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
      int no_random_zi_vars = random_zi_slopes.length();
      
      for (int i = 0; i < no_random_zi_vars; i++) {
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), year2var)) {
          CharacterVector ran_zi_year_names = random_zi_vars[i];
          NumericVector ran_zi_year_slopes = random_zi_slopes[i];
          int no_ran_zi_year_slopes = ran_zi_year_slopes.length();
          
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
          int no_ran_zi_patch_slopes = ran_zi_patch_slopes.length();
          
          for (int j = 0; j < no_ran_zi_patch_slopes; j++) {
            for (int k = 0; k < no_patches; k++) {
              if (stringcompare_hard(as<std::string>(ran_zi_patch_names(j)),
                  as<std::string>(mainpatches(k)))) {
                zeropatch_coefs(k) = ran_zi_patch_slopes(j);
              }
            }
          }
        }
        
        if (stringcompare_hard(as<std::string>(random_zi_names(i)), indcova2var)) {
          CharacterVector ran_zi_inda2_names = random_zi_vars[i];
          NumericVector ran_zi_inda2_slopes = random_zi_slopes[i];
          int no_ran_zi_inda2_slopes = ran_zi_inda2_slopes.length();
        
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
          int no_ran_zi_inda1_slopes = ran_zi_inda1_slopes.length();
        
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
          int no_ran_zi_indb2_slopes = ran_zi_indb2_slopes.length();
        
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
          int no_ran_zi_indb1_slopes = ran_zi_indb1_slopes.length();
        
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
          int no_ran_zi_indc2_slopes = ran_zi_indc2_slopes.length();
        
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
          int no_ran_zi_indc1_slopes = ran_zi_indc1_slopes.length();
        
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
List decomp3sp(arma::mat Amat) {
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

#endif
