// [[Rcpp::depends(RcppArmadillo)]]
#define BOOST_DISABLE_ASSERTS

#include <RcppArmadillo.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using namespace arma;

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
//' 
//' @return A logical value indicating whether \code{str2} occurs within
//' \code{str1}.
//' 
//' @keywords internal
//' @noRd
bool stringcompare_simple(std::string str1, std::string str2) {
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
  
  return same;
}

//' Compares Three Strings for Interaction Notation
//' 
//' This function compares checks to see if one string is composed of the other
//' two string in R's interaction notation.
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

//' Re-index Projection Matrix On Basis of Overwrite Table
//' 
//' Function \code{.ovreplace()} takes matrix indices provided by functions
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
                    } 
                  }
                } else repentry3(currentindex) = 0.0;
                
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
  // if (cont) totalages += 1;
  
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

//' Create Historically Structured Version of ahMPM
//' 
//' Function \code{.thefifthhousemate()} takes an ahistorical MPM as input, and
//' uses the \code{allstages} index to create a historically structured version
//' of it.
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
  Rcpp::IntegerVector new_index = allstages["index321"];
  
  int num_mats = old_Umats.length();
  int index_elems = new_index.length();
  
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
      new_U(new_index(j)) = old_U(old_index(j));
      new_F(new_index(j)) = old_F(old_index(j));
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
//' @param ppy A data frame showing the population, patch, and year of each
//' matrix to create, in order.
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
List motherbalowski(DataFrame ppy, IntegerVector actualages, DataFrame ageframe,
  List survproxy, List fecproxy, NumericVector f2_inda, NumericVector f1_inda,
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
  
  if (no_indcova_names == 0) {
    indcova2s = {0};
    indcova1s = {0};
    
    zeroindcova2s = {0};
    zeroindcova1s = {0};
  }
  if (no_indcovb_names == 0) {
    indcovb2s = {0};
    indcovb1s = {0};
    
    zeroindcovb2s = {0};
    zeroindcovb1s = {0};
  }
  if (no_indcovc_names == 0) {
    indcovc2s = {0};
    indcovc1s = {0};
    
    zeroindcovc2s = {0};
    zeroindcovc1s = {0};
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
//' This function takes the various vital rate models and other parameters and
//' coordinates them as input into the function-based matrix estimation
//' functions.
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
    
    List madsexmadrigal_oneyear = motherbalowski(listofyears, actualages, ageframe,
      surv_proxy, fec_proxy, f2_inda, f1_inda, f2_indb, f1_indb, f2_indc, f1_indc,
      r2_inda, r1_inda, r2_indb, r1_indb, r2_indc, r1_indc, surv_dev, fec_dev, dens,
      fecmod, finalage, negfec, yearnumber, patchnumber, exp_tol, theta_tol, simplicity);
    
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

