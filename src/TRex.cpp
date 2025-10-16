// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]


#include <RcppArmadillo.h>
#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace arma;

//Claim used matrix
arma::vec xct(int N);
arma::vec nct(int N);

//Update_Expectations
arma::vec Update_mu_ij(arma::vec m_ij);
arma::vec Update_z_ijk_case1(arma::vec mu_ij);
arma::vec Update_z_ijk_case0(arma::vec mu_ij);
arma::vec Update_m_ij(arma::vec a_i, arma::vec b_i, arma::vec T_i,
                      arma::vec z_ijk_case1, arma::vec z_ijk_case0, arma::vec nct, arma::vec xct);
arma::vec Update_S_ij(arma::vec a_i, arma::vec b_i, arma::vec nct);
arma::vec Update_T_i(arma::vec Ji, arma::vec group_indicator, arma::vec a_i, arma::vec b_i, arma::vec m_ij,
                     double a_star, double b_star, double mu_star);
arma::vec Update_V_i(arma::vec Ji, arma::vec a_i, arma::vec b_i, double a_star, double b_star);
double Update_mu_star(double I, arma::vec T_i, arma::vec Ji);
double Update_D(double I, double a_star, double b_star);
double Update_a_star(double a_star_0, double I);
double Update_b_star(double I, arma::vec Ji, double b_star_0, double mu_star, double D,
                     arma::vec T_i, arma::vec V_i, arma::vec unique_indicator);
arma::vec Update_a_i(arma::vec Ji_Mc,arma::vec a_i_0);
arma::vec Update_b_i(arma::vec b_i_0, arma::vec T_i, arma::vec V_i,
              arma::vec m_ij, arma::vec S_ij, arma::vec group_indicator);


//Initialization
arma::vec Initialize_mu_ij(arma::vec xct, arma::vec nct);
arma::vec Initialize_z_ijk_case1(arma::vec mu_ij);
arma::vec Initialize_z_ijk_case0(arma::vec mu_ij);
arma::vec Initialize_m_ij(arma::vec xct, arma::vec nct);
arma::vec Initialize_S_ij(arma::vec xct, arma::vec nct);
arma::vec Initialize_T_i(arma::vec m_ij, arma::vec group_indicator, arma::vec Ji);
arma::vec Initialize_V_i(arma::vec Ji);
double Initialize_mu_star(arma::vec T_i, arma::vec unique_indicator);
double Initialize_D(double I);
double Initialize_a_star(arma::vec xct, arma::vec nct);
double Initialize_b_star(arma::vec xct, arma::vec nct);
arma::vec Initialize_a_i(arma::vec xct, arma::vec nct);
arma::vec Initialize_b_i(arma::vec xct, arma::vec nct);




// [[Rcpp::export]]
List Main_Function(int MaxIter, arma::vec xct,arma::vec nct,arma::vec tr_labels, bool display_progress=true) {
  int i;
  double Mc = 0;

  arma::vec unique_elements = arma::unique(tr_labels);
  double I = static_cast<int>(unique_elements.n_elem);

  //Unique identifier
  arma::vec unique_indicator = tr_labels;
  std::unordered_map<double, int> occurrences;
  for (std::size_t i = 0; i < tr_labels.n_elem; ++i) {
    if (occurrences[tr_labels(i)] == 0){
      unique_indicator(i) = 1;
    }else{
      unique_indicator(i) = 0;  // Set to 0 if it appears for the second time
    }
    occurrences[tr_labels(i)]++;
  }

  //Ji, Mc, and Ji_Mc
  arma::vec Ji(tr_labels.size());
  for (std::size_t i = 0; i < tr_labels.n_elem; ++i) {
    Ji(i) = occurrences[tr_labels(i)];
    if(occurrences[tr_labels(i)]==1){
      Mc++;
    }
  }

  arma::vec Ji_Mc = Ji;
  Ji_Mc.elem(arma::find(Ji_Mc==1)).fill(Mc);

  Progress p(MaxIter, display_progress);


  //initialize parameters
  arma::vec mu_ij = Initialize_mu_ij(xct, nct);
  arma::vec z_ijk_case1 = Initialize_z_ijk_case1(mu_ij);
  arma::vec z_ijk_case0 = Initialize_z_ijk_case0(mu_ij);
  arma::vec m_ij = Initialize_m_ij(xct, nct);
  arma::vec S_ij = Initialize_S_ij(xct, nct);
  arma::vec T_i = Initialize_T_i(m_ij, tr_labels, Ji);
  arma::vec V_i = Initialize_V_i(Ji);
  double mu_star = Initialize_mu_star(T_i, unique_indicator);
  double D = Initialize_D(I);
  double a_star_0 = Initialize_a_star(xct, nct);
  double b_star_0 = Initialize_b_star(xct, nct);
  arma::vec a_i_0 = Initialize_a_i(xct, nct);
  arma::vec b_i_0 = Initialize_b_i(xct, nct);

  double a_star = a_star_0;
  double b_star = b_star_0;
  arma::vec a_i = a_i_0;
  arma::vec b_i = b_i_0;

  for(i=0;i<MaxIter;i++){
    p.increment();

    mu_ij = Update_mu_ij(m_ij);
    z_ijk_case1 = Update_z_ijk_case1(mu_ij);
    z_ijk_case0 = Update_z_ijk_case0(mu_ij);
    m_ij = Update_m_ij(a_i, b_i, T_i, z_ijk_case1, z_ijk_case0, nct, xct);
    S_ij = Update_S_ij(a_i, b_i, nct);
    T_i = Update_T_i(Ji, tr_labels, a_i, b_i, m_ij, a_star, b_star, mu_star);
    V_i = Update_V_i(Ji, a_i, b_i, a_star, b_star);
    mu_star = Update_mu_star(I, T_i, Ji);
    D = Update_D(I, a_star, b_star);
    a_star = Update_a_star(a_star_0, I);
    b_star = Update_b_star(I, Ji, b_star_0, mu_star, D, T_i, V_i, unique_indicator);
    a_i = Update_a_i(Ji_Mc, a_i_0);
    b_i = Update_b_i(b_i_0, T_i, V_i, m_ij, S_ij, tr_labels);

  }

  arma::vec theta_ij = m_ij;
  arma::vec theta_i = T_i;
  arma::vec sigma_i = a_i / b_i;
  double tau2 = a_star / b_star;

  return(List::create(Rcpp::Named("z_ijk_case1")=z_ijk_case1,
                      Rcpp::Named("z_ijk_case0")=z_ijk_case0,
                      Rcpp::Named("theta_ij")=theta_ij,
                      Rcpp::Named("theta_i")=theta_i,
                      Rcpp::Named("sigma_i")=sigma_i,
                      Rcpp::Named("mu_star")=mu_star,
                      Rcpp::Named("tau2")=tau2
                      ));
}

//Update_Expectations
arma::vec Update_mu_ij(arma::vec m_ij){
  return(m_ij);
};

arma::vec Update_z_ijk_case1(arma::vec mu_ij){
  return(mu_ij + normpdf(-mu_ij) / (normcdf(mu_ij)));
};

arma::vec Update_z_ijk_case0(arma::vec mu_ij){
  return(mu_ij - normpdf(-mu_ij) / (normcdf(-mu_ij)));
};

arma::vec Update_m_ij(arma::vec a_i, arma::vec b_i, arma::vec T_i,
                      arma::vec z_ijk_case1, arma::vec z_ijk_case0, arma::vec nct, arma::vec xct){
  return((1 / (a_i / b_i + nct)) % (T_i % a_i / b_i + xct % z_ijk_case1 + (nct-xct) % z_ijk_case0));
};

arma::vec Update_S_ij(arma::vec a_i, arma::vec b_i, arma::vec nct){
  return(1 / (a_i / b_i + nct));
};

arma::vec Update_T_i(arma::vec Ji, arma::vec group_indicator, arma::vec a_i, arma::vec b_i, arma::vec m_ij,
                     double a_star, double b_star, double mu_star){
  std::map<int, double> groupSums;

  for (size_t i = 0; i < m_ij.size(); ++i) {
    groupSums[group_indicator[i]] += m_ij[i];
  }

  arma::vec result(m_ij.size());
  for (size_t i = 0; i < m_ij.size(); ++i) {
    result[i] = groupSums[group_indicator[i]];
  }

  return((1 / (a_star / b_star + Ji % a_i / b_i)) % (a_star / b_star * mu_star + result % a_i / b_i));
};

arma::vec Update_V_i(arma::vec Ji, arma::vec a_i, arma::vec b_i, double a_star, double b_star){
  return(1 / (a_star / b_star + Ji % a_i / b_i));
};

double Update_mu_star(double I, arma::vec T_i, arma::vec Ji){
  arma::vec temp1 = T_i / Ji;
  return(arma::sum(temp1) / (I));
};

double Update_D(double I, double a_star, double b_star){
  return(1 / (I * a_star / b_star));
};

double Update_a_star(double a_star_0, double I){
  return(a_star_0+I/2);
};

double Update_b_star(double I, arma::vec Ji, double b_star_0, double mu_star, double D,
                     arma::vec T_i, arma::vec V_i, arma::vec unique_indicator){
  arma::vec temp1 = T_i / Ji;
  arma::vec temp2 = (pow(T_i,2) + V_i) / Ji;

  return(b_star_0 + 0.5 * (I * (pow(mu_star,2) + D) - 2 * mu_star * (sum(temp1)) + sum(temp2)));
};

arma::vec Update_a_i(arma::vec Ji_Mc,arma::vec a_i_0){
  return(a_i_0+Ji_Mc/2);
};

arma::vec Update_b_i(arma::vec b_i_0, arma::vec T_i, arma::vec V_i,
              arma::vec m_ij, arma::vec S_ij, arma::vec group_indicator){

  arma::vec temp = pow(T_i,2) + V_i - 2 * T_i % m_ij + pow(m_ij, 2) + S_ij;

  std::map<int, double> groupSums;

  for (size_t i = 0; i < temp.size(); ++i) {
    groupSums[group_indicator[i]] += temp[i];
  }

  arma::vec result(temp.size());
  for (size_t i = 0; i < temp.size(); ++i) {
    result[i] = groupSums[group_indicator[i]];
  }

  return(b_i_0 + 0.5 * result);

};
//Initialize parameters
//Initialization
arma::vec Initialize_mu_ij(arma::vec xct, arma::vec nct){
  return(log((xct+0.001)/(nct-xct+0.001)));
}

arma::vec Initialize_z_ijk_case1(arma::vec mu_ij){
  return(mu_ij + normpdf(-mu_ij) / (1 - normcdf(-mu_ij)));
}

arma::vec Initialize_z_ijk_case0(arma::vec mu_ij){
  return(mu_ij - normpdf(-mu_ij) / (1 - normcdf(-mu_ij)));
}

arma::vec Initialize_m_ij(arma::vec xct, arma::vec nct){
  return(log((xct+0.001)/(nct-xct+0.001)));
}

arma::vec Initialize_S_ij(arma::vec xct, arma::vec nct){
  return(1 / (1+nct));
}

arma::vec Initialize_T_i(arma::vec m_ij, arma::vec group_indicator, arma::vec Ji){
  arma::vec temp(m_ij.size());

  std::map<int, double> groupSums;
  for (size_t i = 0; i < m_ij.size(); ++i) {
    groupSums[group_indicator[i]] += m_ij[i];
  }

  arma::vec result(m_ij.size());
  for (size_t i = 0; i < m_ij.size(); ++i) {
    result[i] = groupSums[group_indicator[i]];
  }

  return(result / Ji);
}

arma::vec Initialize_V_i(arma::vec Ji){
  return(1 / (1+Ji));
}

double Initialize_mu_star(arma::vec T_i, arma::vec unique_indicator){
  double sum_mu_star = arma::sum(T_i % unique_indicator) / arma::sum(unique_indicator);
  return(sum_mu_star);
}

double Initialize_D(double I){
  return(1 / I);
}

double Initialize_a_star(arma::vec xct, arma::vec nct){
  return(100);
}

double Initialize_b_star(arma::vec xct, arma::vec nct){
  return(100);
}

arma::vec Initialize_a_i(arma::vec xct, arma::vec nct){
  arma::vec a_i_0(xct.size());
  a_i_0.fill(100);
  return(a_i_0);
}

arma::vec Initialize_b_i(arma::vec xct, arma::vec nct){
  arma::vec b_i_0(xct.size());
  b_i_0.fill(100);
  return(b_i_0);
}
