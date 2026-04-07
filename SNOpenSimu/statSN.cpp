#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec compute_stat(arma::vec x, int T, String type, arma::vec weight) {
  int N = x.n_elem;
  int nT = N - T;
  
  // 核心前缀和 S[i] = X_0 + ... + X_{i-1}
  arma::vec S(N + 1, fill::zeros);
  for (int i = 0; i < N; i++) {
    S[i+1] = S[i] + x[i];
  }
  
  arma::vec stat_val(nT, fill::zeros);
  
  if (type == "snemk") {
    // 全局前缀和，用于 O(1) 计算自正则矩阵
    arma::vec P_j(N + 1, fill::zeros), P_j2(N + 1, fill::zeros);
    arma::vec P_S(N + 1, fill::zeros), P_S2(N + 1, fill::zeros), P_jS(N + 1, fill::zeros);
    
    for (int i = 1; i <= N; i++) {
      double j = i, Si = S[i];
      P_j[i]  = P_j[i-1]  + j;
      P_j2[i] = P_j2[i-1] + j * j;
      P_S[i]  = P_S[i-1]  + Si;
      P_S2[i] = P_S2[i-1] + Si * Si;
      P_jS[i] = P_jS[i-1] + j * Si;
    }
    
    // O(N) 预计算左侧样本自正则矩阵 V_L[z]
    arma::vec V_L(N + 1, fill::zeros);
    for (int z = T; z <= N; z++) {
      double Sz = S[z];
      V_L[z] = (double)z * z * P_S2[z-1] - 2.0 * z * Sz * P_jS[z-1] + Sz * Sz * P_j2[z-1];
    }
    
    double T15 = std::pow((double)T, 1.5);
    
    for (int k = T + 1; k <= N; k++) {
      double max_sn = 0.0;
      for (int z = T; z < k; z++) {
        double V_R = 0.0;
        // O(1) 极速计算右侧样本自正则矩阵 V_R(z,k)
        if (k - 1 >= z + 1) {
          double sum_1 = (k - 1) - z;
          double sum_j = P_j[k-1] - P_j[z];
          double sum_j2 = P_j2[k-1] - P_j2[z];
          double sum_S = P_S[k-1] - P_S[z];
          double sum_S2 = P_S2[k-1] - P_S2[z];
          double sum_jS = P_jS[k-1] - P_jS[z];
          
          double Sz = S[z];
          double sum_A2 = sum_S2 - 2.0 * Sz * sum_S + Sz * Sz * sum_1;
          double sum_AB = sum_jS - (double)z * sum_S - Sz * sum_j + (double)z * Sz * sum_1;
          double sum_B2 = sum_j2 - 2.0 * z * sum_j + (double)z * z * sum_1;
          
          double k_z = k - z, S_k_z = S[k] - S[z];
          V_R = k_z * k_z * sum_A2 - 2.0 * k_z * S_k_z * sum_AB + S_k_z * S_k_z * sum_B2;
        }
        
        double v_total = V_L[z] + V_R;
        if (v_total > 1e-12) {
          double diff = S[z] / z - (S[k] - S[z]) / (k - z);
          double num = T15 * (k - z) * std::abs(diff);
          double sn = num / std::sqrt(v_total);
          if (sn > max_sn) max_sn = sn;
        }
      }
      stat_val(k - T - 1) = max_sn * weight(k - T - 1);
    }
  }
  return stat_val;
}

