
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec compute_stat_discounted(arma::vec x, int T, arma::vec weight, double eps = 0.05) {
  int N = x.n_elem;
  int nT = N - T;
  double m = (double)T;
  
  // 1. 核心前缀和
  arma::vec S(N + 1, fill::zeros);
  for (int i = 0; i < N; i++) { S[i+1] = S[i] + x[i]; }
  
  arma::vec stat_val(nT, fill::zeros);
  
  // 2. 预计算：时间折现权重及其衍生前缀和 (O(N))
  arma::vec P_W(N+1, fill::zeros), P_Wi(N+1, fill::zeros), P_Wi2(N+1, fill::zeros);
  arma::vec P_WS(N+1, fill::zeros), P_WiS(N+1, fill::zeros), P_WS2(N+1, fill::zeros);
  
  for (int i = 1; i <= N; i++) {
    // 核心：时间折现权重 W_i = 1 / (1 + i/m)^4
    double W_i = 1.0 / std::pow(1.0 + (double)i / m, 4.0);
    double i_dbl = (double)i;
    double Si = S[i];
    
    P_W[i] = P_W[i-1] + W_i;
    P_Wi[i] = P_Wi[i-1] + W_i * i_dbl;
    P_Wi2[i] = P_Wi2[i-1] + W_i * i_dbl * i_dbl;
    P_WS[i] = P_WS[i-1] + W_i * Si;
    P_WiS[i] = P_WiS[i-1] + W_i * i_dbl * Si;
    P_WS2[i] = P_WS2[i-1] + W_i * Si * Si;
  }
  
  double sqrt_m = std::sqrt(m);
  
  // 3. 核心监测循环
  for (int k = T + 1; k <= N; k++) {
    // 内部修剪：由 eps 决定的候选搜索范围
    int lower_bound = std::max(1, (int)(eps * k));
    int upper_bound = (int)((1.0 - eps) * k);
    
    if (lower_bound >= upper_bound) {
      stat_val[k - T - 1] = 0.0;
      continue;
    }
    
    double max_num = -1.0;
    double min_denom = 1e300;
    
    // 解耦 A：极大化分子 (Numerator)
    for (int z = lower_bound; z <= upper_bound; z++) {
      double num = std::abs(S[z] - ((double)z / k) * S[k]);
      if (num > max_num) max_num = num;
    }
    
    // 解耦 B：极小化分母 (Denominator)
    for (int v = lower_bound; v <= upper_bound; v++) {
      // 左侧加权内部方差 VL
      double C_L = -S[v] / (double)v;
      double VL_raw = P_WS2[v] + C_L * C_L * P_Wi2[v] + 2.0 * C_L * P_WiS[v];
      
      // 右侧加权内部方差 VR
      double C_R = (S[k] - S[v]) / (double)(k - v);
      double A_coef = -C_R; 
      double B_coef = C_R * v - S[v];
      
      double d_WS2 = P_WS2[k] - P_WS2[v];
      double d_Wi2 = P_Wi2[k] - P_Wi2[v];
      double d_W   = P_W[k]   - P_W[v];
      double d_WiS = P_WiS[k] - P_WiS[v];
      double d_WS  = P_WS[k]  - P_WS[v];
      double d_Wi  = P_Wi[k]  - P_Wi[v];
      
      double VR_raw = d_WS2 + A_coef*A_coef*d_Wi2 + B_coef*B_coef*d_W + 
                      2.0*A_coef*d_WiS + 2.0*B_coef*d_WS + 2.0*A_coef*B_coef*d_Wi;
      
      double denom_raw = std::max(0.0, VL_raw) + std::max(0.0, VR_raw);
      
      if (denom_raw > 1e-12) {
        double denom = std::sqrt(denom_raw);
        if (denom < min_denom) min_denom = denom;
      }
    }
    
    // 解耦 C：拼合
    if (min_denom < 1e299 && min_denom > 1e-12) {
      // 这里的 sqrt_m 是物理量级对齐的关键
      stat_val[k - T - 1] = weight[k - T - 1] * sqrt_m * max_num / min_denom;
    } else {
      stat_val[k - T - 1] = 0.0;
    }
  } // <--- 这里是之前缺失的关键闭合括号！
  
  return stat_val;
}

