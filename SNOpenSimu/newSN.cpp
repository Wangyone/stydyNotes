
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// =======================================================
// 模块 A: 纯布朗运动极限泛函模拟器 (The True Asymptotic Limit)
// =======================================================
// [[Rcpp::export]]
arma::vec simulate_limit_functional(int M, int N, double eps = 0.05) {
    arma::vec sup_stats(M, fill::zeros);
    double dt = 1.0 / N;
    arma::vec y = arma::linspace(0.0, 1.0, N + 1);
    
    for(int i = 0; i < M; i++) {
        // 生成纯种连续布朗运动路径 B(y)
        arma::vec dW = arma::randn(N) * std::sqrt(dt);
        arma::vec B(N + 1, fill::zeros);
        for(int j = 0; j < N; j++) B[j+1] = B[j] + dW[j];
        
        // 代数降维：O(1) 积分所需的解析前缀和
        arma::vec P_B2(N+1, fill::zeros), P_y2(N+1, fill::zeros), P_yB(N+1, fill::zeros);
        arma::vec P_B(N+1, fill::zeros), P_y(N+1, fill::zeros);
        for(int j = 1; j <= N; j++) {
            P_B2[j] = P_B2[j-1] + B[j]*B[j]*dt;
            P_y2[j] = P_y2[j-1] + y[j]*y[j]*dt;
            P_yB[j] = P_yB[j-1] + y[j]*B[j]*dt;
            P_B[j]  = P_B[j-1]  + B[j]*dt;
            P_y[j]  = P_y[j-1]  + y[j]*dt;
        }
        
        double path_sup = 0.0;
        int ix_start = N / 2; // y_x = 0.5 对应现实时间起点
        int ix_end = std::round(0.99 * N); // y_x -> 1 对应现实无穷远
        
        for(int ix = ix_start; ix <= ix_end; ix++) {
            double y_x = y[ix];
            double x = y_x / (1.0 - y_x);
            
            // 理论动态时空映射的修剪域
            double u_min = 1.0 + eps * (x - 1.0);
            double u_max = x - eps * (x - 1.0);
            double y_min = u_min / (1.0 + u_min);
            double y_max = u_max / (1.0 + u_max);
            
            int imin = std::ceil(y_min * N);
            int imax = std::floor(y_max * N);
            if(imin >= imax) continue;
            
            double max_num = -1.0;
            double min_denom = 1e300;
            
            // 极限分子：有界布朗桥距离
            for(int j = imin; j <= imax; j++) {
                double y1 = y[j];
                double num = std::abs(B[j] - (y1 / y_x) * B[ix]);
                if(num > max_num) max_num = num;
            }
            
            // 极限分母：化简后的完全无奇点连续积分
            for(int j = imin; j <= imax; j++) {
                double y2 = y[j];
                double B2_y2 = B[j] / y2;
                // V_L: 左极限方差
                double VL = P_B2[j] - 2.0 * B2_y2 * P_yB[j] + B2_y2 * B2_y2 * P_y2[j];
                
                // V_R: 右极限方差 (二次展开提取系数)
                double B_star_x = B[ix]*(1.0 - y2) - B[j]*(1.0 - y_x);
                double c1 = 1.0 - y2;
                double c2 = B[j] - B_star_x / (y_x - y2);
                double c3 = -B[j] + y2 * B_star_x / (y_x - y2);
                
                double d_B2 = P_B2[ix] - P_B2[j];
                double d_y2 = P_y2[ix] - P_y2[j];
                double d_1  = y_x - y2;
                double d_yB = P_yB[ix] - P_yB[j];
                double d_B  = P_B[ix] - P_B[j];
                double d_y  = P_y[ix] - P_y[j];
                
                double int_D2 = c1*c1*d_B2 + c2*c2*d_y2 + c3*c3*d_1 + 2.0*c1*c2*d_yB + 2.0*c1*c3*d_B + 2.0*c2*c3*d_y;
                double VR = int_D2 / ((1.0 - y2)*(1.0 - y2)); // (1-y_2)^{-2} 的雅可比乘积
                
                double denom_var = std::max(0.0, VL) + std::max(0.0, VR);
                if(denom_var > 1e-12) {
                    double denom = std::sqrt(denom_var);
                    if(denom < min_denom) min_denom = denom;
                }
            }
            
            if(min_denom < 1e299 && min_denom > 1e-12) {
                double stat = max_num / min_denom;
                if(stat > path_sup) path_sup = stat;
            }
        }
        sup_stats[i] = path_sup;
    }
    return sup_stats;
}

// =======================================================
// 模块 B: 有限样本实证统计量 (The Finite-Sample Statistic)
// =======================================================
// [[Rcpp::export]]
arma::vec compute_newSN(arma::vec x, int T, double eps = 0.05) {
  int N = x.n_elem;
  int nT = N - T;
  double m = (double)T;
  
  arma::vec S(N + 1, fill::zeros);
  for (int i = 0; i < N; i++) { S[i+1] = S[i] + x[i]; }
  arma::vec stat_val(nT, fill::zeros);
  
  arma::vec P_W(N+1, fill::zeros), P_Wi(N+1, fill::zeros), P_Wi2(N+1, fill::zeros);
  arma::vec P_WS(N+1, fill::zeros), P_WiS(N+1, fill::zeros), P_WS2(N+1, fill::zeros);
  for (int i = 1; i <= N; i++) {
    double W_i = 1.0 / std::pow(1.0 + (double)i / m, 4.0);
    double i_dbl = (double)i;
    double Si = S[i];
    P_W[i]   = P_W[i-1]   + W_i;
    P_Wi[i]  = P_Wi[i-1]  + W_i * i_dbl;
    P_Wi2[i] = P_Wi2[i-1] + W_i * i_dbl * i_dbl;
    P_WS[i]  = P_WS[i-1]  + W_i * Si;
    P_WiS[i] = P_WiS[i-1] + W_i * i_dbl * Si;
    P_WS2[i] = P_WS2[i-1] + W_i * Si * Si;
  }
  
  for (int k = T + 1; k <= N; k++) {
    int margin = std::max(1, (int)(eps * (k - T)));
    int lower_bound = T + margin;
    int upper_bound = k - margin;
    if (lower_bound >= upper_bound) { stat_val[k - T - 1] = 0.0; continue; }
    
    double max_num = -1.0;
    double min_denom = 1e300;
    
    for (int z = lower_bound; z <= upper_bound; z++) {
      double weight_z = 1.0 / (1.0 + (double)z / m);
      double num = weight_z * std::abs(S[z] - ((double)z / k) * S[k]);
      if (num > max_num) max_num = num;
    }
    for (int v = lower_bound; v <= upper_bound; v++) {
      double C_L = -S[v] / (double)v;
      double VL = P_WS2[v] + C_L * C_L * P_Wi2[v] + 2.0 * C_L * P_WiS[v];
      double C_R = (S[k] - S[v]) / (double)(k - v);
      double A_coef = -C_R; 
      double B_coef = C_R * v - S[v];
      
      double d_WS2 = P_WS2[k] - P_WS2[v], d_Wi2 = P_Wi2[k] - P_Wi2[v], d_W = P_W[k] - P_W[v];
      double d_WiS = P_WiS[k] - P_WiS[v], d_WS = P_WS[k] - P_WS[v], d_Wi = P_Wi[k] - P_Wi[v];
      
      double VR = d_WS2 + A_coef*A_coef*d_Wi2 + B_coef*B_coef*d_W + 
                  2.0*A_coef*d_WiS + 2.0*B_coef*d_WS + 2.0*A_coef*B_coef*d_Wi;
      double denom_var = std::max(0.0, VL) + std::max(0.0, VR);
      if (denom_var > 1e-12) {
        double denom = std::sqrt(denom_var);
        if (denom < min_denom) min_denom = denom;
      }
    }
    if (min_denom < 1e299 && min_denom > 1e-12) {
      stat_val[k - T - 1] = std::sqrt(m) * max_num / min_denom;
    }
  } 
  return stat_val;
}

