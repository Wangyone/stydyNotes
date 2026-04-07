# 清理环境变量
rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(foreach)

setwd("D:/Desktop/newPaper/SNOpenSimu")

# ==========================================
# 0. 辅助生成数据函数
# ==========================================
genData <- function(m, Tm, model, cp_frac = 0.5, delta = 0.5) {
    burnin <- ceiling(m / 2)
    obs <- burnin + m + Tm
    x <- numeric(obs)
    if (model == 1 || model == 101) { x <- rnorm(obs, mean = 0, sd = 1) } 
    else if (model == 2 || model == 102) { x <- arima.sim(model = list(ar = 0.3), n = obs) } 
    else if (model == 3 || model == 103) { x <- arima.sim(model = list(ar = 0.5), n = obs) } 
    else if (model == 4 || model == 104) { x <- arima.sim(model = list(ar = 0.7), n = obs) } 
    else if (model == 5 || model == 105) { x <- arima.sim(model = list(ar = 0.9), n = obs) } 
    
    if (model %in% 101:105) {
        change.pt <- burnin + m + floor(Tm * cp_frac)
        x[change.pt:obs] <- x[change.pt:obs] + delta
    }
    return(x[(burnin + 1):obs])
}

# 提取 C++ 代码字符串 (供 Stage 2 独立沙箱编译使用)
cpp_code_string <- "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

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
    P_W[i] = P_W[i-1] + W_i;     P_Wi[i] = P_Wi[i-1] + W_i * i_dbl;
    P_Wi2[i] = P_Wi2[i-1] + W_i * i_dbl * i_dbl;
    P_WS[i] = P_WS[i-1] + W_i * Si; P_WiS[i] = P_WiS[i-1] + W_i * i_dbl * Si;
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
"

# ==========================================
# 1. 纯 R 语言极限泛函模拟 (极致向量化优化)
# ==========================================
cat(">>> [阶段 1/2] 正在用纯 R 语言多线程模拟连续布朗泛函渐近临界值...\n")
M_cv <- 10000 
N_grid <- 2000
eps <- 0.05

cores_to_use <- max(1, parallel::detectCores() - 1)
cl_cv <- makeCluster(cores_to_use)
registerDoParallel(cl_cv)

# 使用纯 R 语言的向量化计算单条轨迹，无编译风险
results_cv <- foreach(i = 1:M_cv, .combine = c) %dopar% {
    dt <- 1.0 / N_grid
    y <- seq(0, 1, length.out = N_grid + 1)
    dW <- rnorm(N_grid, mean = 0, sd = sqrt(dt))
    B <- c(0, cumsum(dW))
    
    # 前缀和计算积分
    P_B2 <- cumsum(B^2 * dt); P_y2 <- cumsum(y^2 * dt); P_yB <- cumsum(y * B * dt)
    P_B  <- cumsum(B * dt);   P_y  <- cumsum(y * dt)
    
    path_sup <- 0
    ix_start <- N_grid / 2 + 1
    ix_end <- round(0.99 * N_grid) + 1
    
    for (ix in ix_start:ix_end) {
        y_x <- y[ix]
        x <- y_x / (1.0 - y_x)
        
        u_min <- 1.0 + eps * (x - 1.0)
        u_max <- x - eps * (x - 1.0)
        y_min <- u_min / (1.0 + u_min)
        y_max <- u_max / (1.0 + u_max)
        
        imin <- ceiling(y_min * N_grid) + 1
        imax <- floor(y_max * N_grid) + 1
        if (imin >= imax) next
        
        idx <- imin:imax
        y_sub <- y[idx]; B_sub <- B[idx]
        
        # 向量化分子极大值
        max_num <- max(abs(B_sub - (y_sub / y_x) * B[ix]))
        
        # 向量化分母极小值
        B2_y2 <- B_sub / y_sub
        VL <- P_B2[idx] - 2.0 * B2_y2 * P_yB[idx] + (B2_y2^2) * P_y2[idx]
        
        B_star_x <- B[ix]*(1.0 - y_sub) - B_sub*(1.0 - y_x)
        c1 <- 1.0 - y_sub
        c2 <- B_sub - B_star_x / (y_x - y_sub)
        c3 <- -B_sub + y_sub * B_star_x / (y_x - y_sub)
        
        int_D2 <- (c1^2)*(P_B2[ix] - P_B2[idx]) + (c2^2)*(P_y2[ix] - P_y2[idx]) + (c3^2)*(y_x - y_sub) + 
            2*c1*c2*(P_yB[ix] - P_yB[idx]) + 2*c1*c3*(P_B[ix] - P_B[idx]) + 2*c2*c3*(P_y[ix] - P_y[idx])
        VR <- int_D2 / ((1.0 - y_sub)^2)
        
        denom <- sqrt(pmax(0, VL) + pmax(0, VR))
        valid <- denom > 1e-10
        if (any(valid)) {
            stat <- max_num / min(denom[valid])
            if (stat > path_sup) path_sup <- stat
        }
    }
    return(path_sup)
}
stopCluster(cl_cv)

cv_95_true <- quantile(results_cv, 0.95)
cat(">>> 数学之美：R 语言布朗近似 95% 临界值 (cv_95) =", round(cv_95_true, 4), "<<<\n\n")

# ==========================================
# 2. 实证仿真评估 (物理隔绝 Rcpp 并行冲突)
# ==========================================
cat(">>> [阶段 2/2] 启动 H0/H1 实证仿真 (安全沙箱编译模式)...\n")
T_0 <- 100
Tm <- 30 * T_0
nsim <- 1000
models <- c(1:5, 101:105)

cl_sim <- makeCluster(cores_to_use)
registerDoParallel(cl_sim)

results <- foreach(mod = models, .combine = rbind, 
                   .packages = c("Rcpp", "RcppArmadillo")) %dopar% {
                       
                       # 【终极防弹机制】：每个核心在自己的临时沙箱写入并编译 C++，绝不串联共享！
                       tmp_file <- tempfile(fileext = ".cpp")
                       writeLines(cpp_code_string, tmp_file)
                       Rcpp::sourceCpp(tmp_file)
                       
                       rejections <- 0
                       for (s in 1:nsim) {
                           yt <- genData(m = T_0, Tm = Tm, model = mod, cp_frac = 0.1, delta = 1)
                           sn_stats <- compute_newSN(yt, T_0, 0.05)
                           if (max(sn_stats) > cv_95_true) rejections <- rejections + 1
                       }
                       
                       # 扫尾清理沙箱临时文件
                       unlink(tmp_file)
                       
                       data.frame(Model = mod, Rejection_Rate = rejections / nsim)
                   }
stopCluster(cl_sim)

# ==========================================
# 3. 输出顶刊级别的最终仿真结果
# ==========================================
results$Type <- ifelse(results$Model <= 5, "Type I Error (H0)", "Power (H1)")
results$AR_Coef <- rep(c(0, 0.3, 0.5, 0.7, 0.9), 2)

cat("\n======================================================\n")
cat("      Decoupled Time-Discounted SN 理论闭环验证表\n")
cat("======================================================\n")
print(results[, c("Type", "AR_Coef", "Model", "Rejection_Rate")], row.names = FALSE)
cat("======================================================\n")
