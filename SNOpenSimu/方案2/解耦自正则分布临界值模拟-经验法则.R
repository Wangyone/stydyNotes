rm(list = ls())

library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(foreach)

# ==========================================
# 利用高精度 C++ 引擎获取真实的解耦自正则临界值
# ==========================================
T_0 <- 300                 # 历史期 (模拟布朗运动的网格密度)
L0 <- 10                   # 监测期倍数
Tm <- L0 * T_0             
N <- T_0 + Tm              
nsim_cv <- 5000            # 模拟 5000 次确保分位数绝对稳定

# 权重函数设定 (假设 gamma = 0，如果你的主程序有衰减，这里必须保持一致)
gamma <- 0
t_vec <- (1:Tm) / T_0
weights <- (1 + t_vec)^(-gamma)

cores <- max(1, detectCores() - 1)
cl <- makeCluster(cores)
registerDoParallel(cl)

cat("正在使用 C++ 引擎进行大样本白噪声模拟 (获取真实临界值)...\n")
startTime <- Sys.time()

# 模拟原假设 H0 下的最大统计量分布
null_dist <- foreach(s = 1:nsim_cv, .combine = 'c', .packages = c("Rcpp", "RcppArmadillo")) %dopar% {
  Rcpp::sourceCpp("D:/Desktop/SNOpenSimu/statSN.cpp")
  
  # 生成纯标准正态白噪声 (在 FCLT 下严格收敛于布朗运动)
  sim_x <- rnorm(N, mean = 0, sd = 1)
  
  # 调用我们最新写的解耦算子
  sn_stats <- compute_stat(sim_x, T_0, "decoupled_snemk", weights)
  
  # 直接取全局最大值
  max(sn_stats)
}

stopCluster(cl)
endTime <- Sys.time()
cat("临界值模拟完成！耗时:", round(difftime(endTime, startTime, units = "mins"), 2), "分钟\n\n")

# ==========================================
# 输出真实临界值
# ==========================================
real_cv <- round(quantile(null_dist, 0.95), 4)
cat("==========================================\n")
cat("!!! 请将此临界值填入你的主模拟程序 !!!\n")
cat("Decoupled Min-Max SN 真实 0.95 临界值: ", real_cv, "\n")
cat("==========================================\n")

hist(null_dist, breaks = 50, probability = TRUE, col="lightgreen", border="white",
     main="True Asymptotic Distribution (via C++)")
lines(density(null_dist), col="darkgreen", lwd=2)
abline(v = real_cv, col="red", lwd=2, lty=2)