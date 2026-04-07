rm(list = ls())

library(foreach)
library(doParallel)
library(doSNOW)
library(Rcpp)
library(RcppArmadillo)

# ==========================================
# 1. 核心参数设定 (与你的主程序必须严格一致)
# ==========================================
T_0 <- 100                 # 历史期样本量 (m)
L0 <- 10                   # 监测期倍数
Tm <- L0 * T_0             # 监测期样本量
N <- T_0 + Tm              # 总样本量 (2200)

nsim_cv <- 5000            # 临界值模拟次数 (推荐 5000 次确保稳定)
weights <- rep(1, Tm)      # 权重函数保持恒定 1

# ==========================================
# 2. 并行生成经验临界值 (Empirical Critical Value)
# ==========================================
cores_to_use <- max(1, parallel::detectCores() - 1)
cat("正在启动并行计算，寻找真实临界值...\n")

cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nsim_cv, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

startTime <- Sys.time()

# 核心逻辑：直接向 C++ 喂入标准正态分布 (代表原假设 H0)
null_stats <- foreach(s = 1:nsim_cv, .combine = 'c', .options.snow = opts, .packages = c("Rcpp", "RcppArmadillo")) %dopar% {
  # 注意：这里必须是你正确的 cpp 文件路径！
  Rcpp::sourceCpp("D:/Desktop/SNOpenSimu/statSN.cpp")
  
  # 生成纯白噪声数据 (严格符合 H0)
  sim_x <- rnorm(N, mean = 0, sd = 1)
  
  # 调用 C++ 计算自正则统计量
  sn_stats <- compute_stat(sim_x, T_0, "snemk", weights)
  
  # 返回整个监测期内的全局最大值
  max(sn_stats)
}

close(pb)
stopCluster(cl)

endTime <- Sys.time()
cat("\n临界值计算完成！耗时:", round(difftime(endTime, startTime, units = "mins"), 2), "分钟\n")

# ==========================================
# 3. 输出真实的 95% 临界值
# ==========================================
final_quantile <- round(quantile(null_stats, 0.95, na.rm = TRUE), 4)

cat("\n========================================\n")
cat("!!! 请将此临界值复制到你的主程序中 !!!\n")
cat("真实 0.95 临界值 (final_quantile) = ", final_quantile, "\n")
cat("========================================\n")

# 绘制真实的极限分布直方图
hist(null_stats, breaks = 50, probability = TRUE, col="lightblue", border="white",
     main="True Empirical Null Distribution of Split-Sample SN", xlab="Statistic Value")
lines(density(null_stats), col="darkblue", lwd=2)
abline(v = final_quantile, col="red", lwd=2, lty=2)