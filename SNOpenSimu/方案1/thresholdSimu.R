rm(list = ls())

library(foreach)
library(doParallel)
library(doSNOW)

# ==========================================
# 纯理论无穷域极限分布：Infinite Horizon Split-Sample SN
# ==========================================
ad_func_infinity <- function(N, delta) {
  # 1. 生成布朗运动路径
  y <- (0:N) / N
  dy <- 1 / N
  B <- c(0, cumsum(rnorm(N, mean = 0, sd = 1/sqrt(N))))
  
  # 【理论与数值的完美妥协】
  # 既然 tau -> 1 时统计量必趋于 0，我们在 0.99 处截断以捕捉全局 supremum
  # tau = 0.99 对应真实监测期是历史期的 99 倍，这已经是事实上的“无穷远”
  tau_max <- 0.99
  i_max <- floor(N * tau_max) 
  
  # 2. 预计算所有需要的前缀和
  w <- 1 / (1 - y)^4
  w[(i_max + 1):(N + 1)] <- 0 # 防止尾部大数吃小数，保护前缀和精度
  
  P_S0  <- cumsum(w) * dy
  P_S1  <- cumsum(w * y) * dy
  P_S2  <- cumsum(w * y^2) * dy
  P_SB  <- cumsum(w * B) * dy
  P_SB2 <- cumsum(w * B^2) * dy
  P_SyB <- cumsum(w * y * B) * dy
  
  max_val <- 0
  
  # 3. 向量化搜索全局上确界
  for (i in 3:i_max) {
    tau <- y[i]
    B_tau <- B[i]
    
    j <- 2:(i-1)
    v <- y[j]
    B_v <- B[j]
    
    # 极速计算左侧与右侧积分，并用 pmax 消除由于浮点展开带来的微小负数误差
    J_L <- pmax(P_SB2[j] - 2 * B_v * P_SyB[j] + B_v^2 * P_S2[j], 0)
    
    A <- tau - v
    B_diff <- B_tau - B_v
    
    I_B2 <- (P_SB2[i] - P_SB2[j]) - 2 * B_v * (P_SB[i] - P_SB[j]) + B_v^2 * (P_S0[i] - P_S0[j])
    I_y2 <- (P_S2[i] - P_S2[j]) - 2 * v * (P_S1[i] - P_S1[j]) + v^2 * (P_S0[i] - P_S0[j])
    I_By <- (P_SyB[i] - P_SyB[j]) - v * (P_SB[i] - P_SB[j]) - B_v * (P_S1[i] - P_S1[j]) + B_v * v * (P_S0[i] - P_S0[j])
    
    J_R <- pmax(A^2 * I_B2 - 2 * A * B_diff * I_By + B_diff^2 * I_y2, 0)
    
    denom <- sqrt((1 - tau)^2 * J_L + J_R)
    
    valid <- denom > 1e-10
    if (any(valid)) {
      val <- ((1 - tau)^delta) * (1 - v[valid]) * abs(B_tau - B_v[valid]) / denom[valid]
      current_max <- max(val)
      if (current_max > max_val) {
        max_val <- current_max
      }
    }
  }
  return(max_val)
}

# ==========================================
# 并行化分位点模拟
# ==========================================
alpha <- 0.05
p <- 1 - alpha
delta <- 0

N0 <- 1500     # 稍微增加网格密度，让无穷域映射更加精细
nrep <- 5000   

cores_to_use <- max(1, parallel::detectCores() - 1)
cat("正在启动纯无穷域渐近极限模拟，使用核心数:", cores_to_use, "\n")

cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

startTime <- Sys.time()

ad <- foreach(s = 1:nrep, .combine = 'c', .options.snow = opts, .export = "ad_func_infinity") %dopar% {
  ad_func_infinity(N0, delta)
}

close(pb)
stopCluster(cl)

endTime <- Sys.time()
cat("\n计算完成！总耗时:", round(difftime(endTime, startTime, units='mins'), 2), "分钟\n")

# ==========================================
# 输出无穷域渐近分布 0.95 临界值
# ==========================================
cat("\n===== 纯无穷域 [0, ∞) 渐近分布模拟结果 =====\n")
final_quantile_inf <- round(quantile(ad, p, na.rm=TRUE), 4)
cat("理论无穷域 0.95 临界值: ", final_quantile_inf, "\n", sep="")

hist(ad, breaks = 50, probability = TRUE, main = "Infinite Horizon Asymptotic Distribution",
     xlab = "Statistic Value", ylab = "Density", col="lightblue", border="white")
lines(density(ad), col = "darkblue", lwd = 2)
abline(v = final_quantile_inf, col = "red", lwd = 2, lty = 2)