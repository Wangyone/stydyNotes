rm(list = ls())

library(foreach)
library(doParallel)
library(doSNOW)

# ==========================================
# 创新：解耦极值自正则 (Decoupled Min-Max SN) 极限分布
# ==========================================
ad_func_decoupled_minmax <- function(N, delta) {
  # 1. 极速生成标准布朗运动路径
  y <- (0:N) / N
  dy <- 1 / N
  B <- c(0, cumsum(rnorm(N, mean = 0, sd = 1/sqrt(N))))
  
  # 理论无穷域映射：tau 搜索到 0.99 即可捕获全局最大值
  tau_max <- 0.99
  i_max <- floor(N * tau_max) 
  
  # 2. 预计算所有需要的前缀和
  w <- 1 / (1 - y)^4
  w[(i_max + 1):(N + 1)] <- 0 # 保护尾部前缀和精度
  
  P_S0  <- cumsum(w) * dy
  P_S1  <- cumsum(w * y) * dy
  P_S2  <- cumsum(w * y^2) * dy
  P_SB  <- cumsum(w * B) * dy
  P_SB2 <- cumsum(w * B^2) * dy
  P_SyB <- cumsum(w * y * B) * dy
  
  global_max_stat <- 0
  
  # 3. 向量化搜索全局上确界
  # 至少需要几个点来计算内部方差，从 i = 5 开始
  for (i in 5:i_max) {
    tau <- y[i]
    B_tau <- B[i]
    
    # 候选游走点 v (对应分子中的 v1 和分母中的 v2，向量化统一计算)
    j_vec <- 2:(i-1)
    v <- y[j_vec]
    B_v <- B[j_vec]
    
    # ==========================================
    # 【解耦步骤 1：独立最大化分子】
    # ==========================================
    num_array <- (1 - v) * abs(B_tau - B_v)
    max_num <- max(num_array) 
    
    # ==========================================
    # 【解耦步骤 2：独立最小化分母】
    # ==========================================
    # 计算左侧积分 J_L(v)
    J_L <- pmax(P_SB2[j_vec] - 2 * B_v * P_SyB[j_vec] + B_v^2 * P_S2[j_vec], 0)
    
    # 计算右侧积分 J_R(v, tau)
    A <- tau - v
    B_diff <- B_tau - B_v
    
    I_B2 <- (P_SB2[i] - P_SB2[j_vec]) - 2 * B_v * (P_SB[i] - P_SB[j_vec]) + B_v^2 * (P_S0[i] - P_S0[j_vec])
    I_y2 <- (P_S2[i] - P_S2[j_vec]) - 2 * v * (P_S1[i] - P_S1[j_vec]) + v^2 * (P_S0[i] - P_S0[j_vec])
    I_By <- (P_SyB[i] - P_SyB[j_vec]) - v * (P_SB[i] - P_SB[j_vec]) - B_v * (P_S1[i] - P_S1[j_vec]) + B_v * v * (P_S0[i] - P_S0[j_vec])
    
    J_R <- pmax(A^2 * I_B2 - 2 * A * B_diff * I_By + B_diff^2 * I_y2, 0)
    
    denom_array <- sqrt((1 - tau)^2 * J_L + J_R)
    
    # 提取有效的最小分母 (兜底保护)
    valid_denom <- denom_array[denom_array > 1e-8]
    
    if (length(valid_denom) > 0) {
      min_denom <- min(valid_denom)
      
      # ==========================================
      # 【解耦步骤 3：组合最新统计算子】
      # ==========================================
      if (min_denom > 1e-8) {
        current_stat <- ((1 - tau)^delta) * max_num / min_denom
        if (current_stat > global_max_stat) {
          global_max_stat <- current_stat
        }
      }
    }
  }
  return(global_max_stat)
}

# ==========================================
# 并行化分位点模拟
# ==========================================
alpha <- 0.05
p <- 1 - alpha
delta <- 0.25

N0 <- 1000     # 维持较高网格密度
nrep <- 2000   

cores_to_use <- max(1, parallel::detectCores() - 1)
cat("正在启动解耦自正则 (Decoupled Min-Max SN) 渐近极限模拟...\n")
cat("启用的核心数:", cores_to_use, "\n")

cl <- makeCluster(cores_to_use)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

startTime <- Sys.time()

# 显式导出解耦函数
ad <- foreach(s = 1:nrep, .combine = 'c', .options.snow = opts, .export = "ad_func_decoupled_minmax") %dopar% {
  ad_func_decoupled_minmax(N0, delta)
}

close(pb)
stopCluster(cl)

endTime <- Sys.time()
cat("\n计算完成！总耗时:", round(difftime(endTime, startTime, units='mins'), 2), "分钟\n")

# ==========================================
# 输出渐近分布 0.95 临界值
# ==========================================
cat("\n===== Decoupled Min-Max SN 渐近分布模拟结果 =====\n")
final_quantile_inf <- round(quantile(ad, p, na.rm=TRUE), 4)
cat("理论无穷域 0.95 临界值: ", final_quantile_inf, "\n", sep="")

hist(ad, breaks = 50, probability = TRUE, main = "Decoupled Min-Max SN Asymptotic Distribution",
     xlab = "Statistic Value", ylab = "Density", col="lightgreen", border="white")
lines(density(ad), col = "darkgreen", lwd = 2)
abline(v = final_quantile_inf, col = "red", lwd = 2, lty = 2)