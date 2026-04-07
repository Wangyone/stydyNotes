rm(list = ls())
library(foreach)
library(doParallel)
library(doSNOW)

# ==========================================
# 终极圣杯：全局开放式 (Open-End) 连续极限泛函
# ==========================================
sim_open_end_functional <- function(N_grid, eps = 0.05, gamma = 0.25) {
  # 1. 在 [0, 1] 域生成标准布朗运动
  y <- (0:N_grid) / N_grid
  dy <- 1 / N_grid
  B <- c(0, cumsum(rnorm(N_grid, mean = 0, sd = 1/sqrt(N_grid))))
  
  # 预计算左侧积分
  P_B2 <- cumsum(B^2) * dy
  P_yB <- cumsum(y * B) * dy
  P_y2 <- cumsum(y^2) * dy
  I_L_vec <- numeric(N_grid + 1)
  for (j in 2:(N_grid + 1)) {
    yj <- y[j]; Bj <- B[j]
    I_L_vec[j] <- P_B2[j] - 2 * (Bj / yj) * P_yB[j] + (Bj / yj)^2 * P_y2[j]
  }
  
  global_max_stat <- 0
  
  # 2. 全局开放式搜索：y_x 从 0.51 搜索到 0.999 (即物理时间 x 从 1.04 搜索到 1000倍!)
  idx_start <- floor(0.51 * N_grid) 
  idx_end   <- floor(0.999 * N_grid)
  
  for (i in idx_start:idx_end) {
    yx <- y[i]
    Bx <- B[i]
    x <- yx / (1 - yx) # 映射到真实的物理相对时间
    
    # ---------------------------------------------------------
    # 物理时间的 Sup-Wald 修剪映射
    # 候选点 u 必须在监测期的内部：[1 + eps*(x-1), x - eps*(x-1)]
    # ---------------------------------------------------------
    u_min <- 1 + eps * (x - 1)
    u_max <- x - eps * (x - 1)
    if (u_max <= u_min) next
    
    # 映射回 y 域寻找网格索引
    y_min <- u_min / (1 + u_min)
    y_max <- u_max / (1 + u_max)
    idx_cand <- max(1, round(y_min * N_grid)) : max(1, round(y_max * N_grid))
    
    if (length(idx_cand) < 1) next
    y_cand <- y[idx_cand]
    B_cand <- B[idx_cand]
    
    # 计算无奇点的分子
    num_array <- (1 / (1 - y_cand)) * abs(B_cand - (y_cand / yx) * Bx)
    max_num <- max(num_array)
    
    # 计算无奇点的右侧积分分母
    I_R_array <- sapply(idx_cand, function(j) {
      y2 <- y[j]; B2 <- B[j]
      y_sub <- y[j:i]; B_sub <- B[j:i]
      
      B_star   <- B_sub * (1 - y2) - B2 * (1 - y_sub)
      B_star_x <- Bx * (1 - y2)    - B2 * (1 - yx)
      integrand <- ( B_star - ((y_sub - y2) / (yx - y2)) * B_star_x )^2
      
      sum(integrand) * dy / (1 - y2)^2
    })
    
    denom_array <- sqrt(I_L_vec[idx_cand] + I_R_array)
    valid_denom <- denom_array[denom_array > 1e-8]
    if (length(valid_denom) == 0) next
    min_denom <- min(valid_denom)
    
    # ---------------------------------------------------------
    # 施加开放式极限的定海神针：时间权重 x^(-gamma)
    # ---------------------------------------------------------
    current_stat <- (x^(-gamma)) * max_num / min_denom
    
    if (current_stat > global_max_stat) {
      global_max_stat <- current_stat
    }
  }
  return(global_max_stat)
}

# ==========================================
# 并行化模拟：提取全局临界值
# ==========================================
N_grid <- 1000    # 网格密度
nrep <- 1000      # 1500次足以稳定
gamma_val <- 0.25 # 为了在无穷域收敛，gamma 必须 > 0

cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(cores)
registerDoSNOW(cl)

cat("正在模拟真正的全局开放式极限 (监控期 -> 无穷)...\n")
pb <- txtProgressBar(max = nrep, style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

startTime <- Sys.time()
dist_universal <- foreach(s = 1:nrep, .combine = 'c', .options.snow = opts) %dopar% {
  sim_open_end_functional(N_grid = N_grid, eps = 0.05, gamma = gamma_val)
}
stopCluster(cl)
endTime <- Sys.time()

# 提取并输出终极临界值
cv_95_universal <- round(quantile(dist_universal, 0.95), 4)
cat("\n=========================================\n")
cat(sprintf("在 Gamma = %s 下，全局开放式 0.95 临界值为: %s\n", gamma_val, cv_95_universal))
cat("=========================================\n")
