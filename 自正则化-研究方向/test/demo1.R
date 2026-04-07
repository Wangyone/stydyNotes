# ==============================================================================
# V3.0: A Randomized Kernel-based Online Monitoring Procedure 
#       via Self-Normalized Stochastic Integrals
# ==============================================================================
library(ggplot2)
library(foreach)
library(doParallel)
library(MASS)

# ==========================================
# 阶段一：数据生成器 (支持强相依性、均值与方差突变)
# ==========================================
generate_online_stream <- function(m, T_max, model = "H0", rho = 0.8, 
                                   shift = 1.0, var_ratio = 2.5, cp_delay = 100, seed_val = NULL) {
    if (!is.null(seed_val)) set.seed(seed_val)
    
    n <- m + T_max
    Z <- numeric(n)
    err <- rnorm(n, 0, 1)
    
    Z[1] <- err[1]
    # 强 AR(1) 序列生成
    for (t in 2:n) {
        Z[t] <- rho * Z[t-1] + err[t]
    }
    
    # 变点注入 (发生在历史期 m 之后，延迟 cp_delay)
    cp_idx <- m + cp_delay
    if (model == "H1_mean") {
        Z[(cp_idx + 1):n] <- Z[(cp_idx + 1):n] + shift
    } else if (model == "H1_var") {
        Z[(cp_idx + 1):n] <- Z[(cp_idx + 1):n] * sqrt(var_ratio)
    }
    
    return(list(Z = Z, cp_idx = cp_idx))
}

# ==========================================
# 阶段二：V3.0 无截断 Itô 投影算子 (核心重构)
# ==========================================
# 1. 投影器初始化（永久抽样，锁定基底）
init_ito_projector <- function(K = 5, M = 500, sigma_w = 1, seed = 999) {
    set.seed(seed)
    # 1. 全空间测度抽样 u ~ w(u) (无截断实数轴)
    u_vec <- rnorm(M, mean = 0, sd = sigma_w) 
    # 2. 白噪声探针矩阵 dV_k (矩阵形式)
    eps_mat <- matrix(rnorm(M * K), nrow = M, ncol = K) 
    
    return(list(u = u_vec, eps = eps_mat, K = K, M = M))
}

# 2. 极速矩阵流式映射 (抛弃积分，拥抱矩阵乘法)
apply_ito_projection <- function(Z, proj) {
    # 外积计算 u * X_t
    uX <- outer(proj$u, Z, "*") 
    # 特征映射: 利用 cos+sin 的奇偶性抵消交叉项，等价于 L2 距离
    Feat <- cos(uX) + sin(uX) 
    
    # Itô 投影: 内积白噪声 (N x K 矩阵)
    Y <- (1 / sqrt(proj$M)) * (t(Feat) %*% proj$eps)
    return(Y)
}

# ==========================================
# 阶段三 & 四：在线自正则 CUSUM 与 动态边界监测
# ==========================================
online_sn_cusum <- function(Y, m, cv = NULL, gamma = 0.25) {
    n <- nrow(Y)
    T_max <- n - m
    K <- ncol(Y)
    
    # -------------------------------------
    # 历史期 (1 到 m) 免疫洗白矩阵 V_m 构建
    # -------------------------------------
    Y_hist <- Y[1:m, , drop = FALSE]
    Y_hist_bar <- colMeans(Y_hist)
    cumsum_hist <- apply(Y_hist, 2, cumsum)
    
    V_m <- matrix(0, K, K)
    for (j in 1:m) {
        S_j <- cumsum_hist[j, ] - j * Y_hist_bar
        V_m <- V_m + outer(S_j, S_j)
    }
    V_m <- V_m / (m^2)
    # 微小正则项防止极端共线性，确保可逆
    V_m_inv <- solve(V_m + diag(1e-8, K)) 
    
    # -------------------------------------
    # 开放期 (m+1 到 m+T_max) 在线监测
    # -------------------------------------
    Y_mon <- Y[(m+1):n, , drop = FALSE]
    cumsum_mon <- apply(Y_mon, 2, cumsum)
    sum_hist <- cumsum_hist[m, ]
    
    max_test_stat <- 0
    alarm_time <- NA
    
    for (k in 1:T_max) {
        t <- m + k
        tau <- t / m  # 监测比率 (大于 1)
        
        # 多维 CUSUM 向量
        Q_t <- cumsum_mon[k, ] - (k / m) * sum_hist
        
        # 自正则二次型统计量
        T_stat <- (1 / m) * as.numeric(t(Q_t) %*% V_m_inv %*% Q_t)
        
        # Horváth 动态非退化边界函数
        g_tau <- (1 + tau) * (tau / (1 + tau))^gamma
        adjusted_stat <- T_stat / g_tau
        
        if (adjusted_stat > max_test_stat) {
            max_test_stat <- adjusted_stat
        }
        
        # 实时报警判断 (一旦突破阈值，立刻终止)
        if (!is.null(cv) && is.na(alarm_time) && adjusted_stat > cv) {
            alarm_time <- t
            break 
        }
    }
    
    # 返回: 若为求临界值返回最大统计量，若为在线监测返回报警时间
    if (is.null(cv)) return(max_test_stat) else return(alarm_time)
}

# ==========================================
# 开启多线程并行引擎
# ==========================================
num_cores <- parallel::detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat(sprintf(" V3.0 多线程在线引擎启动！调用 %d 核心...\n\n", num_cores))

# 全局参数设定
m_hist <- 200      # 历史无污染期长度
T_max <- 200*10       # 监测期限长度 (总长度 = 1000)
cp_delay <- 150    # 变点发生时间 (m + 150 = 350)
K_dim <- 3         # 随机基底维度
sim_cv <- 2000     # CV 模拟次数
sim_mc <- 1000     # 功效验证次数
rho_val <- 0.2     # 强相依性测试

# 初始化唯一的 V3.0 投影器 (永久锁定白噪声与频率基底)
projector <- init_ito_projector(K = K_dim, M = 500)

start_time <- Sys.time()

# ------------------------------------------
# 步骤 1：模拟多元纯净泛函极限，求解临界值 CV
# ------------------------------------------
cat(sprintf("【步骤 1：利用纯高斯向量模拟 K=%d 的无限期边界分布 (%d 次)】\n", K_dim, sim_cv))
cat("-> 理论之美: 得益于自正则化，无需 AR 数据即可模拟极限阈值！\n")

Tn_null_dist <- foreach(i = 1:sim_cv, .combine = c) %dopar% {
    # 【核心】：直接生成纯 i.i.d 降维后的高斯数据，绕过数据相依性，证明 Pivotality
    Y_pure_gaussian <- matrix(rnorm((m_hist + T_max) * K_dim), ncol = K_dim)
    online_sn_cusum(Y_pure_gaussian, m = m_hist)
}
cv_95 <- quantile(Tn_null_dist, 0.95)
cat(sprintf("=> 模拟完成！95%% 置信水平动态边界乘数 c_alpha = %.3f\n\n", cv_95))

# ------------------------------------------
# 步骤 2：并行执行真实流数据在线监测模拟
# ------------------------------------------
cat(sprintf("【步骤 2：强相依(rho=%.1f)真实流数据测试，对比 H0 与 H1】\n", rho_val))

# 并行运行 H0, H1_mean, H1_var
results <- foreach(i = 1:sim_mc, .combine = rbind) %dopar% {
    
    # 1. 生成三种流数据
    data_H0   <- generate_online_stream(m_hist, T_max, "H0", rho_val, seed_val = i)
    data_H1_M <- generate_online_stream(m_hist, T_max, "H1_mean", rho_val, shift = 1.0, cp_delay = cp_delay, seed_val = i*2)
    data_H1_V <- generate_online_stream(m_hist, T_max, "H1_var", rho_val, var_ratio = 2.5, cp_delay = cp_delay, seed_val = i*3)
    
    # 2. Itô 特征投影映射
    Y_H0   <- apply_ito_projection(data_H0$Z, projector)
    Y_H1_M <- apply_ito_projection(data_H1_M$Z, projector)
    Y_H1_V <- apply_ito_projection(data_H1_V$Z, projector)
    
    # 3. 在线监测 (返回报警时间，NA 为未报警)
    alarm_H0   <- online_sn_cusum(Y_H0, m_hist, cv_95)
    alarm_H1_M <- online_sn_cusum(Y_H1_M, m_hist, cv_95)
    alarm_H1_V <- online_sn_cusum(Y_H1_V, m_hist, cv_95)
    
    c(
        !is.na(alarm_H0), !is.na(alarm_H1_M), !is.na(alarm_H1_V),
        alarm_H1_M - data_H1_M$cp_idx,  # 均值变点检测延迟
        alarm_H1_V - data_H1_V$cp_idx   # 方差变点检测延迟
    )
}

stopCluster(cl)
end_time <- Sys.time()

# ------------------------------------------
# 步骤 3：数据统计与 V3.0 实验报告打印
# ------------------------------------------
size_emp  <- mean(results[, 1])
power_M   <- mean(results[, 2])
power_V   <- mean(results[, 3])

# 计算成功报警条件下的平均检测延迟 (Average Detection Delay, ADD)
add_M <- mean(results[results[, 2] == 1, 4])
add_V <- mean(results[results[, 3] == 1, 5])

cat("\n=======================================================\n")
cat("V3.0 架构在线监测实验报告 (Open-ended SN-Itô-CUSUM)\n")
cat("=======================================================\n")
cat(sprintf("历史窗口 (m)      : %d\n", m_hist))
cat(sprintf("监测期限 (T_max)  : %d\n", T_max))
cat(sprintf("真实变点 (CP)     : %d (m + %d)\n", m_hist + cp_delay, cp_delay))
cat(sprintf("强相依系数 (rho)  : %.2f\n", rho_val))
cat(sprintf("临界边界 (c_alpha): %.3f\n", cv_95))
cat("-------------------------------------------------------\n")
cat(sprintf("Type-I Error (误报率)   : %.3f  (理论目标: 0.050)\n", size_emp))
cat(sprintf("Power (均值变点功效)    : %.3f  (理论目标: 1.000)\n", power_M))
cat(sprintf(" ADD   (均值检测延迟)    : %.1f 个样本\n", add_M))
cat(sprintf("Power (方差变点功效)    : %.3f  (理论目标: 1.000)\n", power_V))
cat(sprintf("ADD   (方差检测延迟)    : %.1f 个样本\n", add_V))
cat("-------------------------------------------------------\n")
cat(sprintf("总耗时: %.2f 秒\n", as.numeric(difftime(end_time, start_time, units="secs"))))
cat("=======================================================\n")

