# 导入所需的包
library(ggplot2)
library(foreach)
library(doParallel)

# ==========================================
# 核心函数定义区 (保持不变，确保并行节点均可调用)
# ==========================================

rho <- 0.0
# 1：数据生成 (AR(1) 模型)
generate_data <- function(n, model = "H0", shift_size = 2, seed_val = NULL) {
    if (!is.null(seed_val)) set.seed(seed_val)
    Z <- numeric(n)
    err <- rnorm(n, 0, 1)
    Z[1] <- err[1]
    for (t in 2:n) {
        Z[t] <- rho * Z[t-1] + err[t]
    }
    if (model == "H1") {
        cp_idx <- floor(n / 2)
        # Z[(cp_idx + 1):n] <- Z[(cp_idx + 1):n] + shift_size
        Z[(cp_idx + 1):n] <- Z[(cp_idx + 1):n]*(1+shift_size)
    }
    return(Z)
}

# 2：特征函数 (ECF) 的 K 维随机投影
project_ecf <- function(Z, K = 5, U = 5, M = 500, proj_seed = 999) {
    n <- length(Z)
    u_grid <- seq(-U, U, length.out = M)
    du <- 2 * U / M
    
    # 保证所有的并行核使用完全一致的随机投影方向
    set.seed(proj_seed)
    V_rand <- matrix(0, nrow = M, ncol = K)
    for (k in 1:K) {
        V_rand[, k] <- cumsum(rnorm(M)) * sqrt(du)
    }
    
    P_matrix <- matrix(0, nrow = n, ncol = K)
    for (t in 1:n) {
        X_t_u <- cos(u_grid * Z[t]) + sin(u_grid * Z[t])
        P_matrix[t, ] <- as.vector(t(X_t_u) %*% V_rand) * du
    }
    return(P_matrix)
}

# 3：严格标准化的多元自正则变点检验
sn_cusum_test <- function(P, trim = 0.15) {
    n <- nrow(P)
    K <- ncol(P)
    
    cumsum_P <- apply(P, 2, cumsum)
    P_bar <- cumsum_P[n, ] / n
    
    start_idx <- floor(n * trim)
    end_idx <- floor(n * (1 - trim))
    
    max_Tn <- 0
    
    for (k in start_idx:end_idx) {
        S_k <- cumsum_P[k, ] - k * P_bar
        mean_L <- cumsum_P[k, ] / k
        mean_R <- (cumsum_P[n, ] - cumsum_P[k, ]) / (n - k)
        
        L_matrix <- cumsum_P[1:k, , drop=FALSE] - matrix(1:k, ncol=1) %*% matrix(mean_L, nrow=1)
        R_cumsum <- sweep(cumsum_P[(k+1):n, , drop=FALSE], 2, cumsum_P[k, ], "-")
        R_matrix <- R_cumsum - matrix(1:(n-k), ncol=1) %*% matrix(mean_R, nrow=1)
        
        V_k <- (t(L_matrix) %*% L_matrix + t(R_matrix) %*% R_matrix) / (n^2)
        V_k <- V_k + diag(1e-8, K) 
        
        Tn_k <- (1/n) * sum(S_k * solve(V_k, S_k))
        if (Tn_k > max_Tn) max_Tn <- Tn_k
    }
    return(max_Tn)
}

# ==========================================
# 开启多线程引擎
# ==========================================
# 自动检测 CPU 核心数，保留 1 个核心给操作系统避免卡顿
num_cores <- parallel::detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat(sprintf("多线程引擎已启动！调用了 %d 个 CPU 核心并行计算...\n\n", num_cores))

# 实验参数设定
K_dim <- 5
n_samples <- 400
sim_cv <- 2000    
sim_mc <- 1000    

start_time_total <- Sys.time()

# ------------------------------------------
# 阶段一：并行模拟渐近枢轴分布临界值
# ------------------------------------------
cat(sprintf("【阶段一：并行模拟 K=%d 维多元自正则渐近分布临界值 (%d 次)】\n", K_dim, sim_cv))

# 使用 %dopar% 将 for 循环分发到各个核心，.combine='c' 会自动将结果合并为向量
Tn_null_dist <- foreach(i = 1:sim_cv, .combine = c) %dopar% {
    # 节点内部独立生成正态数据
    P_iid <- matrix(rnorm(n_samples * K_dim), nrow = n_samples, ncol = K_dim)
    sn_cusum_test(P_iid)
}

cv_95 <- quantile(Tn_null_dist, 0.95)
cat(sprintf("=> 模拟完成！95%% 置信水平临界值 CV = %.3f\n\n", cv_95))

# ------------------------------------------
# 阶段二：并行执行蒙特卡洛拒绝率模拟
# ------------------------------------------
cat(sprintf("【阶段二：并行执行 %d 次蒙特卡洛拒绝率模拟 (n=%d)】\n", sim_mc, n_samples))

# .combine='rbind' 会自动将每次循环返回的 c(rej_H0, rej_H1) 拼成一个大矩阵
results <- foreach(i = 1:sim_mc, .combine = rbind, .packages = c()) %dopar% {
    # 设置动态的随机数种子以保证每组数据独立且可复现
    Z_H0 <- generate_data(n_samples, model = "H0", seed_val = i * 10)
    Z_H1 <- generate_data(n_samples, model = "H1", shift_size = 0.8, seed_val = i * 10 + 1)
    
    Tn_H0 <- sn_cusum_test(project_ecf(Z_H0, K = K_dim))
    Tn_H1 <- sn_cusum_test(project_ecf(Z_H1, K = K_dim))
    
    # 返回布尔值：1 表示拒绝，0 表示不拒绝
    c(ifelse(Tn_H0 > cv_95, 1, 0), ifelse(Tn_H1 > cv_95, 1, 0))
}

# 汇总并行计算的结果
reject_H0 <- sum(results[, 1])
reject_H1 <- sum(results[, 2])

size_emp <- reject_H0 / sim_mc
power_emp <- reject_H1 / sim_mc

# 关闭集群释放资源
stopCluster(cl)
end_time_total <- Sys.time()

# ------------------------------------------
# 打印最终报告
# ------------------------------------------
cat("\n==========================================\n")
cat("多线程蒙特卡洛实验报告 (Parallel MC Results)\n")
cat("==========================================\n")
cat(sprintf("调用核心数 (CPU Cores)  : %d\n", num_cores))
cat(sprintf("重复次数 (Replications) : %d\n", sim_mc))
cat(sprintf("样本容量 (Sample Size)  : %d\n", n_samples))
cat(sprintf("投影维度 (Random Proj)  : K = %d\n", K_dim))
cat(sprintf("模拟临界值 (CV at 0.95) : %.3f\n", cv_95))
cat("---------\n")


cat(sprintf("AR(1)系数 : %.3f \n", rho))
cat(sprintf("第一类错误率 (Empirical Size) : %.3f  (理论目标: 0.05)\n", size_emp))
cat(sprintf("检验功效 (Empirical Power)    : %.3f  (理论目标: 1.00)\n", power_emp))
cat(sprintf("多线程总耗时 (Total Time)     : %.2f 秒\n", as.numeric(difftime(end_time_total, start_time_total, units="secs"))))
cat("==========================================\n")

