# 清理环境，重新加载
rm(list = ls())
library(doParallel)
library(foreach)

M <- 10000        # 万次模拟，确保临界值精准
N <- 2000         
epsilon <- 0.05   
tau_min <- 0.5    
tau_max <- 0.99   

cl <- makeCluster(max(1, detectCores() - 1))
registerDoParallel(cl)

dt <- 1 / N
grid_tau <- seq(0, 1, length.out = N + 1)
idx_monitor <- which(grid_tau >= tau_min & grid_tau <= tau_max)

cat("正在模拟绝对纯净的 Global SN 极限泛函临界值...\n")
start_time <- Sys.time()

results_cv <- foreach(i = 1:M, .combine = c) %dopar% {
    # 1. 生成布朗运动
    dW <- rnorm(N, mean = 0, sd = sqrt(dt))
    B <- c(0, cumsum(dW))
    
    # 2. 预计算分母积分项的前缀和
    P_B2 <- cumsum(B^2) * dt
    P_By <- cumsum(B * grid_tau) * dt
    P_y2 <- (grid_tau^3) / 3   # y^2 的精确积分解析解
    
    sup_stat <- 0 
    
    # 3. 遍历监测域
    for (t_idx in idx_monitor) {
        tau <- grid_tau[t_idx]
        B_tau <- B[t_idx]
        c_val <- B_tau / tau  # 布朗桥的斜率 (y/y_x)*B(y_x)
        
        # 全局方差分母：一步到位，没有任何发散项！
        denom_var <- P_B2[t_idx] - 2 * c_val * P_By[t_idx] + (c_val^2) * P_y2[t_idx]
        if (denom_var <= 1e-10) next
        denom <- sqrt(denom_var)
        
        idx_v <- which(grid_tau >= epsilon & grid_tau <= (tau - epsilon))
        if (length(idx_v) == 0) next
        
        # 修复的分子：严格的纯布朗桥距离
        numerators <- abs(B[idx_v] - c_val * grid_tau[idx_v])
        
        stat <- max(numerators) / denom
        if (stat > sup_stat) sup_stat <- stat
    }
    return(sup_stat)
}
stopCluster(cl)

# 打印真实的临界值 (应该在 1.5 到 2.5 之间)
cv_95_true <- quantile(results_cv, 0.95)
cat(">>> 真实的 95% 临界值 (cv_95) =", round(cv_95_true, 4), "<<<\n")


# 绘图：直方图 + 密度曲线
hist(results_cv, 
     prob = TRUE,        # 关键：把直方图转为概率密度（才能和密度图叠加）
     col = "lightblue",  # 直方图颜色
     border = "white",   # 边框白色
     main = "结果分布：直方图 + 密度曲线", 
     xlab = "数值")

# 叠加密度曲线
lines(density(results_cv), 
      col = "red", 
      lwd = 3)           # 线条粗细

# 可选：加一条均值线
abline(v = mean(results), col = "blue", lwd = 2, lty = 2)