rm(list = ls())

library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(foreach)

# 设置工作路径 (请根据实际情况修改)
setwd("D:/Desktop/SNOpenSimu")
source("D:/Desktop/SNOpenSimu/genDataFunc.R", echo=TRUE)
# ==========================================
# 1. 核心参数设定
# ==========================================
T_0 <- 300                 # 历史期样本量 (m)
L0 <- 10                   # 监测期倍数
Tm <- L0 * T_0             # 监测期样本量 (nT = 1500)
N <- T_0 + Tm              # 总样本量 (1800)
nsim <- 1000               # 模拟次数

# SNEmk 检验参数
final_quantile <- 8.6833 # sqrt(1861.676)   # 给定的临界值
weights <- rep(1, Tm)      # 权重函数 w(t) = 1 (长度必须为 nT)

# 变点参数 (用于 H1)
cp_frac <- 0.1             # 变点位置
delta <- 1                 # 变点均值偏移量

# 待模拟的模型编号: 1-5 (H0/Type I error), 101-105 (H1/Power)
models <- c(1:5, 101:105)



# ==========================================
# 3. 开启并行环境
# ==========================================
# 留出 1 个核心以防系统卡顿
cores_to_use <- max(1, detectCores() - 1)
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

cat("开始并行模拟，启用的核心数:", cores_to_use, "\n")
startTime <- Sys.time()

# ==========================================
# 4. 并行执行模拟实验
# ==========================================
results <- foreach(mod = models, .combine = rbind, .packages = c("Rcpp", "RcppArmadillo")) %dopar% {
  
  # 关键步骤：Windows环境下，并行子进程是独立且空白的R环境。
  # 必须在每个子进程内部重新 sourceCpp 编译/加载 C++ 函数。
  # 请确保子进程能找到该 .cpp 文件，建议使用绝对路径。
  Rcpp::sourceCpp("D:/Desktop/SNOpenSimu/statSN.cpp")
  
  rejections <- 0
  
  for (s in 1:nsim) {
    # Step 1: 生成数据
    yt <- genData(m = T_0, Tm = Tm, model = mod, cp_frac = cp_frac, delta = delta)
    
    # Step 2: 计算 SNEmk 统计量序列
    sn_stats <- compute_stat(yt, T_0, "snemk", weights)
    
    # Step 3: 提取监控期内的最大值，并与临界值比较
    # 注意：自正则统计量在刚开始监控的前几个样本处可能有极端值（矩阵接近奇异）
    # 在实际工程中建议切除前几个点，此处保留严格的全局最大值判断
    burn_in <- floor(Tm * 0.01) 
    
    # Step 3: Split-Sample SN 不存在启动期发散，直接取全局最大值！
    if (max(sn_stats) > final_quantile) {
      rejections <- rejections + 1
    }
  }
  
  # 返回当前模型的结果
  rate <- rejections / nsim
  data.frame(Model = mod, Rejection_Rate = rate)
}

# 关闭并行集群
stopCluster(cl)

endTime <- Sys.time()
cat("模拟完成！总耗时:", round(difftime(endTime, startTime, units = "mins"), 2), "分钟\n\n")

# ==========================================
# 5. 整理与输出结果
# ==========================================
results$Type <- ifelse(results$Model <= 5, "Type I Error (H0)", "Power (H1)")
results$AR_Coef <- rep(c(0, 0.3, 0.5, 0.7, 0.9), 2)

cat("============== 模拟结果 ==============\n")
print(results[, c("Type", "Model", "AR_Coef", "Rejection_Rate")], row.names = FALSE)