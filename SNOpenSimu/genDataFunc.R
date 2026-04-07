# rm(list = ls())
# setwd("D:\\Desktop\\SNOpenSimu")

genData <- function(m, Tm, model, cp_frac = 0.5, delta = 1) {
  # m: 历史期样本量 (无变点)
  # Tm: 监测期样本量 (变点发生区间)
  # model: 模型编号 (1-5为H0, 101-105为H1)
  # cp_frac: 变点在监测期内的相对位置 (例如0.25, 0.5)
  # 变点发生在 $m + cp_frac * Tm$ 处
  # delta: 变点发生后的均值偏移量
  
  # 设置烧链期，用于消除初始值对时间序列的影响
  burnin <- ceiling(m / 2)
  obs <- burnin + m + Tm
  
  # 初始化序列
  x <- numeric(obs)
  
  # ==========================================
  # 生成基础序列（涵盖 H0 和 H1 的 5 种时间依赖性）
  # ==========================================
  if (model == 1 || model == 101) {          # beta = 0 (白噪声)
    x <- rnorm(obs, mean = 0, sd = 1)
  } 
  else if (model == 2 || model == 102) {     # beta = 0.3
    x <- arima.sim(model = list(ar = 0.3), n = obs)
  } 
  else if (model == 3 || model == 103) {     # beta = 0.5
    x <- arima.sim(model = list(ar = 0.5), n = obs)
  } 
  else if (model == 4 || model == 104) {     # beta = 0.7
    x <- arima.sim(model = list(ar = 0.7), n = obs)
  } 
  else if (model == 5 || model == 105) {     # beta = 0.9
    x <- arima.sim(model = list(ar = 0.9), n = obs)
  }
  else {
    stop("未定义的模型编号 (Undefined model index)")
  }
  
  # ==========================================
  # H1 备择假设：在监测期内施加均值变点
  # ==========================================
  if (model %in% 101:105) {
    # 变点绝对位置：烧链期 + 历史期 + 监测期内的相对位置
    change.pt <- burnin + m + floor(Tm * cp_frac)
    
    # 在变点及之后的数据上施加均值偏移 delta
    x[change.pt:obs] <- x[change.pt:obs] + delta
  }
  
  # 截去烧链期的数据，仅返回长度为 m + Tm 的目标时间序列
  x_final <- x[(burnin + 1):obs]
  
  return(x_final)
}