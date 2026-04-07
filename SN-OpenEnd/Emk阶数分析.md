# $\hat{E}_m(k)$ 的阶数分析（样本均值替换参数估计，正态分布假设）
当参数估计用**样本均值**替换（即 $\hat{\theta}_{a}^{b} = \bar{X}_{a:b} = \frac{1}{b-a+1}\sum_{t=a}^b X_t$），且样本 $X_t \sim i.i.d. N(\mu, \sigma^2)$（正态分布，独立同分布）时，$\hat{E}_m(k)$ 的阶数分析可大幅简化——利用正态分布下样本均值的精确收敛性质（无需渐近近似即可推导阶数），核心结论为：  
- 原假设 $H_0$（无变点）下：$\hat{E}_m(k) = O_p(t)$（与 $m$ 无关，随标准化监测时长 $t$ 线性增长）；  
- 备择假设 $H_1$（有变点）下：$\hat{E}_m(k) = O_p(t\sqrt{m}\Delta)$（主导阶由变点幅度 $\Delta$ 驱动，随 $m$ 和 $t$ 同步增长，保证一致性）。

---

## 一、前置定义与核心前提（简化版，适配样本均值）
### 1. 变量标准化
- $m$：初始训练样本量（固定值，渐近分析中 $m \to \infty$）；  
- $k$：监测期累计观测数，定义 **标准化监测时长** $t = \frac{k}{m}$（开放端场景 $t \to \infty$）；  
- $j$：分割点（$0 \leq j < k$），两段样本量：  
  - 前段（用于估计 $\bar{X}_{1:m+j}$）：$n_1 = m+j = m(1+s)$（$s = \frac{j}{m}$，$0 \leq s < t$）；  
  - 后段（用于估计 $\bar{X}_{m+j+1:m+k}$）：$n_2 = k-j = m(t-s)$。

### 2. 样本均值的核心性质（正态分布下）
因 $X_t \sim i.i.d. N(\mu, \sigma^2)$，样本均值满足：  
- 无偏性：$\mathbb{E}[\bar{X}_{a:b}] = \mu$；  
- 方差精确性：$\text{Var}(\bar{X}_{a:b}) = \frac{\sigma^2}{b-a+1}$，即 $\bar{X}_{a:b} - \mu \sim N\left(0, \frac{\sigma^2}{b-a+1}\right)$；  
- 阶数简化：$\bar{X}_{a:b} - \mu = O_p\left(\frac{1}{\sqrt{b-a+1}}\right)$（与渐近框架一致，但正态分布下为**精确阶数**，非近似）。

### 3. $\hat{E}_m(k)$ 的简化公式（样本均值替换）
论文中 $\hat{E}_m(k)$ 的长期方差矩阵 $\hat{\Sigma}_m$ 对应正态分布的方差 $\sigma^2$（标量，因参数为均值），且 $\hat{\Sigma}_m \xrightarrow{p} \sigma^2$（一致估计），因此简化为：  
$$\hat{E}_m(k) = m^{-1/2} \cdot \max_{j=0}^{k-1} \underbrace{(k-j)}_{\text{后段样本量权重}} \cdot \underbrace{\left| \bar{X}_{1:m+j} - \bar{X}_{m+j+1:m+k} \right|}_{\text{样本均值差异}} \cdot \sigma^{-1}$$  
（注：马氏距离退化为绝对值，$\sigma^{-1}$ 为常数，不影响阶数，后续分析可忽略）。

---

## 二、原假设 $H_0$（无变点，$\mu_1 = \mu_2 = \mu$）下的阶数
$H_0$ 下两段样本均值均为 $\mu$ 的无偏估计，差异仅来自抽样随机波动，核心逻辑：**拆解各部分阶数 → 取最大值的主导阶 → 合并得到 $\hat{E}_m(k)$ 的总阶数**。

### 1. 样本均值差异的阶数
由正态分布样本均值的方差性质：  
$$\left| \bar{X}_{1:m+j} - \bar{X}_{m+j+1:m+k} \right| = O_p\left( \sqrt{\frac{\sigma^2}{n_1} + \frac{\sigma^2}{n_2}} \right) = O_p\left( \sqrt{\frac{1}{m(1+s)} + \frac{1}{m(t-s)}} \right)$$  
（$\sigma^2$ 为常数，阶数由 $\sqrt{\frac{1}{n_1} + \frac{1}{n_2}}$ 决定）。

### 2. 权重因子 $(k-j)$ 的阶数
$k-j = m(t-s)$，因此阶数为 $O(m(t-s))$。

### 3. 单一项 $(k-j) \cdot |\bar{X}_{1:m+j} - \bar{X}_{m+j+1:m+k}|$ 的阶数
将两部分阶数相乘：  
$$O(m(t-s)) \cdot O_p\left( \sqrt{\frac{1}{m(1+s)} + \frac{1}{m(t-s)}} \right) = O_p\left( \sqrt{m} \cdot (t-s) \cdot \sqrt{\frac{1}{1+s} + \frac{1}{t-s}} \right)$$  
化简后：  
$$= O_p\left( \sqrt{m} \cdot \left( \sqrt{\frac{(t-s)^2}{1+s}} + \sqrt{t-s} \right) \right)$$

### 4. 取最大值 $\max_{j=0}^{k-1}$ 的主导阶
遍历所有 $j$（即 $s \in [0,t)$），找到使上式最大的 $s$：  
- 当 $s \to 0$（分割点接近监测起点，对应传统CUSUM的固定基准）：$(t-s) \approx t$，$1+s \approx 1$，因此项的阶数为：  
  $$O_p\left( \sqrt{m} \cdot \left( t + \sqrt{t} \right) \right) = O_p\left( \sqrt{m} \cdot t \right)$$  
- 当 $s \to t$（分割点接近监测终点）：$(t-s) \to 0$，项的阶数为 $O_p(\sqrt{m} \cdot \sqrt{t-s}) \to 0$，无影响；  
- 当 $s = t/2$（分割点在中部）：项的阶数为 $O_p\left( \sqrt{m} \cdot \sqrt{t} \right)$，远小于 $s \to 0$ 时的 $O_p(\sqrt{m}t)$。

因此，$\max_{j=0}^{k-1} \cdot (k-j) \cdot |\bar{X}_{1:m+j} - \bar{X}_{m+j+1:m+k}| = O_p\left( \sqrt{m} \cdot t \right)$。

### 5. $\hat{E}_m(k)$ 的总阶数
代入 $\hat{E}_m(k)$ 的公式，乘以前置因子 $m^{-1/2}$：  
$$\hat{E}_m(k) = m^{-1/2} \cdot O_p\left( \sqrt{m} \cdot t \right) = O_p(t)$$

#### 关键结论（$H_0$ 下）：
$\hat{E}_m(k)$ 的阶数仅与标准化监测时长 $t$ 相关，为 $O_p(t)$，与初始训练样本量 $m$ 无关——这意味着开放端场景下 $t \to \infty$ 时，$\hat{E}_m(k)$ 随 $t$ 线性增长，但结合论文的权重函数 $w(t) = O(1/t)$，加权后统计量 $w(t)\hat{E}_m(k) = O_p(1)$，完美解决了发散问题，可控制误报率。

---

## 三、备择假设 $H_1$（有变点，$\mu_1 \neq \mu_2$）下的阶数
$H_1$ 下参数在 $m+k^*$ 时刻跳变（$\mu_1$ 变 $\mu_2$，跳变幅度 $\Delta = |\mu_1 - \mu_2| > 0$），核心逻辑：**最优分割点主导阶数 → 变点信号放大 → 验证一致性**。

### 1. 最优分割点的选择
当分割点 $j = k^*-1$（恰好落在变点前），两段样本完全分离：  
- 前段 $1:m+j$：全为变点前观测，$\bar{X}_{1:m+j} \xrightarrow{p} \mu_1$，估计误差 $O_p\left( \frac{1}{\sqrt{m(1+t^*)}} \right)$（$t^* = \frac{k^*}{m}$ 为标准化变点时间）；  
- 后段 $m+j+1:m+k$：全为变点后观测，$\bar{X}_{m+j+1:m+k} \xrightarrow{p} \mu_2$，估计误差 $O_p\left( \frac{1}{\sqrt{m(t-t^*)}} \right)$。

因此，样本均值差异的分解式为：  
$$\left| \bar{X}_{1:m+j} - \bar{X}_{m+j+1:m+k} \right| = \Delta + O_p\left( \frac{1}{\sqrt{m(1+t^*)}} + \frac{1}{\sqrt{m(t-t^*)}} \right)$$  
主导阶为 **确定性常数 $\Delta$**（抽样误差可忽略）。

### 2. 单一项的阶数（最优分割点处）
权重因子 $(k-j) = k - (k^*-1) \approx k = mt$（$t \to \infty$ 时），因此：  
$$(k-j) \cdot \left| \bar{X}_{1:m+j} - \bar{X}_{m+j+1:m+k} \right| = O\left( mt \cdot \Delta \right) + O_p\left( \frac{mt}{\sqrt{m}} \right) = O\left( mt\Delta \right) + O_p\left( t\sqrt{m} \right)$$

### 3. $\hat{E}_m(k)$ 的总阶数
乘以前置因子 $m^{-1/2}$：  
$$\hat{E}_m(k) = m^{-1/2} \cdot \left( O\left( mt\Delta \right) + O_p\left( t\sqrt{m} \right) \right) = O\left( t\sqrt{m}\Delta \right) + O_p(t)$$

### 4. 一致性验证（核心性质）
论文要求变点幅度满足 $\sqrt{m}\Delta \to \infty$（$m \to \infty$ 时），因此主导阶为 $O\left( t\sqrt{m}\Delta \right)$，且：  
- 当 $m \to \infty$ 时，$\sqrt{m}\Delta \to \infty$，无论 $t$ 多大（变点多晚），$\hat{E}_m(k)$ 均依概率发散到 $\infty$；  
- 结合权重函数 $w(t) = O(1/t)$，加权后统计量 $w(t)\hat{E}_m(k) = O\left( \sqrt{m}\Delta \right) + O_p(1) \to \infty$，保证以概率1检测到变点，满足一致性。

#### 关键结论（$H_1$ 下）：
$\hat{E}_m(k)$ 的主导阶为 $O\left( t\sqrt{m}\Delta \right)$，随 $m$（初始样本量）和 $t$（监测时长）同步增长，且不受 $t$ 衰减影响——完美解决了传统方法对晚期变点检测功效不足的问题。

---

## 四、总结：阶数分析的核心价值（样本均值替换+正态分布）
| 假设场景 | $\hat{E}_m(k)$ 的阶数 | 与 $m$、$t$ 的关系 | 核心意义 |
|----------|-----------------------|--------------------|----------|
| $H_0$（无变点） | $O_p(t)$ | 与 $m$ 无关，随 $t$ 线性增长 | 结合 $w(t)=O(1/t)$ 可压平噪声，控制误报率 |
| $H_1$（有变点） | $O(t\sqrt{m}\Delta)$ | 随 $m$ 平方根增长、随 $t$ 线性增长 | 变点信号持续放大，保证对晚期变点的检测一致性 |

### 简化后的关键洞察：
样本均值替换+正态分布的设定，让 $\hat{E}_m(k)$ 的阶数更清晰——没有复杂的参数估计余项，仅保留核心的“噪声阶数”和“信号阶数”，直接验证了论文设计的合理性：  
- 替换 $(m+j)$ 为固定 $m$，将 $H_0$ 下的发散阶从原闭端统计量的 $O_p(t/\sqrt{m})$ 降至 $O_p(t)$，适配开放端权重；  
- 保留 $(k-j)$ 权重，让 $H_1$ 下的信号阶随 $t$ 线性增长，同时利用 $\sqrt{m}\Delta$ 的发散性，确保一致性。
