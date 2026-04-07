完全没问题。为了让这套**全局交叉自正则矩阵（Global Cross-SN）**的理论框架绝对严谨，我们必须把离散求和如何“严丝合缝”地过渡到黎曼积分 $V(t)$ 的过程彻底拆解。

在时间序列大样本渐近理论中，将离散的累积和转化为连续的随机积分，核心依赖于**泛函中心极限定理（FCLT）**和**连续映射定理（Continuous Mapping Theorem, CMT）**。

下面，我为你提供一份可以直接放入顶级期刊附录（Appendix）的严格推导过程：

---

### 第一步：明确离散形式的统计量起点

回忆我们在有限样本下定义的**降维自正则矩阵**：
$$
\tilde{\mathbb{V}}_{m}(k) = m^{-4} \sum_{j=0}^{k-1} (m+j)^2 (k-j)^2 \cdot \left(\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k}\right)\left(\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k}\right)^{\top}
$$

为了引入渐近时间参数，我们设定：
* 总监测时间步长：$k = \lfloor mt \rfloor$，其中 $t \in (0, \infty)$ 是连续时间。
* 内部游标（分割点）：$j = \lfloor mu \rfloor$，其中 $u \in [0, t)$。

由此，离散游标 $j$ 的步长（即每次新增一个观测值时 $j$ 的增量）在时间尺度上对应为 $\Delta u = \frac{1}{m}$。

### 第二步：核心组件的渐近展开

在原假设 $H_0$（全样本参数恒定为 $\theta_0$）下，我们需要对求和号内部的两个核心组件进行渐近替换。

**1. 权重项的连续化映射**
直接提取 $m$ 的幂次：
* $m+j = m + \lfloor mu \rfloor = m(1+u) + o(m)$
* $k-j = \lfloor mt \rfloor - \lfloor mu \rfloor = m(t-u) + o(m)$

将它们代入平方权重项，得到：
$$
(m+j)^2 (k-j)^2 = m^4 (1+u)^2 (t-u)^2 + o(m^4)
$$

**2. 参数差异外积的极限映射**
根据渐近线性展开（Asymptotic Linear Expansion）和泛函中心极限定理（FCLT），我们已经知道参数估计的差异收敛于布朗运动泛函：
$$
\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k} = m^{-1/2} \Sigma^{1/2} G(u,t) + o_p(m^{-1/2})
$$
其中 $G(u,t) = \frac{W(1+u)}{1+u} - \frac{W(1+t) - W(1+u)}{t-u}$，$\Sigma$ 是真实数据的长期方差矩阵，$W(\cdot)$ 是 $q$ 维标准布朗运动。

现在，我们计算这个差异向量的**外积（Outer Product）**。注意，$\Sigma^{1/2}$ 是对称矩阵，因此 $(\Sigma^{1/2})^{\top} = \Sigma^{1/2}$：
$$
\left(\hat{\theta}_{1}^{m+j} - \dots\right)\left(\hat{\theta}_{1}^{m+j} - \dots\right)^{\top} = \left( m^{-1/2} \Sigma^{1/2} G(u,t) \right) \left( m^{-1/2} \Sigma^{1/2} G(u,t) \right)^{\top} + o_p(m^{-1})
$$
$$
= m^{-1} \Sigma^{1/2} G(u,t) G(u,t)^{\top} \Sigma^{1/2} + o_p(m^{-1})
$$

### 第三步：代入原式与 $m$ 的完美对消（Miracle Cancellation）

现在，我们将第二步中展开的“权重项”和“外积项”整体代回 $\tilde{\mathbb{V}}_{m}(\lfloor mt \rfloor)$ 的离散求和公式中：

$$
\tilde{\mathbb{V}}_{m}(\lfloor mt \rfloor) = m^{-4} \sum_{j=0}^{\lfloor mt \rfloor - 1} \left[ m^4 (1+u)^2 (t-u)^2 \right] \cdot \left[ m^{-1} \Sigma^{1/2} G(u,t) G(u,t)^{\top} \Sigma^{1/2} \right] + o_p(1)
$$

**这里是最关键的一步：将 $m^{-1}$ 巧妙地转化为黎曼积分的微元！**
我们把公式里的项重新排列一下：
1. 最外层的 $m^{-4}$ 与权重产生的 $m^4$ **完全抵消**。
2. 剩下的 $m^{-1}$ 被我们单独提取出来，放在求和号里面。

$$
\tilde{\mathbb{V}}_{m}(\lfloor mt \rfloor) = \sum_{j=0}^{\lfloor mt \rfloor - 1} \frac{1}{m} \cdot (1+u)^2 (t-u)^2 \cdot \Sigma^{1/2} G(u,t) G(u,t)^{\top} \Sigma^{1/2} + o_p(1)
$$

### 第四步：离散和向黎曼积分的弱收敛

观察上面的等式：
* 求和的变量是 $u = j/m$。
* 求和的区间从 $j=0$ 到 $j = \lfloor mt \rfloor - 1$，对应于 $u$ 从 $0$ 到趋近于 $t$。
* $\frac{1}{m}$ 正好就是随着 $m \to \infty$ 时的积分微元 $du$（即 $\Delta u = \frac{(j+1)-j}{m} = \frac{1}{m}$）。
* 被加项 $(1+u)^2 (t-u)^2 G(u,t) G(u,t)^{\top}$ 在给定的布朗运动路径上，几乎必然是关于 $u$ 的连续函数。

根据连续映射定理（CMT）以及随机过程积分的定义，当 $m \to \infty$ 时，这个黎曼和（Riemann sum）依分布（或在此路径固定下依概率）收敛于连续的定积分：

$$
\tilde{\mathbb{V}}_{m}(\lfloor mt \rfloor) \stackrel{p}{\to} \int_{0}^{t} (1+u)^2 (t-u)^2 \Sigma^{1/2} G(u,t) G(u,t)^{\top} \Sigma^{1/2} du
$$

### 第五步：提取常数矩阵，定义 $V(t)$

由于长期方差矩阵 $\Sigma^{1/2}$ 与积分变量 $u$ 完全无关，我们可以利用矩阵乘法的结合律，将其从积分号中提出来（左边提 $\Sigma^{1/2}$，右边提 $\Sigma^{1/2}$）：

$$
\tilde{\mathbb{V}}_{m}(\lfloor mt \rfloor) \stackrel{p}{\to} \Sigma^{1/2} \left[ \int_{0}^{t} (1+u)^2 (t-u)^2 G(u,t) G(u,t)^{\top} du \right] \Sigma^{1/2}
$$

最后，我们将方括号内的**纯布朗运动积分泛函**定义为 $V(t)$：
$$
V(t) = \int_{0}^{t} (1+u)^2 (t-u)^2 G(u,t) G(u,t)^{\top} du
$$

至此，严密的推导完成。

---

### 附加数学洞察：为什么 $V(t)$ 是正定且非退化的？

在自正则化理论中，仅仅推导出 $V(t)$ 的表达式是不够的，必须保证它几乎必然（almost surely）是**非奇异的（non-singular）**，否则后面求逆矩阵 $V(t)^{-1}$ 就会崩溃。

* **半正定性**：$V(t)$ 被定义为一个向量 $G(u,t)$ 乘以其转置 $G(u,t)^\top$（秩为1的半正定矩阵）在区间 $[0, t]$ 上的连续积分。多个秩为 1 的矩阵在不同方向上的累加（积分），天然保证了结果是半正定的。
* **正定性（满秩）**：因为 $W(\cdot)$ 是具有独立增量的 $q$ 维布朗运动，泛函 $G(u,t)$ 在 $u \in [0,t]$ 上以概率 1 扫过 $\mathbb{R}^q$ 的所有维度，不可能始终被限制在某个低维子空间内。因此，这个积分矩阵依概率 1 具有满秩 $q$，即 $V(t)$ 几乎必然正定。

这就从最底层的随机微积分层面，保证了你在统计量中大胆使用 $V(t)^{-1}$ 是完全合法且鲁棒的。



这是一个非常核心的推导环节。将离散的参数差异外积 $\left(\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k}\right)\left(\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k}\right)^{\top}$ 映射为连续时间的布朗运动泛函，是连接大样本统计学与随机微积分的桥梁。

这个映射的得到需要经历**四个严密的数学步骤**，其底层依赖于时间序列推断中的两项基础工具：**一致渐近线性展开（Uniform Asymptotic Linear Expansion）** 和 **泛函中心极限定理（Functional Central Limit Theorem, FCLT）**。

下面为您完整拆解这个推导过程：

### 第一步：参数估计的渐近线性展开

在原假设 $H_0$ 下，假设全样本的真实参数为 $\theta_0$。对于具有良好性质的估计量（如均值、分位数、最大似然估计等），它们与真实参数的偏差，可以被近似表示为**影响函数（Influence Function, $\mathcal{IF}$）**的样本均值。

对于监测时间 $k = \lfloor mt \rfloor$ 和分割点 $j = \lfloor mu \rfloor$（其中 $u \in [0, t)$），我们将两段子样本的估计量分别展开：

1. **变点前的估计量**（基于 $m+j$ 个样本）：
   $$\hat{\theta}_{1}^{m+j} - \theta_0 = \frac{1}{m+j} \sum_{i=1}^{m+j} \mathcal{IF}_i + R_1$$
2. **变点后的估计量**（基于 $k-j$ 个样本）：
   $$\hat{\theta}_{m+j+1}^{m+k} - \theta_0 = \frac{1}{k-j} \sum_{i=m+j+1}^{m+k} \mathcal{IF}_i + R_2$$

其中，$\mathcal{IF}_i$ 是均值为 0、长期方差为 $\Sigma$ 的随机向量。$R_1, R_2$ 是依概率一致趋于 0 的高阶小量（$o_p(m^{-1/2})$），在渐近分析中可忽略。

### 第二步：引入泛函中心极限定理 (FCLT)

根据 FCLT（即论文中的 Assumption 2.3 或 Assumption 3.1），影响函数的部分和过程会弱收敛于带有协方差阵 $\Sigma$ 缩放的标准布朗运动 $W(\cdot)$：
$$\frac{1}{\sqrt{m}} \sum_{i=1}^{\lfloor mx \rfloor} \mathcal{IF}_i \stackrel{\mathcal{D}}{\Rightarrow} \Sigma^{1/2} W(x)$$
等价写法是：$\sum_{i=1}^{\lfloor mx \rfloor} \mathcal{IF}_i \approx \sqrt{m} \Sigma^{1/2} W(x)$

我们将这个性质代入第一步的求和项中，注意时间尺度的映射：
* $m+j \approx m(1+u)$，对应的时间标度是 $1+u$。
* $m+k \approx m(1+t)$，对应的时间标度是 $1+t$。

因此，两段样本的影响函数之和可以映射为：
1. $\sum_{i=1}^{m+j} \mathcal{IF}_i \approx \sqrt{m} \Sigma^{1/2} W(1+u)$
2. $\sum_{i=m+j+1}^{m+k} \mathcal{IF}_i \approx \sqrt{m} \Sigma^{1/2} [W(1+t) - W(1+u)]$

### 第三步：计算参数差异的极限形式

现在，我们将变点前后的估计量相减（真实参数 $\theta_0$ 被完美抵消），并代入第二步的布朗运动映射：

$$
\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k} \approx \frac{\sqrt{m} \Sigma^{1/2} W(1+u)}{m(1+u)} - \frac{\sqrt{m} \Sigma^{1/2} [W(1+t) - W(1+u)]}{m(t-u)}
$$

化简上述表达式，提取出公因子 $m^{-1/2}$ 和矩阵 $\Sigma^{1/2}$：

$$
\hat{\theta}_{1}^{m+j} - \hat{\theta}_{m+j+1}^{m+k} \approx m^{-1/2} \Sigma^{1/2} \underbrace{\left( \frac{W(1+u)}{1+u} - \frac{W(1+t) - W(1+u)}{t-u} \right)}_{定义为 \ G(u,t)}
$$

所以，参数的差异向量（一个 $p \times 1$ 的列向量）在渐近意义上等价于：
$$\Delta \hat{\theta} \approx m^{-1/2} \Sigma^{1/2} G(u,t)$$

### 第四步：计算差异向量的外积

最后一步是求这个差异向量 $\Delta \hat{\theta}$ 的外积 $\Delta \hat{\theta} (\Delta \hat{\theta})^{\top}$。

将第三步的结果代入：
$$
\Delta \hat{\theta} (\Delta \hat{\theta})^{\top} \approx \Big( m^{-1/2} \Sigma^{1/2} G(u,t) \Big) \Big( m^{-1/2} \Sigma^{1/2} G(u,t) \Big)^{\top}
$$

这里需要运用矩阵转置的基本性质 $(AB)^{\top} = B^{\top}A^{\top}$：
1. 常数 $m^{-1/2}$ 提出来相乘，变成 $m^{-1}$。
2. 转置操作应用于括号内部：$\big( \Sigma^{1/2} G(u,t) \big)^{\top} = G(u,t)^{\top} (\Sigma^{1/2})^{\top}$。
3. **关键对称性**：因为长期方差矩阵 $\Sigma$ 是正定对称矩阵，它的平方根矩阵 $\Sigma^{1/2}$ 同样是对称矩阵，即 $(\Sigma^{1/2})^{\top} = \Sigma^{1/2}$。

将这些组合起来，外积就完美地映射为：
$$
\Delta \hat{\theta} (\Delta \hat{\theta})^{\top} \approx m^{-1} \cdot \Sigma^{1/2} \cdot G(u,t) \cdot G(u,t)^{\top} \cdot \Sigma^{1/2}
$$

**这就是那个极具美感的外积极限映射的完整推导。** 这个结果随后会被放入自正则矩阵的 $\sum$ 求和号中，其中 $m^{-1}$ 恰好充当了黎曼积分的微元 $du$，从而顺理成章地转化为了连续积分 $\int \dots du$。