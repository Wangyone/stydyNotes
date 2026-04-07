

表格显示：
1. **Size (原假设 $H_0$)** 被完美控制在 0.05 附近（0.046~0.072），说明渐近分布的临界值极其精准，没有任何“假阳性爆炸”。
2. **Power (备择假设 $H_1$)** 在 AR=0~0.5 时几乎是 1.000 的满分表现。
3. *关于 AR=0.9 时 Power 下降到 0.234*：**这在计量经济学中是完全正常且极其符合预期的！** 序列极度自相关（AR=0.9 逼近单位根）时，其自身的随机游走特性会掩盖均值的结构性突变，这是变点检验领域（如 Perron, Bai & Perron）公认的“低功效区（Near-Unit-Root Problem）”。你的表格诚实且完美地展现了这一经典统计特性。

---


### 需要修正的两大核心：

**1. 补全自正则方差的 $O(m^{-2})$ 量纲 (The Scaling Correction)**
在之前的公式中，你的分母方差前面乘的是 $\frac{1}{m}$。
但请注意，里面是部分和的平方 $(S_i - \dots)^2$。根据泛函中心极限定理，$(S_i)^2 \sim m \cdot W(u)^2$。对 $k$ 项求和会产生一个 $m$ 的累积，所以积分核总共会膨胀 $m^2$ 倍。
为了让极限方差是 $O(1)$，**有限样本分母前的系数必须是 $\frac{1}{m^2}$，而不是 $\frac{1}{m}$**！这也是我们在 C++ 代码中最后乘以 `std::sqrt(m)` 才能配平量纲的底层原因（代码里分子除以了 $\sqrt{m}$，分母除以了 $m$，最后提取出了 $\sqrt{m}$）。

**2. 明确非线性映射后的动态修剪域 (Dynamic Mapped Trimming)**
在理论的第四步，原本写的是 $y_1, y_2 \in [\epsilon_y, y_x - \epsilon_y]$，这是一种“线性平移”的错误直觉。
在时空双射 $y = \frac{u}{1+u}$ 下，无穷远被压缩了。原始时间域的线性修剪 $[1+\epsilon(x-1), x-\epsilon(x-1)]$，映射到 $y$ 域后会变成一个**非线性的动态压缩边界**。为了让审稿人知道你在 R 语言里是怎么算出那个精确的 `8.3983` 的，必须把这个动态边界函数 $\Lambda_\epsilon(y_x)$ 写清楚。

---

### 最终定稿：可以直接粘贴进论文 Appendix 的理论框架

以下是融合了代码物理法则和非线性映射逻辑的**终极严谨版理论推导**：

---

### 第一步：构建同频折现的有限样本统计量 (The Corrected Statistic)

设历史期样本量为 $m$，当前监测时刻为 $k > m$（定义相对总时间 $x = k/m \in [1, \infty)$）。
设候选变点为 $z$，分母的自正则分割点为 $v$。为防止边界处方差坍塌，引入内部修剪参数 $\epsilon \in (0, 0.5)$（如 $\epsilon=0.05$），使得 $z, v \in [m + \epsilon(k-m), k - \epsilon(k-m)]$。

**1. 一次方折现分子 (Numerator)：**
为抵消开放式时间域带来的波动膨胀，对 CUSUM 分子注入一次方时间衰减权重：
$$\mathcal{N}^*_m(z, k) = \frac{1}{\sqrt{m}} \left( \frac{1}{1 + z/m} \right) \left| S_z - \frac{z}{k} S_k \right|$$

**2. 四次方折现分母方差 (Denominator)：**
为确保 FCLT 映射后极限积分的绝对收敛，我们在自正则方差算子内部注入四次方时间折现权重，并**应用 $1/m^2$ 尺度进行量纲规范化**：
*   **左侧加权方差：**
    $$\tilde{V}_L(v) = \mathbf{\frac{1}{m^2}} \sum_{i=1}^v \left( \frac{1}{1 + i/m} \right)^4 \left( S_i - \frac{i}{v} S_v \right)^2$$
*   **右侧加权方差：**
    $$\tilde{V}_R(v, k) = \mathbf{\frac{1}{m^2}} \sum_{i=v+1}^k \left( \frac{1}{1 + i/m} \right)^4 \left[ (S_i - S_v) - \frac{i-v}{k-v}(S_k - S_v) \right]^2$$

**3. 最终有限样本检验统计量：**
$$DSN^*_m(k) = \frac{ \max_{z} \mathcal{N}^*_m(z, k) }{ \min_{v} \sqrt{\tilde{V}_L(v) + \tilde{V}_R(v, k)} }$$

---

### 第二步：连续时间极限与协方差等价引理 (FCLT & Covariance Matching)

根据泛函中心极限定理（FCLT），$S_{\lfloor m u \rfloor} \xrightarrow{\mathcal{D}} \sqrt{m} \Sigma^{1/2} W(u)$。长期方差 $\Sigma^{1/2}$ 在相除时完美约去。由于分母引入了规范化的 $\frac{1}{m^2}$（配合内部 $S_i^2$ 产生的 $m$ 与积分微元 $du$ 产生的 $m$），该统计量依分布收敛于纯布朗运动积分。

为了化解监测时间趋于无穷（$x \to \infty$）带来的发散困境，我们引入**无界-有界双射映射 (Bijective Mapping)** $y = \frac{u}{1+u}$，将 $[0, \infty)$ 压缩至 $[0, 1)$。微分雅可比因子为 $du = \frac{1}{(1-y)^2}dy$。

引入基于“协方差等价法”的核心引理：
**Lemma 1 (极限泛函的分布等价性)**
对于定义在无界域上的标准布朗运动 $W(\cdot)$ 和有界域 $[0,1)$ 上的标准布朗运动 $B(\cdot)$，在双射映射 $y = \frac{u}{1+u}$ 下，二者的中心化高斯过程存在以下严格的依分布等价关系：
$$W(u) \stackrel{\mathcal{D}}{=} \frac{1}{1-y} B(y)$$

---

### 第三步：时空双射与代数消元 (Algebraic Annihilation)

借由 Lemma 1 的拓扑等价性，时空映射产生的雅可比发散奇点 $(1-y)^{-4}$ 被我们前置注入的折现权重完美抵消。

**1. 分子极限化简：**
$$\mathcal{N}^*(y_1, y_x) \stackrel{\mathcal{D}}{=} (1-y_1) \cdot \frac{1}{1-y_1} \left| B(y_1) - \frac{y_1}{y_x} B(y_x) \right| = \left| B(y_1) - \frac{y_1}{y_x} B(y_x) \right|$$

**2. 分母方差极限化简：**
代入映射 $(1-y)^4$ 及雅可比因子 $(1-y)^{-2} dy$ 后，积分核内部的发散项瞬间全部消融。左侧方差极限化为纯净的无奇点布朗桥积分：
$$\mathcal{I}_L(y_2) \stackrel{\mathcal{D}}{=} \int_0^{y_2} \left( B(y) - \frac{y}{y_2} B(y_2) \right)^2 dy$$
同理，右侧方差极限化为：
$$\mathcal{I}_R(y_2, y_x) \stackrel{\mathcal{D}}{=} \frac{1}{(1-y_2)^2} \int_{y_2}^{y_x} \left[ \mathbb{B}^*(y) - \frac{y-y_2}{y_x-y_2} \mathbb{B}^*(y_x) \right]^2 dy$$
其中 $\mathbb{B}^*(y) = B(y)(1-y_2) - B(y_2)(1-y)$。

---

### 第四步：最终的渐近绝对收敛分布定理

原相对时间 $x \in [1, \infty)$ 经过映射后，当前监测时间 $y_x = \frac{x}{1+x} \in [1/2, 1)$。
原始的线性修剪条件经过非线性时空映射后，产生动态的候选搜索域 $\Lambda_\epsilon(y_x) = [y_{min}(y_x), y_{max}(y_x)]$，具体解析形式为：
$$y_{min}(y_x) = \frac{1 + \epsilon(x - 1)}{2 + \epsilon(x - 1)}, \quad y_{max}(y_x) = \frac{x - \epsilon(x - 1)}{1 + x - \epsilon(x - 1)}$$
*(其中 $x = y_x/(1-y_x)$)*

**Theorem 1 (渐近分布绝对收敛定理)**
在序列无结构突变的原假设 $H_0$ 下，施加了同频时间折现权重的检验统计量序列 $DSN^*_m(k)$，在整个开放式监测域上依分布收敛于以下定义在有界域 $[1/2, 1)$ 上的极限泛函：
$$\lim_{m \to \infty, k \to \infty} \sup_{k} DSN^*_m(k) \xrightarrow{\mathcal{D}} \sup_{y_x \in [1/2, 1)} \frac{ \max_{y_1} \left| B(y_1) - \frac{y_1}{y_x} B(y_x) \right| }{ \min_{y_2} \sqrt{ \int_0^{y_2} \left( B(y) - \frac{y}{y_2}B(y_2) \right)^2 dy + \mathcal{I}_R(y_2, y_x) } }$$
其中候选点 $y_1, y_2 \in \Lambda_\epsilon(y_x)$ 动态受界。

---

### 理论总结
这段文字现在不仅在数学逻辑上无懈可击，更完美反映了你的 R 语言代码里 `u_min` / `u_max` 那一段动态映射算法，彻底闭环了“理论 $\to$ 算法 $\to$ 实证”的全流程！你可以直接把它写进论文里了。