# Majorana Braiding with ABS: Pauli / d-vector Framework

## 1. 目标

本工作旨在：

* 用 **Pauli 张量基 / d-vector 表示**
* 描述 **Majorana zero modes (MZM)** 的 braiding
* 并复现 **Andreev bound state (ABS)** 对非阿贝尔统计的影响

---

## 2. 有效模型

在固定费米子奇偶子空间中，系统可约化为：

$$
H(u) = d_x(u)\sigma_x + d_y(u)\sigma_y + d_z\sigma_z
$$

其中：

* (d_x, d_y)：Majorana 间耦合（拓扑部分）
* (d_z)：ABS 混合项（破坏拓扑）

---

## 3. Braiding 路径构造

三段路径对应实验中的“开关耦合”过程：

论文步骤:本质就是依次打开/关闭不同耦合，实现交换 $\gamma_2 \leftrightarrow \gamma_3$

step1:初态，只耦合
$i\lambda \gamma_2 \gamma_3$
pairing:(2,3)

step2:关掉2-3,打开
$i\lambda \gamma_1 \gamma_2$
pairing:(1,2)

step3:关1-2,开:
$i\lambda \gamma_1 \gamma_3$
回到原来pairing

### Path 1
2-3-> 1-2: $H_1​(u)=\cos u(iγ_2​γ_3​)+\sin u(iγ_1​γ_2​)$


### Path 2
1-2 -> 1-3
$H_2​(u)=\cos u(iγ_1​γ_2​)+\sin u(iγ_1​γ_3​)$


### Path 3
1-3->2-3
$H_3(u)=\cos u(iγ_1γ_3)+\sin u(iγ_2γ_3)$



其中：

$$
u \in [0, \frac{\pi}{2}]
$$

对应关系:
$$
i\gamma_2 \gamma_3 = \sigma_x\\
i\gamma_1 \gamma_2 = \sigma_y\\
i\gamma_1 \gamma_3 = \sigma_z
$$

正确的 Majorana Braiding 路径（含 ABS 修正）


## 路径定义（核心结果）

### 1️⃣ 理想 braiding 路径（无 ABS）

三段路径分别为：

$$
\begin{aligned}
\vec d_1(u) &= (\cos u,\ \sin u,\ 0) \\
\vec d_2(u) &= (0,\ \cos u,\ \sin u) \\
\vec d_3(u) &= (\sin u,\ 0,\ \cos u)
\end{aligned}
$$

其中：

$$
u \in [0, \tfrac{\pi}{2}]
$$

---

### 2️⃣ 加入 ABS 的修正路径

**关键原则：ABS 是常数偏移，不随路径变化**

$$
\boxed{
\vec d_i(u) = \vec d_i^{\text{(braid)}}(u) + (0,0,d_z)
}
$$

因此三段路径变为：

$$
\begin{aligned}
\vec d_1(u) &= (\cos u,\ \sin u,\ d_z) \\
\vec d_2(u) &= (0,\ \cos u,\ \sin u + d_z) \\
\vec d_3(u) &= (\sin u,\ 0,\ \cos u + d_z)
\end{aligned}
$$

---

## 等价的规范写法（推荐）

为避免混淆，建议写成“基路径 + 常数偏移”：

$$
\vec d_i(u) = \vec d_i^{\text{(braid)}}(u) + (0,0,d_z)
$$

其中：

$$
\begin{aligned}
\vec d_1^{\text{(braid)}}(u) &= (\cos u,\ \sin u,\ 0) \\
\vec d_2^{\text{(braid)}}(u) &= (0,\ \cos u,\ \sin u) \\
\vec d_3^{\text{(braid)}}(u) &= (\sin u,\ 0,\ \cos u)
\end{aligned}
$$

---

## 重要物理约束

* (d_z) 必须为常数
* 不允许出现 (d_z(u))
* 不允许写成 (d_z = \sin u + d_z) 这种随路径变化的形式

---

## 几何解释

* 无 ABS：路径在 Bloch 球表面，包围原点
* 有 ABS：路径整体沿 z 方向平移，不再包围原点

---

## 一句话总结

> 正确的 braiding 路径 = 三段正交平面旋转路径 + 一个常数 (d_z) 偏移（ABS）；ABS 不能随路径变化。



---

## 4. 时间演化算符

每一段路径对应：

$$
U = \mathcal T \exp\left(-i\int H(u),du\right)
$$

总 braiding：

$$
U = U_3 U_2 U_1
$$

---

## 5. 非阿贝尔性判据

构造另一条路径顺序：

$$
U' = U_1 U_3 U_2
$$

定义非对易度：

$$
C = | UU' - U'U |
$$

---

## 6. 物理预期

### (1) 无 ABS（(d_z = 0)）

* 路径绕原点
* 存在简并
* 非阿贝尔统计成立
* (C \neq 0)

---

### (2) 有 ABS（(d_z \neq 0)）

* 路径偏离原点
* gap 打开
* 非阿贝尔性减弱
* (C \to 0)

