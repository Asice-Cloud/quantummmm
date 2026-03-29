### Quaternion Representation of Quantum Gates



#### content

- 单比特门的四元数表示
- 推广到多比特，四元数表示的一般结论
- 多比特门的表示，复杂门的表示，以及结论



量子门->复矩阵：

- 单比特门：$SU(2)$
- 多比特门：$SU(2^n)$

SU(2) 与 **单位四元数群**是同构的，因此四元数可以提供一种更几何化、更紧凑的表示方式。

形成群

\[
Sp(1) \cong SU(2)
\]



##### 四元数与 Pauli 矩阵

Pauli 矩阵

\[
I =
\begin{pmatrix}
1 & 0 \\
0 & 1
\end{pmatrix}
\]

\[
X =
\begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix}
\]

\[
Y =
\begin{pmatrix}
0 & -i \\
i & 0
\end{pmatrix}
\]

\[
Z =
\begin{pmatrix}
1 & 0 \\
0 & -1
\end{pmatrix}
\]

对应关系
$$
1 \rightarrow I , i\rightarrow X , j\rightarrow Y, k\rightarrow Z 
$$
因此

\[
a + bi + cj + dk
\]

对应矩阵

\[
aI + bX + cY + dZ
\]



### 单比特门的四元数表示

单比特量子门

\[
U \in SU(2)
\]

可以写为

\[
U = aI + bX + cY + dZ
\]

等价于四元数

\[
q = a + bi + cj + dk
\]

其中

\[
a^2 + b^2 + c^2 + d^2 = 1
\]

这表明 $SU(2)\cong Unit Quaternion$





######  旋转门

Bloch 球旋转

\[
U = e^{-i\theta \vec{n}\cdot\vec{\sigma}/2}
\]

对应四元数

\[
q =
\cos(\theta/2)
+
(n_x i + n_y j + n_z k)
\sin(\theta/2)
\]

这是标准的三维旋转四元数。



##### 常见单比特门

###### X 门

\[
X = i
\]

###### Y 门

\[
Y = j
\]

###### Z 门

\[
Z = k
\]

###### Hadamard 门

\[
H = \frac{1}{\sqrt{2}}(X+Z)
\]

四元数

\[
H = \frac{1}{\sqrt{2}}(i + k)
\]

###### 相位门 S

\[
S =
\begin{pmatrix}
1 & 0 \\
0 & i
\end{pmatrix}
\]

四元数形式

\[
S = \cos(\pi/4) + k \sin(\pi/4)
\]



### 多比特系统一般结论

$n $比特 Hilbert 空间

\[
\mathcal H = (\mathbb{C}^2)^{\otimes n}
\]

维度

\[
2^n
\]

量子门

\[
U \in SU(2^n)
\]



单比特门可以用四元数表示，因为

\[
SU(2) \cong Sp(1)
\]

其中 \(Sp(1)\) 是单位四元数群。

然而多比特系统中量子门属于

\[
SU(2^n)
\]

普通四元数已经不足以表示更高维矩阵，因此推广为 **matrix-valued quaternion**,将四元数系数改造为矩阵



#####  matrix-valued quaternion

定义

\[
Q = A + Bi + Cj + Dk
\]

其中

\[
A,B,C,D \in M_{2^{n-1}}(\mathbb C)
\]

这里：

- \(i,j,k\) 是四元数基底
- 系数是复矩阵

这种结构为

\[
M_{2^{n-1}}(\mathbb C) \otimes \mathbb H
\]



######  矩阵表示

利用 Pauli 矩阵与四元数的对应

可以将

\[
Q = A + Bi + Cj + Dk
\]

写为矩阵形式

\[
Q =
I\otimes A
+
X\otimes B
+
Y\otimes C
+
Z\otimes D
\]

因此

\[
Q \in M_{2^n}(\mathbb C)
\]

这表明这种表示等价于Pauli decomposition:
$$
U = \sum_{a_1,...a_k} C_{a_1...a_k}\sigma_{a1}\otimes ... \otimes \sigma_{a_k}
$$


---

#####  物理意义

为了理解该表示的物理含义，我们将 Hilbert 空间分解为

\[
\mathcal H
=
\mathbb C^2
\otimes
\mathbb C^{2^{n-1}}
\]

其中

- 因子1：第一比特
- 因子2：剩余 \(n-1\) 比特



###### 任意算符的分块结构

任意算符

\[
U \in M_{2^n}(\mathbb C)
\]

可以写成 2×2 分块矩阵

\[
U =
\begin{pmatrix}
U_{00} & U_{01} \\
U_{10} & U_{11}
\end{pmatrix}
\]

其中

\[
U_{ij} \in M_{2^{n-1}}(\mathbb C)
\]

###### Pauli 基底展开

Pauli 矩阵构成 2×2 矩阵空间的基底：

\[
\{I,X,Y,Z\}
\]

因此任意算符都可以展开为

\[
U =
I\otimes A
+
X\otimes B
+
Y\otimes C
+
Z\otimes D
\]

其中

\[
A,B,C,D \in M_{2^{n-1}}(\mathbb C)
\]



###### Quaternion 解释

由于 $X->i ,Y->j , Z->k$

因此

\[
U =
A + Bi + Cj + Dk
\]

这说明, 第一个因子(第一比特)的自由度 = 四元数基底, 第二个因子(其余比特) = 系数矩阵

所以matrix-valued quaternion其实是量子系统的分层表示



#### 多比特门示例

##### 1.两比特门示例

两比特 Hilbert 空间：

\[
\mathbb C^4
\]

任意门

\[
U \in SU(4)
\]

可以写为

\[
U =
I\otimes A
+
X\otimes B
+
Y\otimes C
+
Z\otimes D
\]

其中

\[
A,B,C,D \in M_2(\mathbb C)
\]

---

###### 1.1 CNOT 门

CNOT 的矩阵为

\[
CNOT =
\begin{pmatrix}
1&0&0&0 \\
0&1&0&0 \\
0&0&0&1 \\
0&0&1&0
\end{pmatrix}
\]

它可以写成

\[
CNOT =
|0\rangle\langle0|\otimes I
+
|1\rangle\langle1|\otimes X
\]

使用

\[
|0\rangle\langle0| = \frac12(I+Z)
\]

\[
|1\rangle\langle1| = \frac12(I-Z)
\]

得到

\[
CNOT =
\frac12
(I\otimes I
+
Z\otimes I
+
I\otimes X
-
Z\otimes X)
\]

因此
$$
A = \frac{1}{2} (I + X)\\
D = \frac{1}{2} (I - X)\\
B = 0\\
C = 0
$$
所以 quaternion 形式

\[
CNOT = A + kD
\]

###### 1.2 CZ 门

CZ 门矩阵

\[
CZ =
\begin{pmatrix}
1&0&0&0 \\
0&1&0&0 \\
0&0&1&0 \\
0&0&0&-1
\end{pmatrix}
\]

写为

\[
CZ =
|0\rangle\langle0|\otimes I
+
|1\rangle\langle1|\otimes Z
\]

展开得到

\[
CZ =
\frac12
(I\otimes I
+
Z\otimes I
+
I\otimes Z
-
Z\otimes Z)
\]

因此
$$
A = \frac{1}{2} (I + Z)\\
D = \frac{1}{2} (I - Z)\\
B = 0\\
C = 0
$$

###### 1.3 SWAP 门

SWAP 门矩阵

\[
SWAP =
\begin{pmatrix}
1&0&0&0 \\
0&0&1&0 \\
0&1&0&0 \\
0&0&0&1
\end{pmatrix}
\]

Pauli 分解为

\[
SWAP =
\frac12
(I\otimes I
+
X\otimes X
+
Y\otimes Y
+
Z\otimes Z)
\]

因此
$$
A = \frac{1}{2} I\\
B = \frac{1}{2} X\\
C = \frac{1}{2} Y\\
D = \frac{1}{2} Z
$$
对应 quaternion 表示

\[
SWAP =
\frac12
(I + iX + jY + kZ)
\]



###### 1.4 iSWAP 门

iSWAP 门矩阵

\[
iSWAP =
\begin{pmatrix}
1&0&0&0 \\
0&0&i&0 \\
0&i&0&0 \\
0&0&0&1
\end{pmatrix}
\]



\[
\text{iSWAP}
=
\exp
\left(
i\frac{\pi}{4}(X\otimes X + Y\otimes Y)
\right)
\]

利用指数公式

\[
e^{i\theta A}
=
\cos\theta\,I
+
i\sin\theta\,A
\]

得到

\[
\text{iSWAP}
=
\cos(\pi/4)I
+
i\sin(\pi/4)(X\otimes X + Y\otimes Y)
\]

由于

\[
\cos(\pi/4)=\sin(\pi/4)=\frac{1}{\sqrt2}
\]

最终得到

\[
\text{iSWAP}
=
\frac{1}{\sqrt2}I
+
\frac{i}{\sqrt2}(X\otimes X + Y\otimes Y)
\]



现在将其写成

\[
U =
I\otimes A
+
X\otimes B
+
Y\otimes C
+
Z\otimes D
\]

的形式。

\[
A = \frac{1}{\sqrt2}I
\]

\[
B = \frac{i}{\sqrt2}X
\]

\[
C = \frac{i}{\sqrt2}Y
\]

iswap里面没有 \(Z\otimes*\) 项，因此

\[
D = 0
\]

最终 Matrix-valued Quaternion 表达
\[
\text{iSWAP}
=
A + Bi + Cj + Dk
\]





##### 2.三比特门的 表示

三比特 Hilbert 空间为

\[
\mathcal H = (\mathbb C^2)^{\otimes 3}
\]

维度

\[
2^3 = 8
\]

任意三比特门

\[
U \in SU(8)
\]

都可以用 matrix-valued quaternion 形式表示。





###### 2.1 一般形式

将 Hilbert 空间分解为

\[
\mathcal H =
\mathbb C^2 \otimes \mathbb C^4
\]

因此任意三比特算符都可以写成

\[
U =
I\otimes A
+
X\otimes B
+
Y\otimes C
+
Z\otimes D
\]

其中

\[
A,B,C,D \in M_4(\mathbb C)
\]

等价 quaternion 表示

\[
U = A + Bi + Cj + Dk
\]

这里 A,B,C,D都是4*4的复矩阵



---

###### 2.2 Toffoli 门

Toffoli（CCNOT）门矩阵为

\[
|a,b,c\rangle
\rightarrow
|a,b,c \oplus ab\rangle
\]

其算符形式

\[
\text{Toffoli}
=
|0\rangle\langle0|\otimes I_4
+
|1\rangle\langle1|\otimes U
\]

其中

\[
U =
\begin{pmatrix}
1&0&0&0\\
0&1&0&0\\
0&0&0&1\\
0&0&1&0
\end{pmatrix}
\]



Pauli 表示:

使用

\[
|0\rangle\langle0| = \frac12(I+Z)
\]

\[
|1\rangle\langle1| = \frac12(I-Z)
\]

得到

\[
\text{Toffoli}
=
\frac12
(I\otimes(I_4+U))
+
\frac12
(Z\otimes(I_4-U))
\]

Quaternion 表示

因此
$$
A = \frac{1}{2} (I₄ + U)\\
D = \frac{1}{2} (I₄ - U)\\
B = 0\\
C = 0
$$


所以

\[
\text{Toffoli}
=
A + kD
\]





##### 3.复杂例子和意义

通过√iSWAP 和 fSim 门的 Matrix-Valued Quaternion 表示可以看到使用这样的表示的控制作用

###### 3.1 √iSWAP 门

 矩阵形式
\[
\sqrt{\text{iSWAP}} =
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & \frac{1+i}{2} & \frac{1-i}{2} & 0 \\
0 & \frac{1-i}{2} & \frac{1+i}{2} & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
\]

它在 |01⟩ 和 |10⟩ 子空间内产生旋转：

\[
|01\rangle \rightarrow \frac{1+i}{2}|01\rangle + \frac{1-i}{2}|10\rangle
\]

\[
|10\rangle \rightarrow \frac{1-i}{2}|01\rangle + \frac{1+i}{2}|10\rangle
\]



 Pauli/Quaternion 分解
\[
\sqrt{\text{iSWAP}} = \cos(\pi/8) I + i \sin(\pi/8)(X\otimes X + Y\otimes Y)
\]

因此 quaternion 分解：

\[
U = A + Bi + Cj + Dk
\]

- \(A = \cos(\pi/8) I\)  
- \(B = i \sin(\pi/8) X\)  
- \(C = i \sin(\pi/8) Y\)  
- \(D = 0\)

注意到:

- 只包含 \(i,j\) 分量  
- 对应 XY 相互作用平面旋转

---

###### 3.2 fSim 门

\[
\text{fSim}(\theta,\phi) =
\begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & \cos\theta & -i\sin\theta & 0 \\
0 & -i\sin\theta & \cos\theta & 0 \\
0 & 0 & 0 & e^{-i\phi}
\end{pmatrix}
\]

其中 θ 控制 |01⟩ ↔ |10⟩ 的交换，φ 控制 |11⟩ 的相位。



将 Hilbert 空间分为

\[
\mathbb C^2 \otimes \mathbb C^2
\]

得到 quaternion 分解：

\[
U = I \otimes A + X \otimes B + Y \otimes C + Z \otimes D
\]

对应系数：

- \(A = \begin{pmatrix} 1 & 0 \\ 0 & \cos\theta \end{pmatrix}\)  
- \(B = \begin{pmatrix} 0 & 0 \\ 0 & -i\sin\theta \end{pmatrix}\)  
- \(C = \begin{pmatrix} 0 & 0 \\ 0 & -i\sin\theta \end{pmatrix}\)  
- \(D = \begin{pmatrix} 0 & 0 \\ 0 & e^{-i\phi} - \cos\theta \end{pmatrix}\)

检查：

\[
U = A + Bi + Cj + Dk
\]

- XY 子空间旋转 → B, C  
- Z 控制相位 → D  
- I 基底 → A

也有：

- 可以调节 θ 和 φ 来实现不同 entangling 门  
- 适合实现任意 XY 交互 + Z 相位门  
- 用 quaternion 语言自然分离 第一比特旋转 和 剩余比特系数矩阵



总结为：

√iSWAP 和 fSim 门在 matrix-valued quaternion 框架下：

- 都可写作

\[
U = A + Bi + Cj + Dk
\]

- B, C 对应 XY 双线性项  
- D 对应 Z 基底控制  
- A 对应恒等操作分量

因此该表示统一了：

- 两比特 entangling 门  
- XY 平面旋转  
- Majorana braid 的 quaternion 结构  
