## 拓扑量子计算的理论学习



**主要学习内容** 

1.Clifford代数的基本入门和辫子群等

2.Anyous的应用



参考文献：

``` 
近世代数 ---丘维声
AN INTRODUCTION TO CLIFFORD ALGEBRAS AND SPINORS --Jayme Vaz, Jr ; Roldão da Rocha, Jr
微分几何引论 --陈维桓
Topological Quantum Computation ---Zhenghan Wang
Non-Abelian Anyons and Topological Quantum Computation
```



#### 1.Clifford代数和辨群

Clifford代数内容很宽泛，建立了一套完整的几何变换的描述。比如上面提到的旋转，诸如$R^2,R^3$上的旋转都可写成更大代数的乘积，比如$R^3$中可以使用单位四元数来描述旋转,$h\longmapsto qh\bar q$。而物理中最常见的对称化和反对称化运算，在Clifford中也可以由内积和外积来描述:

$Def: (V,<\cdot,\cdot>)$是一个$n$维内积空间，接着定义Clifford代数：
$$
CL(V):=\frac{TV}{(x\otimes x+<x,x>\cdot 1|x\in V)}
$$
TV为V的张量代数

即Clifford满足**生成关系**：$x\otimes x=-<x,x>\cdot 1=-|x|^2\cdot1$ 



物理上的联系：

如：Dirac 方程 $(iγ^μ ∂_μ - m)ψ = 0$，γ 矩阵满足 ${γ^μ,γ^ν}=2η^{μν}$

以及**Majorana fermion**:

首先是费米子关系式：

$ f =  \tfrac{1}{\sqrt{2}} (\gamma_1+i\gamma_2),$
$ f^{\dagger} =  \tfrac{1}{\sqrt{2}} (\gamma_1-i\gamma_2) ~.$

而Majorana费米子满足: $\{γ_i, γ_j\} = γ_i γ_j + γ_j γ_i = 2 δ_{ij}, γ_i† = γ_i.$

这是实 Clifford 代数，由一组满足上述反对易关系的生成元构成。所以从代数角度看，Majorana 系统就是在操作一个 Clifford 代数的表示

从 Majorana 到复费米子与 Hilbert 空间维度

- 把两个 Majorana 配对成一个常规复杂费米子：
    $c = (γ_1 + i γ_2)/2, c† = (γ_1 - i γ_2)/2.$
- 若有 $2n$ 个 Majorana，则可构造 n 个复费米子，整个费米子希尔伯特空间大小为 $2^n$。通过张量积乘在一起连接起来。



**变换关系：辫子群braid group**

设一组Majorana算子$\{c_n\}$

构造辫子算子：$\sigma_k=(1+c_{k+1}c_k)\sqrt2, $   $\sigma_k^{-1}=(1-c_{k+1}c_k)\sqrt2 $    ,$k=1,2...n$,定义$c_{n+1}=c_1$

则有$σ_i σ_j = σ_j σ_i$, 若$ |i-j| \ge 2$，  $σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}$   ------> 辨群$B_n$

辫子操作： $T_k(x)=\sigma_k x \sigma_k^{-1}$ , $T_k(c_k)=c_{k+1}, T_k(c_{l+1})=-c_k ;T_k(c_i)=c_i ,\space if\space i \ne k, k+1 $

![img](https://picx.zhimg.com/v2-ae7fa28d440a51f7449f54241e769069_r.jpg?source=1def8aca)



**推广**： 对于一组Clifford代数里的算子$M_i$，都可以用$\sigma_i=(1+M_i)/\sqrt2$来构造Artin辫子群，

比如用四元数 $A=\frac{1}{\sqrt2}(1+i),...$  , 现在就在做四元数等更多的构造角度来看待Majorana.



#### 2.Anyous 和拓扑量子计算

个人理解： 拓扑量子计算利用二维的任意子，利用辨群结构充当逻辑门来进行计算。（图片来自wiki）

![undefined](https://upload.wikimedia.org/wikipedia/commons/0/0f/Topological_quantum_computer.jpg)

**The unitary operation corresponding to exchanging anyons depends only on the topology of the braid**



**Anyons**：

任意子是对费米子和玻色子概念的广义化，其交换可以产生任意相位或更一般的幺正变换。



Abel anyons：$\mid\psi_1\psi_2>=e^{i\theta}\mid\psi_2\psi_1>$

Non abel anyons： 不同交换顺序不一定可交换;  

可以根据融合代数：

 $a×b=∑_cN^a_{bc} c（N^a_{bc}为非负整数）$

辫子表示：给定 $n$个任何子，辫子群 $B_n$在结合空间 $H$上有表示 $ρ:Bn→U(H)$；每个生成元 $σi$（交换第 i 与 i+1粒子）映为幺正矩阵 $ρ(σi)$。

计算 $σi$的常用步骤：若在当前 fusion tree 下直接对角，则 $σi=diag⁡(R)$；否则用 F变基再施加 diag⁡(R)，再变回原基：$σi=F−1 diag⁡(R)$ 。 其中：

- F是多粒子结合顺序变换的矩阵,即基变换：决定不同融合顺序下态的展开系数
- R给出交换两个任何子的局部相位或矩阵元：决定把两个任何子绕动一圈的局部效果

例如Fibonacci anyon: {1,τ}，$τ×τ=1+τ$，其辫子表示在逻辑子空间上是稠密的，单靠编辫子即可近似任意单比特幺正（即拓扑上普适）。在结合关系基础上构建辫子群，形成门电路...

![image-20260124005830555](/home/asice-cloud/.config/Typora/typora-user-images/image-20260124005830555.png)





和拓扑的关系( 基本群表示出的空间性质):

在任何二维以上的空间里，自旋统计定理规定任何多粒子状态都必须要么遵循费米-狄拉克统计，要么遵循玻色-爱因斯坦统计。这与$n>2$的$SO(n,1)$基本群有关，其值为$Z_2$。因此这里只有两个可能性

在二维空间里情况发生了变化，这里SO(2,1)的基本群是Z。详细地说特殊正交群SO(2,1)的射影表示不仅仅有SO(2,1)或者其二重复盖群旋量群Spin(2,1)的线性表示。而这些额外的表示被称为任意子。

这个概念对非相对论系统也有效。关键是空间旋量群是有无限基本群的SO(2)。

在二维中两个粒子的排列群不再是对称群$S_2$，而是辫子群$B_2$了。
