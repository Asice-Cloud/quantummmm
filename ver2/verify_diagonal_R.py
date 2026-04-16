#!/usr/bin/env python3
"""
verify_diagonal_R.py

验证对角两体生成元 H_P = Jx σx⊗σx + Jy σy⊗σy + Jz σz⊗σz
构造 R = exp(i c00) exp(i Jx σxσx) exp(i Jy σyσy) exp(i Jz σzσz)
并符号化/数值化检验 Yang–Baxter 方程 R12 R13 R23 = R23 R13 R12。

输出：对每类非零矩阵元给出因式分解并列出可行的约束族。

"""
from sympy import (
    symbols,
    Matrix,
    eye,
    cos,
    sin,
    I as symI,
    kronecker_product,
    simplify,
    factor,
)
import numpy as np

# 符号变量
c00, Jx, Jy, Jz = symbols('c00 Jx Jy Jz')

# Pauli 矩阵（SymPy）
sx = Matrix([[0, 1], [1, 0]])
sy = Matrix([[0, -symI], [symI, 0]])
sz = Matrix([[1, 0], [0, -1]])
I2 = eye(2)

# 两体算符
A = kronecker_product(sx, sx)
B = kronecker_product(sy, sy)
C = kronecker_product(sz, sz)
I4 = eye(4)

# 由于 A,B,C 互相对易，使用乘积展开闭式：exp(iJ A) = cos(J) I + i sin(J) A
Rx = cos(Jx) * I4 + symI * sin(Jx) * A
Ry = cos(Jy) * I4 + symI * sin(Jy) * B
Rz = cos(Jz) * I4 + symI * sin(Jz) * C
from sympy import exp as sympy_exp
phase = sympy_exp(symI * c00)  # exp(i c00)

R = simplify(phase * (Rx * Ry * Rz))

# 嵌入三格空间
R12 = kronecker_product(R, I2)
R23 = kronecker_product(I2, R)
# 交换算符 P（两比特）
P = Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
P23 = kronecker_product(I2, P)
R13 = simplify(P23 * R12 * P23)

# YBE 差
expr = simplify(R12 * R13 * R23 - R23 * R13 * R12)

# 收集非零条目并因式分解
nonzero = []
for i in range(expr.rows):
    for j in range(expr.cols):
        val = simplify(expr[i, j])
        if val != 0:
            nonzero.append(simplify(val))

uniq = []
seen = set()
for v in nonzero:
    s = str(v)
    if s not in seen:
        seen.add(s)
        uniq.append(v)

print('=== Yang–Baxter 差矩阵：非零条目数 =', len(nonzero), '，去重后 =', len(uniq))
print('\n代表性因子化表达（取若干唯一条目）：')
for k, v in enumerate(uniq[:8]):
    print('\n--- 条目 #', k + 1)
    fv = factor(v)
    print(fv)

# 从典型的因式中提取常见因子形式（人工解析）
# 为便于阅读，下面寻找几个常见因子模式：sin(Jy-Jx)、复数因子 (cos*... + i sin*...)
# 取第一个非零条目来显示其因子分解（若存在）
if uniq:
    first = factor(uniq[0])
    print('\n第一个非零条目的因式分解（用于解析约束的参考）:\n', first)

# 给出实数解的常见族（数值验证）
print('\n=== 数值验证（快速扫描）===')

# 数值构造 R（NumPy）
def R_num(Jxv, Jyv, Jzv, c00v=0.0):
    sx_n = np.array([[0, 1], [1, 0]], dtype=complex)
    sy_n = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz_n = np.array([[1, 0], [0, -1]], dtype=complex)
    I2_n = np.eye(2, dtype=complex)
    A_n = np.kron(sx_n, sx_n)
    B_n = np.kron(sy_n, sy_n)
    C_n = np.kron(sz_n, sz_n)
    I4_n = np.eye(4, dtype=complex)
    Rx_n = np.cos(Jxv) * I4_n + 1j * np.sin(Jxv) * A_n
    Ry_n = np.cos(Jyv) * I4_n + 1j * np.sin(Jyv) * B_n
    Rz_n = np.cos(Jzv) * I4_n + 1j * np.sin(Jzv) * C_n
    R_n = np.exp(1j * c00v) * (Rx_n @ Ry_n @ Rz_n)
    return R_n


def ybe_holds_num(Rn, tol=1e-9):
    I2_n = np.eye(2, dtype=complex)
    S = np.zeros((4, 4), dtype=complex)
    for i in range(2):
        for j in range(2):
            S[i * 2 + j, j * 2 + i] = 1
    R12 = np.kron(Rn, I2_n)
    R23 = np.kron(I2_n, Rn)
    R13 = np.kron(S, I2_n) @ R23 @ np.kron(S, I2_n)
    lhs = R12 @ R13 @ R23
    rhs = R23 @ R13 @ R12
    return np.allclose(lhs, rhs, atol=tol)

# 验证若干典型族：
# 1) 各向同性 Jx=Jy=Jz
print('各向同性 Jx=Jy=Jz 的数值检查:')
for Jv in [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]:
    Rn = R_num(Jv, Jv, Jv)
    print(' J =', float(Jv), ' -> YBE:', ybe_holds_num(Rn))

# 2) 粗网格扫描（0..pi/2）用于找出明显格点解
print('\n粗网格扫描（0..pi/2）找出满足 YBE 的格点（坐标显示）...')
sols = []
grid = np.linspace(0, np.pi / 2, 9)
for Jxv in grid:
    for Jyv in grid:
        for Jzv in grid:
            Rn = R_num(Jxv, Jyv, Jzv)
            if ybe_holds_num(Rn, tol=1e-8):
                sols.append((float(Jxv), float(Jyv), float(Jzv)))
print('找到解数量:', len(sols))
# 只打印部分例子
for s in sols[:30]:
    print(' ', s)

print('\n脚本完成。解析输出请参看以上因式分解，数值输出给出典型解样本。')
