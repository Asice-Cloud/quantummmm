#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
构造并验证原方程组的通解族：
(f_x,f_y,f_z) = s(u)*v, 并把 h 分解为该方向上的分量加上任意 gauge 函数。

运行该脚本会做符号验证并给出一个具体替换示例。
"""
from sympy import symbols, Function, diff, simplify, Matrix, exp, sin, cos

u = symbols('u')
# s(u) 为任意标量函数，v_x,v_y,v_z 为常向量分量
s = Function('s')(u)
v_x, v_y, v_z = symbols('v_x v_y v_z')
# gauge 任意函数
g1 = Function('g1')(u)
g2 = Function('g2')(u)
g3 = Function('g3')(u)

# 定义 f 分量
f_x = s * v_x
f_y = s * v_y
f_z = s * v_z

# h 的分解（任取 g1,g2,g3）
h_xx = f_x/2 + g1
h_yy = f_x/2 - g1
h_xy = -f_y/2 + g2
h_yx = f_y/2 + g2
h_zi = f_z/2 + g3
h_iz = -f_z/2 + g3

# 原方程的三个标量表达式
eq1 = (diff(h_xx, u) + diff(h_yy, u)) * (-h_xy + h_yx) - (-diff(h_xy, u) + diff(h_yx, u)) * (h_xx + h_yy)

eq2 = (-diff(h_xy, u) + diff(h_yx, u)) * (h_zi - h_iz) - (diff(h_zi, u) - diff(h_iz, u)) * (-h_xy + h_yx)

eq3 = (diff(h_zi, u) - diff(h_iz, u)) * (h_xx + h_yy) - (diff(h_xx, u) + diff(h_yy, u)) * (h_zi - h_iz)

# 化简（理论上应为 0）
eq1s = simplify(eq1)
eq2s = simplify(eq2)
eq3s = simplify(eq3)


def main():
    print("符号化简结果（应为 0）：")
    print("eq1:", eq1s)
    print("eq2:", eq2s)
    print("eq3:", eq3s)

    # 交叉积检验 f × f' == 0
    f = Matrix([f_x, f_y, f_z])
    fp = Matrix([diff(f_x, u), diff(f_y, u), diff(f_z, u)])
    print("f × f' (符号表达式):", simplify(f.cross(fp)))

    # 具体示例替换：s=exp(2u), v=(1,2,3), g1=sin(u), g2=0, g3=cos(u)
    s_ex = exp(2 * u)
    subs = {s: s_ex, v_x: 1, v_y: 2, v_z: 3, g1: sin(u), g2: 0, g3: cos(u)}
    eq1_sub = simplify(eq1.subs(subs))
    eq2_sub = simplify(eq2.subs(subs))
    eq3_sub = simplify(eq3.subs(subs))

    print("\n示例替换（s=exp(2u), v=(1,2,3), g1=sin(u), g2=0, g3=cos(u)）:")
    print("eq1_sub:", eq1_sub)
    print("eq2_sub:", eq2_sub)
    print("eq3_sub:", eq3_sub)

    # 数值检查（取 u=0）
    numeric_vals = {u: 0}
    print("\n数值检查（u=0）:")
    print("eq1_sub(u=0)=", float(eq1_sub.subs(numeric_vals)))
    print("eq2_sub(u=0)=", float(eq2_sub.subs(numeric_vals)))
    print("eq3_sub(u=0)=", float(eq3_sub.subs(numeric_vals)))


if __name__ == '__main__':
    main()
