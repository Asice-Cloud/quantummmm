#!/usr/bin/env python3
"""
MZM 极限验证：E₁=0, t₁=0
Fig 1(b) 基准测试
"""
import numpy as np
pi = np.pi; tc = 0.3; E0 = 0.3

def fp(t, tau): return 0.5*(1.0 + np.cos(pi*t/tau))
def fm(t, tau): return 0.5*(1.0 - np.cos(pi*t/tau))

def b1(t, tau):  # Step 1, no E₁, no t₁
    A = np.zeros((5,5))
    t2v = tc*fm(t,tau); A[1,3] = -2*t2v; A[3,1] = 2*t2v
    Edv = E0*fp(t,tau); A[3,4] = 2*Edv; A[4,3] = -2*Edv
    return A

def b2(t, tau):  # Step 2
    A = np.zeros((5,5))
    t2v = tc*fp(t,tau); A[1,3] = -2*t2v; A[3,1] = 2*t2v
    t3v = tc*fm(t,tau); A[2,3] = 2*t3v; A[3,2] = -2*t3v
    return A

def b3(t, tau):  # Step 3
    A = np.zeros((5,5))
    t3v = tc*fp(t,tau); A[2,3] = 2*t3v; A[3,2] = -2*t3v
    Edv = E0*fm(t,tau); A[3,4] = 2*Edv; A[4,3] = -2*Edv
    return A

def prop(bld, tau):
    n = max(500, int(2*tau))
    dt = tau/n; R = np.eye(5)
    for s in range(n):
        t = s*dt
        k1 = bld(t,tau) @ R
        k2 = bld(t+0.5*dt,tau) @ (R+0.5*dt*k1)
        k3 = bld(t+0.5*dt,tau) @ (R+0.5*dt*k2)
        k4 = bld(t+dt,tau) @ (R+dt*k3)
        R += dt/6.0*(k1+2*k2+2*k3+k4)
    return R

def so5_to_spinor_4d(R):
    """
    从 SO(5) 矩阵提取 4×4 复 spinor 表示中的 ancilla qubit 作用
    计算 U(6τ) 在 |ψ₁⁺⟩ = (|0₁₂⟩ + i|1₁₂⟩)/√2 上的投影
    返回 |⟨ψ₁⁻|U|ψ₁⁺⟩|
    """
    ov = 0.5*(R[0,0] + 1j*R[1,0] + 1j*R[0,1] - R[1,1])
    return np.abs(ov)

# ═══════════════════════════════════════════════════════════
print("="*60)
print("MZM 极限验证: E₁=0, t₁=0")
print("预期: 双编织后 |⟨ψ₁⁻|U(6τ)|ψ₁⁺⟩|² = 1")
print("      交换规则: γ₂→γ₃, γ₃→-γ₂")
print("="*60)

tau_list = [0.5, 1.0, 2.0, 5.0, 10.0]
for tp in tau_list:
    tc_phys = tp * 100.0
    R1 = prop(b1, tc_phys)
    R2 = prop(b2, tc_phys)
    R3 = prop(b3, tc_phys)
    U3 = R3 @ R2 @ R1
    U6 = U3 @ U3  # double swap
    
    fid_single = so5_to_spinor_4d(U3)
    fid_double = so5_to_spinor_4d(U6)
    
    # Check exchange rule on γ₂, γ₃
    # In SO(5), γ₂=index 1, γ₃=index 2
    # After single braid: γ₂→γ₃, γ₃→-γ₂
    # This means U3 maps basis vector e₁→e₂, e₂→-e₁
    e1 = np.array([0,1,0,0,0])  # γ₂
    e2 = np.array([0,0,1,0,0])  # γ₃
    e1_out = U3 @ e1
    e2_out = U3 @ e2
    
    dot_1 = np.dot(e2, e1_out)  # should be 1 (γ₂→γ₃)
    dot_2 = np.dot(e1, e2_out)  # should be -1 (γ₃→-γ₂)
    
    print(f"  τ={tp:.1f}: single fid={fid_single:.4f}, double fid={fid_double:.4f}")
    print(f"         γ₂→γ₃: {dot_1:+.4f} (expect +1), γ₃→γ₂: {dot_2:+.4f} (expect -1)")

# Also check the ancilla return
print("\n--- Ancilla qubit return check ---")
for tp in tau_list:
    tc_phys = tp * 100.0
    R1 = prop(b1, tc_phys)
    R2 = prop(b2, tc_phys)
    R3 = prop(b3, tc_phys)
    U6 = R3 @ R2 @ R1 @ R3 @ R2 @ R1
    
    # The ancilla is encoded in γ_a,γ_b (indices 3,4)
    # Check if ancilla subspace is unchanged
    ancilla_block = U6[3:5, 3:5]
    print(f"  τ={tp:.1f}: ancilla block = [[{ancilla_block[0,0]:+.4f}, {ancilla_block[0,1]:+.4f}], [{ancilla_block[1,0]:+.4f}, {ancilla_block[1,1]:+.4f}]]")

print("\n✓ Done")
