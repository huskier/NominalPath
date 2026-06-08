#!/usr/bin/env python3
"""
DGSS 坐标变换链数值验证脚本 (verification artifact for DGSS_Model v4)
====================================================================
依据附录给出的四坐标系往返变换矩阵, 验证:
  1. PFCS ↔ MCS   (α=28°)  互逆
  2. DGSCS ↔ MCS  (β=78°)  互逆
  3. DGSCS → PFCS 复合     = 链式 T_po · T_os, 且平移列 = L_ps
  4. DGSCS ↔ PFCS          互逆 (⑤·⑥ = I)
  5. 完整链往返闭合
  6. 标定一致性约束 (附录 T 分量公式) 精确再现 L_ps

坐标系角色 (已与设计确认):
  α/β 描述 DGSCS 之上的固定安装关系; DGS 命令的实时运动叠加在 DGSLCS→DGSCS 段.
  五环节: DGSLCS -(命令6DOF)-> DGSCS -(β=78°)-> MCS -(α=28°)-> PFCS
"""
import numpy as np

deg = np.pi / 180.0
alpha = 28 * deg      # PFCS ↔ MCS
beta  = 78 * deg      # MCS ↔ DGSCS
ca, sa = np.cos(alpha), np.sin(alpha)
cb, sb = np.cos(beta),  np.sin(beta)

# ---- 标定常量 (附录底部) ----
# L_op: PF 点在机器系(origin); L_os: DGS 原点在机器系; L_ps: DGS 原点相对 PF 在 PFCS
Lopx, Lopy, Lopz =    0.000000, -41.524524, 1248.672688
Losx, Losy, Losz =  916.328834, -41.524524, 1341.435967
Lpsx, Lpsy, Lpsz =  916.328834,  43.549721,   81.905114
L_ps_given = np.array([Lpsx, Lpsy, Lpsz])

ATOL = 1e-5

def H(R, t):
    M = np.eye(4); M[:3, :3] = R; M[:3, 3] = t; return M

def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    return cond

all_ok = True

print("="*64)
print("1. PFCS ↔ MCS  (α = 28°)")
print("="*64)
# ① PF系(p) -> 机器系(o):  Rx(α) + 平移 L_op
R_op = np.array([[1,0,0],[0,ca,-sa],[0,sa,ca]])
T_o_p = H(R_op, [Lopx, Lopy, Lopz])
# ② 机器系(o) -> PF系(p)
R_po = np.array([[1,0,0],[0,ca,sa],[0,-sa,ca]])
t_po = [-Lopx, -Lopy*ca - Lopz*sa, Lopy*sa - Lopz*ca]
T_p_o = H(R_po, t_po)
all_ok &= check("② = ①^{-1}", np.allclose(T_p_o, np.linalg.inv(T_o_p), atol=ATOL))

print("="*64)
print("2. DGSCS ↔ MCS  (β = 78°)")
print("="*64)
# ③ DGS系(s) -> 机器系(o)
R_os = np.array([[0,cb,sb],[-1,0,0],[0,-sb,cb]])
T_o_s = H(R_os, [Losx, Losy, Losz])
# ④ 机器系(o) -> DGS系(s)
R_so = np.array([[0,-1,0],[cb,0,-sb],[sb,0,cb]])
t_so = [Losy, -Losx*cb + Losz*sb, -Losx*sb - Losz*cb]
T_s_o = H(R_so, t_so)
all_ok &= check("④ = ③^{-1}", np.allclose(T_s_o, np.linalg.inv(T_o_s), atol=ATOL))

print("="*64)
print("3. DGSCS → PFCS 复合一致性 + 平移列 = L_ps")
print("="*64)
T_p_s_chain = T_p_o @ T_o_s
R_ps = np.array([[0,cb,sb],[-ca,-sa*sb,sa*cb],[sa,-ca*sb,ca*cb]])
Tx = Losx - Lopx
Ty = (Losy-Lopy)*ca + (Losz-Lopz)*sa
Tz = (Lopy-Losy)*sa + (Losz-Lopz)*ca
T_p_s_formula = H(R_ps, [Tx, Ty, Tz])
all_ok &= check("链式 = 附录公式", np.allclose(T_p_s_chain, T_p_s_formula, atol=ATOL))
all_ok &= check("平移列 T = L_ps (附录给定值)",
                np.allclose([Tx,Ty,Tz], L_ps_given, atol=1e-4))
print(f"      T = ({Tx:.6f}, {Ty:.6f}, {Tz:.6f})")
print(f"  L_ps = ({Lpsx:.6f}, {Lpsy:.6f}, {Lpsz:.6f})")
print(f"  最大残差: {np.max(np.abs(np.array([Tx,Ty,Tz]) - L_ps_given)):.2e} mm")

print("="*64)
print("4. DGSCS ↔ PFCS 互逆 (⑤·⑥ = I)")
print("="*64)
R_sp = np.array([[0,-ca,sa],[cb,-sa*sb,-ca*sb],[sb,sa*cb,ca*cb]])
Tpx = Losy - Lopy
Tpy = (Lopx-Losx)*cb - (Lopz-Losz)*sb
Tpz = (Lopx-Losx)*sb + (Lopz-Losz)*cb
T_s_p = H(R_sp, [Tpx, Tpy, Tpz])
all_ok &= check("⑥ = ⑤^{-1}", np.allclose(T_s_p, np.linalg.inv(T_p_s_formula), atol=ATOL))
all_ok &= check("⑤·⑥ = I", np.allclose(T_p_s_formula @ T_s_p, np.eye(4), atol=ATOL))

print("="*64)
print("5. 完整链往返闭合 (任取 DGS 点, s→o→p→o→s)")
print("="*64)
p_s = np.array([12.3, -4.5, 7.8, 1.0])
p_p = T_p_o @ T_o_s @ p_s
back = np.linalg.inv(T_p_o @ T_o_s) @ p_p
all_ok &= check("往返还原", np.allclose(back, p_s, atol=ATOL))
print(f"      DGS点 {p_s[:3]}  ->  PF系 ({p_p[0]:.4f}, {p_p[1]:.4f}, {p_p[2]:.4f})  ->  还原 {back[:3]}")

print("="*64)
print("6. 标定一致性约束 (= SysML translationConsistency)")
print("="*64)
c1 = np.isclose(Lpsx, Losx - Lopx, atol=1e-4)
c2 = np.isclose(Lpsy, (Losy-Lopy)*ca + (Losz-Lopz)*sa, atol=1e-4)
c3 = np.isclose(Lpsz, (Lopy-Losy)*sa + (Losz-Lopz)*ca, atol=1e-4)
all_ok &= check("Lpsx == Losx - Lopx", c1)
all_ok &= check("Lpsy == (Losy-Lopy)cosα + (Losz-Lopz)sinα", c2)
all_ok &= check("Lpsz == (Lopy-Losy)sinα + (Losz-Lopz)cosα", c3)

print()
print("="*64)
print(f"总体结果: {'全部通过 ✓' if all_ok else '存在失败 ✗'}")
print("="*64)
exit(0 if all_ok else 1)
