#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DGSS_Model v5.2 数学层 SymPy 符号验证
=====================================

对应模型 : DGSS_Model.sysml  (package DGSS_Model, Version 5.2, 2026-06-10)
依赖     : Python >= 3.8, SymPy >= 1.9        安装: pip install sympy
用法     : python verify_dgss_v52_math.py
退出码   : 0 = 全部通过; 1 = 存在 FAIL

验证项与模型元素的对应关系
--------------------------------------------------------------------
  [V1]  buildTransform_DGSLCS_to_DGSCS : R == Rz(α)·Ry(β)·Rx(γ)  (ZYX, 右手定则)
  [V2]  全部旋转矩阵正交归一性          : RᵀR = I 且 det R = +1
  [V3]  composeTransforms              : 分量公式 == 分块乘法 [A|t1]·[B|t2]
  [V4]  invertTransform                : 分量公式 == (Rᵀ, −Rᵀt);  且 T⁻¹∘T = I₄
  [V5]  applyTransform                 : 分量公式 == R·p + t
  [V6]  buildTransform_MCS_to_PFCS     == invert( buildTransform_PFCS_to_MCS )
  [V7]  buildTransform_PFCS_to_DGSCS   == invert( buildTransform_DGSCS_to_PFCS )
  [V8]  buildTransform_DGSCS_to_PFCS   == compose( MCS_to_PFCS, DGSCS_to_MCS )
  [V9]  往返复合 == 恒等变换 (PFCS↔MCS, DGSCS↔PFCS)
  [V10] 标定数据数值一致性: displacementDefinitionConsistency 与
        translationConsistency (FIX⑩ 容差 1.0e-4 mm), 高精度残差报告
  [V11] DGSLCS→PFCS→MCS 全链数值抽样: 复合矩阵 vs 逐级应用 vs 逆变换还原
  [V12] TC 复合开口几何: getIdealHitLocalZ 与 pointInsideTCAperture 判定逻辑,
        理想落点与 idealCollectionPosition 自洽
  [V13] 弹道穿越点 crossingPointAtX 闭式解符号推导 + DropletBallisticBudget
        预算带与开口几何的包含关系
--------------------------------------------------------------------
约定: 所有模型字面量以 sympy.Rational 精确十进制载入, 符号比较一律
      simplify/trigsimp 至零矩阵; 数值比较使用 50 位有效数字 evalf.
"""

import sys
import sympy as sp
from sympy import Rational as Q   # 精确十进制字面量

# ====================================================================
# 0. 测试harness
# ====================================================================
_results = []

def check(name, ok, detail=""):
    _results.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    line = f"  [{tag}] {name}"
    if detail:
        line += f"   {detail}"
    print(line)

def section(title):
    print("\n" + "=" * 68)
    print(title)
    print("=" * 68)

def is_zero_matrix(M):
    """符号矩阵是否恒为零 (经 expand + trigsimp)."""
    Ms = sp.trigsimp(sp.expand_trig(sp.expand(M)))
    return all(sp.simplify(e) == 0 for e in Ms)

def mats_equal(A, B):
    return is_zero_matrix(sp.Matrix(A) - sp.Matrix(B))

# ====================================================================
# 1. 基础旋转/仿射工具 (右手定则, 与模型 RotationSignConvention 一致)
# ====================================================================
def Rx(t):
    return sp.Matrix([[1, 0, 0],
                      [0, sp.cos(t), -sp.sin(t)],
                      [0, sp.sin(t),  sp.cos(t)]])

def Ry(t):
    return sp.Matrix([[ sp.cos(t), 0, sp.sin(t)],
                      [ 0,         1, 0        ],
                      [-sp.sin(t), 0, sp.cos(t)]])

def Rz(t):
    return sp.Matrix([[sp.cos(t), -sp.sin(t), 0],
                      [sp.sin(t),  sp.cos(t), 0],
                      [0,          0,         1]])

def affine(R, t):
    """3×3 R + 3×1 t → 4×4 齐次矩阵."""
    M = sp.eye(4)
    M[:3, :3] = R
    M[:3, 3] = sp.Matrix(t)
    return M

# ====================================================================
# 2. 模型 calc def 的逐字转写 (来源: DGSS_Model.sysml v5.2)
#    每个函数返回 (R, t), 元素表达式与 .sysml 文本一一对应。
# ====================================================================

def model_DGSLCS_to_DGSCS(a, b, g, tx, ty, tz):
    """calc def buildTransform_DGSLCS_to_DGSCS
       in: alpha_Rz(yaw)=a, beta_Ry(pitch)=b, gamma_Rx(roll)=g"""
    ca, sa = sp.cos(a), sp.sin(a)
    cb, sb = sp.cos(b), sp.sin(b)
    cg, sg = sp.cos(g), sp.sin(g)
    R = sp.Matrix([
        [ca*cb,  ca*sb*sg - sa*cg,  ca*sb*cg + sa*sg],
        [sa*cb,  sa*sb*sg + ca*cg,  sa*sb*cg - ca*sg],
        [ -sb,   cb*sg,             cb*cg           ]])
    return R, sp.Matrix([tx, ty, tz])

def model_PFCS_to_MCS(a, Lopx, Lopy, Lopz):
    """calc def buildTransform_PFCS_to_MCS  —  R = Rx(α), T = L_op"""
    ca, sa = sp.cos(a), sp.sin(a)
    R = sp.Matrix([[1,  0,   0],
                   [0, ca, -sa],
                   [0, sa,  ca]])
    return R, sp.Matrix([Lopx, Lopy, Lopz])

def model_MCS_to_PFCS(a, Lopx, Lopy, Lopz):
    """calc def buildTransform_MCS_to_PFCS  —  声称 = invert(PFCS_to_MCS)"""
    ca, sa = sp.cos(a), sp.sin(a)
    R = sp.Matrix([[1,   0,  0],
                   [0,  ca, sa],
                   [0, -sa, ca]])
    t = sp.Matrix([0 - Lopx,
                   0 - Lopy*ca - Lopz*sa,
                   Lopy*sa - Lopz*ca])
    return R, t

def model_DGSCS_to_MCS(b, Losx, Losy, Losz):
    """calc def buildTransform_DGSCS_to_MCS  —  β=78° 单段, T = L_os"""
    cb, sb = sp.cos(b), sp.sin(b)
    R = sp.Matrix([[ 0,  cb,  sb],
                   [-1,   0,   0],
                   [ 0, -sb,  cb]])
    return R, sp.Matrix([Losx, Losy, Losz])

def model_DGSCS_to_PFCS(a, b, Lopx, Lopy, Lopz, Losx, Losy, Losz):
    """calc def buildTransform_DGSCS_to_PFCS  —  复合树边, 平移列应 = L_ps"""
    ca, sa = sp.cos(a), sp.sin(a)
    cb, sb = sp.cos(b), sp.sin(b)
    R = sp.Matrix([[  0,    cb,      sb   ],
                   [-ca, -sa*sb,  sa*cb   ],
                   [ sa, -ca*sb,  ca*cb   ]])
    t = sp.Matrix([ (Losx - Lopx),
                    (Losy - Lopy)*ca + (Losz - Lopz)*sa,
                    (Lopy - Losy)*sa + (Losz - Lopz)*ca ])
    return R, t

def model_PFCS_to_DGSCS(a, b, Lopx, Lopy, Lopz, Losx, Losy, Losz):
    """calc def buildTransform_PFCS_to_DGSCS  —  声称 = invert(DGSCS_to_PFCS)"""
    ca, sa = sp.cos(a), sp.sin(a)
    cb, sb = sp.cos(b), sp.sin(b)
    R = sp.Matrix([[ 0,   -ca,     sa   ],
                   [cb, -sa*sb, -ca*sb  ],
                   [sb,  sa*cb,  ca*cb  ]])
    t = sp.Matrix([ (Losy - Lopy),
                    (Lopx - Losx)*cb - (Lopz - Losz)*sb,
                    (Lopx - Losx)*sb + (Lopz - Losz)*cb ])
    return R, t

def model_composeTransforms(A, t1, B, t2):
    """calc def composeTransforms — 分量公式逐字转写: result = T1·T2"""
    a11, a12, a13 = A[0, 0], A[0, 1], A[0, 2]
    a21, a22, a23 = A[1, 0], A[1, 1], A[1, 2]
    a31, a32, a33 = A[2, 0], A[2, 1], A[2, 2]
    t1x, t1y, t1z = t1[0], t1[1], t1[2]
    b11, b12, b13 = B[0, 0], B[0, 1], B[0, 2]
    b21, b22, b23 = B[1, 0], B[1, 1], B[1, 2]
    b31, b32, b33 = B[2, 0], B[2, 1], B[2, 2]
    t2x, t2y, t2z = t2[0], t2[1], t2[2]
    R = sp.Matrix([
        [a11*b11+a12*b21+a13*b31, a11*b12+a12*b22+a13*b32, a11*b13+a12*b23+a13*b33],
        [a21*b11+a22*b21+a23*b31, a21*b12+a22*b22+a23*b32, a21*b13+a22*b23+a23*b33],
        [a31*b11+a32*b21+a33*b31, a31*b12+a32*b22+a33*b32, a31*b13+a32*b23+a33*b33]])
    t = sp.Matrix([a11*t2x + a12*t2y + a13*t2z + t1x,
                   a21*t2x + a22*t2y + a23*t2z + t1y,
                   a31*t2x + a32*t2y + a33*t2z + t1z])
    return R, t

def model_invertTransform(Rm, tm):
    """calc def invertTransform — 分量公式逐字转写: R' = Rᵀ, t' = −Rᵀt"""
    r11, r12, r13 = Rm[0, 0], Rm[0, 1], Rm[0, 2]
    r21, r22, r23 = Rm[1, 0], Rm[1, 1], Rm[1, 2]
    r31, r32, r33 = Rm[2, 0], Rm[2, 1], Rm[2, 2]
    tx, ty, tz = tm[0], tm[1], tm[2]
    R = sp.Matrix([[r11, r21, r31],
                   [r12, r22, r32],
                   [r13, r23, r33]])
    t = sp.Matrix([-(r11*tx + r21*ty + r31*tz),
                   -(r12*tx + r22*ty + r32*tz),
                   -(r13*tx + r23*ty + r33*tz)])
    return R, t

def model_applyTransform(Rm, tm, p):
    """calc def applyTransform — 分量公式逐字转写: p' = R·p + t"""
    r11, r12, r13 = Rm[0, 0], Rm[0, 1], Rm[0, 2]
    r21, r22, r23 = Rm[1, 0], Rm[1, 1], Rm[1, 2]
    r31, r32, r33 = Rm[2, 0], Rm[2, 1], Rm[2, 2]
    tx, ty, tz = tm[0], tm[1], tm[2]
    px, py, pz = p[0], p[1], p[2]
    return sp.Matrix([r11*px + r12*py + r13*pz + tx,
                      r21*px + r22*py + r23*pz + ty,
                      r31*px + r32*py + r33*pz + tz])

# ====================================================================
# 3. 通用符号
# ====================================================================
al, be, ga = sp.symbols('alpha beta gamma', real=True)        # 命令角 / 安装角
Lopx, Lopy, Lopz = sp.symbols('L_opx L_opy L_opz', real=True)
Losx, Losy, Losz = sp.symbols('L_osx L_osy L_osz', real=True)
txs, tys, tzs = sp.symbols('t_x t_y t_z', real=True)

# ====================================================================
# [V1] DGSLCS→DGSCS 旋转 == Rz(α)·Ry(β)·Rx(γ)
# ====================================================================
section("[V1] buildTransform_DGSLCS_to_DGSCS : R == Rz(α)·Ry(β)·Rx(γ)")
R_model, t_model = model_DGSLCS_to_DGSCS(al, be, ga, txs, tys, tzs)
R_ref = Rz(al) * Ry(be) * Rx(ga)
check("旋转 3×3 与 ZYX 欧拉参考矩阵符号相等", mats_equal(R_model, R_ref))
check("平移列 = (tx, ty, tz) 直通",
      mats_equal(t_model, sp.Matrix([txs, tys, tzs])))

# ====================================================================
# [V2] 全部旋转矩阵正交归一性 + det = +1
# ====================================================================
section("[V2] 正交归一性 RᵀR = I 与 det R = +1 (全部构造函数)")
builders = {
    "DGSLCS_to_DGSCS (ZYX 命令位姿)":
        model_DGSLCS_to_DGSCS(al, be, ga, txs, tys, tzs)[0],
    "PFCS_to_MCS  (Rx(α))":
        model_PFCS_to_MCS(al, Lopx, Lopy, Lopz)[0],
    "MCS_to_PFCS":
        model_MCS_to_PFCS(al, Lopx, Lopy, Lopz)[0],
    "DGSCS_to_MCS (β 单段)":
        model_DGSCS_to_MCS(be, Losx, Losy, Losz)[0],
    "DGSCS_to_PFCS (复合树边)":
        model_DGSCS_to_PFCS(al, be, Lopx, Lopy, Lopz, Losx, Losy, Losz)[0],
    "PFCS_to_DGSCS":
        model_PFCS_to_DGSCS(al, be, Lopx, Lopy, Lopz, Losx, Losy, Losz)[0],
}
for name, R in builders.items():
    ortho = mats_equal(R.T * R, sp.eye(3))
    det1 = sp.simplify(sp.trigsimp(R.det())) == 1
    check(f"{name}: RᵀR = I", ortho)
    check(f"{name}: det R = +1 (右手系, 无反射)", det1)

# ====================================================================
# [V3] composeTransforms 分量公式 == 分块矩阵乘法
# ====================================================================
section("[V3] composeTransforms : 分量公式 == [A|t1]·[B|t2] 分块乘法")
Asym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'a{i+1}{j+1}', real=True))
Bsym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'b{i+1}{j+1}', real=True))
t1s = sp.Matrix(sp.symbols('p1 p2 p3', real=True))
t2s = sp.Matrix(sp.symbols('q1 q2 q3', real=True))
Rc, tc = model_composeTransforms(Asym, t1s, Bsym, t2s)
check("旋转块: R = A·B  (24 个独立符号, 完全一般性)",
      mats_equal(Rc, Asym * Bsym))
check("平移块: t = A·t2 + t1", mats_equal(tc, Asym * t2s + t1s))
M4 = affine(Asym, t1s) * affine(Bsym, t2s)
check("与 4×4 齐次矩阵乘法逐元素一致",
      mats_equal(affine(Rc, tc), M4))

# ====================================================================
# [V4] invertTransform : 公式结构 + 刚体逆的真逆性
# ====================================================================
section("[V4] invertTransform : R' = Rᵀ, t' = −Rᵀt ;  且 T⁻¹∘T = I₄")
Rsym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'r{i+1}{j+1}', real=True))
ts = sp.Matrix(sp.symbols('s1 s2 s3', real=True))
Ri, ti = model_invertTransform(Rsym, ts)
check("R' = Rᵀ (一般 9 符号)", mats_equal(Ri, Rsym.T))
check("t' = −Rᵀ·t (一般 12 符号)", mats_equal(ti, -Rsym.T * ts))
# 真逆性需要 R 正交: 用 ZYX 欧拉参数化的一般旋转验证 T⁻¹∘T = I
R_e = Rz(al) * Ry(be) * Rx(ga)
Rinv_e, tinv_e = model_invertTransform(R_e, ts)
Rrt, trt = model_composeTransforms(Rinv_e, tinv_e, R_e, ts)
check("T⁻¹∘T = I₄ (R 取一般 ZYX 旋转, t 一般符号)",
      mats_equal(Rrt, sp.eye(3)) and mats_equal(trt, sp.zeros(3, 1)))
Rrt2, trt2 = model_composeTransforms(R_e, ts, Rinv_e, tinv_e)
check("T∘T⁻¹ = I₄ (左右互逆)",
      mats_equal(Rrt2, sp.eye(3)) and mats_equal(trt2, sp.zeros(3, 1)))

# ====================================================================
# [V5] applyTransform 分量公式 == R·p + t
# ====================================================================
section("[V5] applyTransform : 分量公式 == R·p + t")
psym = sp.Matrix(sp.symbols('px py pz', real=True))
check("p' = R·p + t (一般 15 符号)",
      mats_equal(model_applyTransform(Rsym, ts, psym), Rsym * psym + ts))

# ====================================================================
# [V6] MCS_to_PFCS == invert(PFCS_to_MCS)   (符号, α 与 L_op 全为符号)
# ====================================================================
section("[V6] buildTransform_MCS_to_PFCS == invert( buildTransform_PFCS_to_MCS )")
R_pm, t_pm = model_PFCS_to_MCS(al, Lopx, Lopy, Lopz)
R_mp, t_mp = model_MCS_to_PFCS(al, Lopx, Lopy, Lopz)
R_inv, t_inv = model_invertTransform(R_pm, t_pm)
check("旋转块相等  R_MCS→PFCS == Rx(α)ᵀ", mats_equal(R_mp, R_inv))
check("平移块相等  t == −Rx(α)ᵀ·L_op", mats_equal(t_mp, t_inv))

# ====================================================================
# [V7] PFCS_to_DGSCS == invert(DGSCS_to_PFCS)   (符号)
# ====================================================================
section("[V7] buildTransform_PFCS_to_DGSCS == invert( buildTransform_DGSCS_to_PFCS )")
R_sp_, t_sp_ = model_DGSCS_to_PFCS(al, be, Lopx, Lopy, Lopz, Losx, Losy, Losz)
R_ps_, t_ps_ = model_PFCS_to_DGSCS(al, be, Lopx, Lopy, Lopz, Losx, Losy, Losz)
R_inv2, t_inv2 = model_invertTransform(R_sp_, t_sp_)
check("旋转块相等", mats_equal(R_ps_, R_inv2))
check("平移块相等", mats_equal(t_ps_, t_inv2))

# ====================================================================
# [V8] DGSCS_to_PFCS == compose(MCS_to_PFCS, DGSCS_to_MCS)   (符号)
#      —— 验证 "经 MCS 并入复合" 的注释声明
# ====================================================================
section("[V8] buildTransform_DGSCS_to_PFCS == compose( MCS_to_PFCS, DGSCS_to_MCS )")
R_sm, t_sm = model_DGSCS_to_MCS(be, Losx, Losy, Losz)
R_cmp, t_cmp = model_composeTransforms(R_mp, t_mp, R_sm, t_sm)
check("旋转块相等  (α, β 全符号)", mats_equal(R_sp_, R_cmp))
check("平移块相等  (= L_ps 解析式)", mats_equal(t_sp_, t_cmp))

# ====================================================================
# [V9] 往返复合 == 恒等
# ====================================================================
section("[V9] 往返复合 == 恒等变换")
Ra, ta = model_composeTransforms(R_mp, t_mp, R_pm, t_pm)
check("MCS_to_PFCS ∘ PFCS_to_MCS = I₄",
      mats_equal(Ra, sp.eye(3)) and mats_equal(ta, sp.zeros(3, 1)))
Rb, tb = model_composeTransforms(R_ps_, t_ps_, R_sp_, t_sp_)
check("PFCS_to_DGSCS ∘ DGSCS_to_PFCS = I₄",
      mats_equal(Rb, sp.eye(3)) and mats_equal(tb, sp.zeros(3, 1)))

# ====================================================================
# [V10] 标定数据数值一致性 (DGSInstallationCalibration)
# ====================================================================
section("[V10] 标定数值: displacementDefinitionConsistency / translationConsistency")
DEG = sp.pi / 180
alpha_v = 28 * DEG                  # alpha_install = 28°
beta_v  = 78 * DEG                  # beta_install  = 78°
P_MCS = sp.Matrix([Q('0.000000'),   Q('-41.524524'), Q('1248.672688')])
S_MCS = sp.Matrix([Q('916.328834'), Q('-41.524524'), Q('1341.435967')])
Lop = sp.Matrix([Q('0.000000'),   Q('-41.524524'), Q('1248.672688')])
Los = sp.Matrix([Q('916.328834'), Q('-41.524524'), Q('1341.435967')])
Lps = sp.Matrix([Q('916.328834'), Q('43.549721'),  Q('81.905114')])
TOL_MM = Q('1.0e-4')                # FIX⑩ 容差 0.1 µm

check("displacementDefinitionConsistency: L_op == P_in_MCS (精确)",
      Lop == P_MCS)
check("displacementDefinitionConsistency: L_os == S_in_MCS (精确)",
      Los == S_MCS)

ca28, sa28 = sp.cos(alpha_v), sp.sin(alpha_v)
expect = sp.Matrix([
    Los[0] - Lop[0],
    (Los[1] - Lop[1]) * ca28 + (Los[2] - Lop[2]) * sa28,
    (Lop[1] - Los[1]) * sa28 + (Los[2] - Lop[2]) * ca28])
res = (Lps - expect).applyfunc(lambda e: sp.Abs(e).evalf(50))
labels = ("Lpsx", "Lpsy", "Lpsz")
for i in range(3):
    nm = (res[i] * Q('1e6')).evalf(6)     # mm → nm
    check(f"translationConsistency: |Δ{labels[i]}| <= 1.0e-4 mm",
          res[i] < sp.Float(TOL_MM, 50),
          f"残差 = {nm} nm")
print("  注: FIX⑩ 注释称 L_ps.y 精确值 ≈ 43.549722, 与 6 位手填字面量差约 0.6 nm,")
print(f"      本脚本高精度复算 L_ps.y(精确) = {expect[1].evalf(12)} mm")

# 复合树边的平移列应在容差内重现 L_ps
R_num, t_num = model_DGSCS_to_PFCS(alpha_v, beta_v, *Lop, *Los)
res2 = (sp.Matrix(Lps) - t_num).applyfunc(lambda e: sp.Abs(e).evalf(50))
check("buildTransform_DGSCS_to_PFCS 平移列 == L_ps (容差 0.1 µm)",
      all(r < sp.Float(TOL_MM, 50) for r in res2))

# ====================================================================
# [V11] DGSLCS→PFCS→MCS 全链数值抽样
# ====================================================================
section("[V11] 全链数值抽样: 复合矩阵 vs 逐级应用 vs 逆变换还原")
# 取行程极限内的命令位姿: |tx|<=9, |ty|<=10, |tz|<=4.5 mm;
#                         |roll|<=2.1°, |pitch|<=2.4°, |yaw|<=4.3°
yaw_c, pitch_c, roll_c = Q('2.0') * DEG, Q('-1.2') * DEG, Q('0.5') * DEG
tx_c, ty_c, tz_c = Q('1.5'), Q('-2.0'), Q('0.7')
P_l = sp.Matrix([Q('3'), Q('-4'), Q('5')])           # DGSLCS 下测试点 (mm)

R_sl, t_sl = model_DGSLCS_to_DGSCS(yaw_c, pitch_c, roll_c, tx_c, ty_c, tz_c)
R_ps2, t_ps2 = model_DGSCS_to_PFCS(alpha_v, beta_v, *Lop, *Los)
# buildTransform_DGSLCS_to_PFCS = composeTransforms(T_ps, T_sl)
R_pl, t_pl = model_composeTransforms(R_ps2, t_ps2, R_sl, t_sl)

p_via_composed = model_applyTransform(R_pl, t_pl, P_l)
p_stepwise = model_applyTransform(
    R_ps2, t_ps2, model_applyTransform(R_sl, t_sl, P_l))
diff1 = (p_via_composed - p_stepwise).applyfunc(
    lambda e: sp.Abs(e.evalf(50)))
check("复合矩阵一次映射 == 两级逐次映射 (DGSLCS→DGSCS→PFCS)",
      all(d < sp.Float('1e-30', 40) for d in diff1),
      f"max|Δ| = {float(max(diff1)):.1e} mm")

# 继续上行至 MCS, 再经逆变换全部还原
R_pm2, t_pm2 = model_PFCS_to_MCS(alpha_v, *Lop)
p_MCS = model_applyTransform(R_pm2, t_pm2, p_via_composed)
R_mp2, t_mp2 = model_MCS_to_PFCS(alpha_v, *Lop)
R_sp2, t_sp2 = model_PFCS_to_DGSCS(alpha_v, beta_v, *Lop, *Los)
R_ls, t_ls = model_invertTransform(R_sl, t_sl)
p_back = model_applyTransform(
    R_ls, t_ls, model_applyTransform(
        R_sp2, t_sp2, model_applyTransform(R_mp2, t_mp2, p_MCS)))
diff2 = (p_back - P_l).applyfunc(lambda e: sp.Abs(e.evalf(50)))
check("MCS → PFCS → DGSCS → DGSLCS 逆链还原原始点",
      all(d < sp.Float('1e-30', 40) for d in diff2),
      f"max|Δ| = {float(max(diff2)):.1e} mm")
print(f"  抽样点 [P]_DGSLCS = (3, -4, 5) mm  →  [P]_PFCS = "
      f"({p_via_composed[0].evalf(10)}, {p_via_composed[1].evalf(10)}, "
      f"{p_via_composed[2].evalf(10)}) mm")

# ====================================================================
# [V12] TC 复合开口几何 与 理想落点
# ====================================================================
section("[V12] TCApertureGeometry / getIdealHitLocalZ / pointInsideTCAperture")
R_semi = Q('33')          # upperSemicircleRadius  [mm]
W_rect = Q('66')          # lowerRectangleWidth    [mm]
H_rect = Q('20.8873')     # lowerRectangleHeight   [mm]
offset = Q('14.9')        # idealHitOffsetFromBottom [mm]
tc_center = sp.Matrix([Q('-310.0930'), Q('0'), Q('5.9873')])  # PFCS
ideal_collect = sp.Matrix([Q('-310.0930'), Q('0'), Q('0')])   # PFCS

ideal_local_z = offset - H_rect           # calc def getIdealHitLocalZ
check("getIdealHitLocalZ = offset − H = −5.9873 mm (精确)",
      ideal_local_z == Q('-5.9873'),
      f"值 = {sp.nsimplify(ideal_local_z)} mm")
check("理想落点 PFCS-z = center.z + localZ == idealCollectionPosition.z (= 0)",
      tc_center[2] + ideal_local_z == ideal_collect[2])
check("开口中心 x 与 idealCollectionPosition.x 一致 (−310.0930 mm)",
      tc_center[0] == ideal_collect[0])
check("矩形宽度 == 半圆直径 (W = 2R, 共享上边)", W_rect == 2 * R_semi)

def point_inside_tc(yp, zp):
    """calc def pointInsideTCAperture 逐字转写 (局部坐标, mm)."""
    upper = (zp >= 0) and (yp*yp + zp*zp <= R_semi*R_semi)
    lower = (zp < 0) and (zp >= -H_rect) and \
            (yp >= -W_rect/2) and (yp <= W_rect/2)
    return upper or lower

cases = [
    # (yp, zp, 期望, 说明)
    (Q(0),        ideal_local_z, True,  "理想落点 (0, −5.9873) ∈ 下矩形"),
    (Q(0),        Q(0),          True,  "交界中心 (0, 0) ∈ 上半圆边界"),
    (Q(0),        R_semi,        True,  "半圆顶点 (0, 33) 边界含"),
    (Q(0),        Q('33.0001'),  False, "半圆顶点外 0.1 µm → 拒"),
    (W_rect/2,    -Q('0.0001'),  True,  "矩形右边界 (33, 0⁻) 含"),
    (Q('33.0001'),-Q(1),         False, "矩形右边界外 → 拒"),
    (Q(0),        -H_rect,       True,  "矩形底边 (0, −20.8873) 含"),
    (Q(0),        Q('-20.8874'), False, "底边下方 0.1 µm → 拒"),
    (Q('32.9'),   Q(1),          True,  "上区域圆内点 (32.9, 1)"),
    (Q('32.9'),   Q(3),          False, "上区域圆外角隙 (32.9, 3) → 拒"),
    (-W_rect/2,   -H_rect,       True,  "矩形左下角 (−33, −20.8873) 含"),
]
for yp, zp, exp, msg in cases:
    check(f"pointInsideTCAperture: {msg}", point_inside_tc(yp, zp) == exp)

# ====================================================================
# [V13] 弹道穿越点闭式解 + DropletBallisticBudget 几何自洽
# ====================================================================
section("[V13] crossingPointAtX 闭式解推导 与 弹道预算带几何包含性")
# 模型: P(t) = P0 + V0·t + ½·g·t²,  g = gravity·(0,0,−1)  (EnvironmentCondition)
t_ = sp.symbols('t', positive=True)
x0, y0, z0, vx, vy, vz, g_, xp_ = sp.symbols(
    'x0 y0 z0 v_x v_y v_z g x_p', real=True)
x_t = x0 + vx * t_                       # 重力无 x 分量 → x 线性
y_t = y0 + vy * t_
z_t = z0 + vz * t_ - sp.Rational(1, 2) * g_ * t_**2
t_star = sp.solve(sp.Eq(x_t, xp_), t_)[0]
y_star = sp.simplify(y_t.subs(t_, t_star))
z_star = sp.simplify(z_t.subs(t_, t_star))
check("穿越时间闭式解 t* = (x_p − x0)/v_x",
      sp.simplify(t_star - (xp_ - x0) / vx) == 0)
check("y* = y0 + v_y·(x_p−x0)/v_x  (代回轨迹方程残差为零)",
      sp.simplify(y_star - (y0 + vy * (xp_ - x0) / vx)) == 0)
check("z* = z0 + v_z·τ − ½gτ²,  τ=(x_p−x0)/v_x  (代回残差为零)",
      sp.simplify(
          z_star - (z0 + vz*(xp_-x0)/vx
                    - g_*(xp_-x0)**2 / (2*vx**2))) == 0)
print("  → 该闭式解即模型中抽象 calc def crossingPointAtX 的解析落地形式,")
print("    可直接用于 DropletPassesTarget / DropletPassesMonitor 约束求值。")

# DropletBallisticBudget:
#   nominalDropBudget : |z_TC − center.z| <= H/2
#   dispersionBudget  : 3σ_y <= W/6,  3σ_z <= H/6
band_lo = tc_center[2] - H_rect / 2
band_hi = tc_center[2] + H_rect / 2
aper_lo = tc_center[2] - H_rect          # 开口竖直下界 (y=0 截面)
aper_hi = tc_center[2] + R_semi          # 开口竖直上界 (y=0 截面)
check("名义下坠预算带 [center−H/2, center+H/2] ⊂ 开口竖直范围 (y=0 截面)",
      (band_lo >= aper_lo) and (band_hi <= aper_hi),
      f"带 = [{band_lo.evalf(8)}, {band_hi.evalf(8)}] mm ⊂ "
      f"[{aper_lo.evalf(8)}, {aper_hi.evalf(8)}] mm")
check("理想落点位于预算带内", band_lo <= ideal_collect[2] <= band_hi)
sigma_y_max = W_rect / 18                # 3σ_y <= W/6  →  σ_y <= W/18
sigma_z_max = H_rect / 18
print(f"  散布预算上限: 3σ_y ≤ W/6 = {(W_rect/6).evalf(6)} mm "
      f"(σ_y ≤ {sigma_y_max.evalf(6)} mm); "
      f"3σ_z ≤ H/6 = {(H_rect/6).evalf(6)} mm "
      f"(σ_z ≤ {sigma_z_max.evalf(6)} mm)")
check("预算带宽 + 3σ_z 余量仍不越下边界: H/2 + H/6 <= H",
      H_rect/2 + H_rect/6 <= H_rect)

# ====================================================================
# 汇总
# ====================================================================
section("汇总")
total = len(_results)
failed = [n for n, ok in _results if not ok]
print(f"  共 {total} 项检查, 通过 {total - len(failed)} 项, 失败 {len(failed)} 项.")
if failed:
    print("  失败项:")
    for n in failed:
        print(f"    - {n}")
    sys.exit(1)
print("  ✔ DGSS_Model v5.2 数学层全部符号/数值验证通过。")
sys.exit(0)
