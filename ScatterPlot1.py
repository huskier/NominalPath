import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. 系统参数定义 (单位统一使用国际标准单位 m, s, kg)
# ---------------------------------------------------------
#N_droplets = 36000          # 模拟360颗液滴，刚好对应空心圆锥周向每1度一颗
N_droplets = 1          # 模拟360颗液滴，刚好对应空心圆锥周向每1度一颗
v0 = 40.0                 # 初始速度: 40 m/s
g = 9.81                  # 重力加速度: 9.81 m/s^2

# 关键截面的 X 坐标 (米)
x_nozzle = 410.453510 / 1000.0   # 0.41045 m
x_target = 0.0                   # 0 m
x_tc = -310.0 / 1000.0           # -0.310 m

# 喷嘴初始空间位置 (严格处于标称中心轴线，无DGS补偿)
y0_nozzle = 0.0
z0_nozzle = 0.0

# 锥角新设置：所有液滴严格位于半锥角 theta 边界上
half_cone_angle_deg = 0.75
theta = np.radians(half_cone_angle_deg) # 固定的半锥角

# ---------------------------------------------------------
# 2. 空心锥样本生成 (液滴严格散布在 0.75° 边界圆周上)
# ---------------------------------------------------------
# 在 0 到 2*pi 之间均匀生成方位角，模拟完整的空心圆锥壳
phi = np.linspace(0, 2 * np.pi, N_droplets, endpoint=False)

# 喷嘴向 -X 方向喷射，液滴速度向量的三个分量
vx = -v0 * np.cos(theta) * np.ones_like(phi) # 所有液滴的 X 速度相同
vy = v0 * np.sin(theta) * np.cos(phi)
vz = v0 * np.sin(theta) * np.sin(phi)

# ---------------------------------------------------------
# 3. 轨迹交点解析计算（纯重力抛物线影响）
# ---------------------------------------------------------
# ---- A. 计算在 Target 截面 (X = 0) 处的落点 ----
t_target = (x_target - x_nozzle) / vx
y_target = y0_nozzle + vy * t_target
z_target = z0_nozzle + vz * t_target - 0.5 * g * (t_target**2)

# ---- B. 计算在 TC 截面 (X = -310mm) 处的落点 ----
t_tc = (x_tc - x_nozzle) / vx
y_tc = y0_nozzle + vy * t_tc
z_tc = z0_nozzle + vz * t_tc - 0.5 * g * (t_tc**2)

# ---------------------------------------------------------
# 4. 数据可视化 (绘图)
# ---------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# 图 1：Target 截面落点
ax1.scatter(y_target * 1000, z_target * 1000, s=8, c='darkblue', alpha=0.8, label='Droplets (Hollow Ring)')
ax1.plot(0, 0, 'k+', markersize=12, label='Nominal Target (0,0)')
mean_z_target = np.mean(z_target) * 1000
ax1.plot(0, mean_z_target, 'ro', markersize=6, label=f'Ring Center (0, {mean_z_target:.2f})')

ax1.set_title("Droplet Distribution at Target Position (X = 0 mm)\n(No Compensation - Pure Physics)")
ax1.set_xlabel("Y Position (mm)")
ax1.set_ylabel("Z Position (mm)")
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.axis('equal')
ax1.legend()

# 图 2：TC 截面落点 (散点图)
ax2.scatter(y_tc * 1000, z_tc * 1000, s=8, c='crimson', alpha=0.8, label='Droplets (Hollow Ring)')
ax2.plot(0, 0, 'k+', markersize=12, label='Nominal Axis (0,0)')
mean_z_tc = np.mean(z_tc) * 1000
ax2.plot(0, mean_z_tc, 'go', markersize=6, label=f'Ring Center (0, {mean_z_tc:.2f})')

# 计算空心圆环在 TC 处的半径
ring_radius_tc = np.max(y_tc) * 1000

ax2.set_title("Droplet Distribution at TC Location (X = -310 mm)\n(Final Scatter Plot)")
ax2.set_xlabel("Y Position (mm)")
ax2.set_ylabel("Z Position (mm)")
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.axis('equal')
ax2.legend()

plt.tight_layout()
plt.show()

# ---------------------------------------------------------
# 5. 定量输出
# ---------------------------------------------------------
print(f"--- 修正版数字化仿真定量评估简报 ---")
print(f"模拟模式: 空心圆锥壳（Hollow Cone）边界模拟")
print(f"半锥角 (Half Cone Angle): {half_cone_angle_deg}°")
print(f"液滴从喷嘴到 TC 的飞行总时间: {t_tc[0]*1000:.2f} ms")
print(f"---------------------------------------------------")
print(f"【Target 截面 (X=0) 特征】:")
print(f"  - 几何形态: 完美的圆形环")
print(f"  - 圆环半径: {np.max(y_target)*1000:.2f} mm")
print(f"  - 圆环中心下沉量 (因重力): {mean_z_target:.4f} mm")
print(f"【TC 截面 (X=-310mm) 特征】:")
print(f"  - 几何形态: 完美的圆形环")
print(f"  - 圆环半径: {ring_radius_tc:.2f} mm")
print(f"  - 圆环中心下沉量 (因重力): {mean_z_tc:.4f} mm")
print(f"  - 液滴绝对包络边界 (最低点): Z = {(mean_z_tc - ring_radius_tc):.2f} mm")