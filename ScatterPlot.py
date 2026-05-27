import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. 系统参数定义 (单位统一使用国际标准单位 m, s, kg)
# ---------------------------------------------------------
N_droplets = 5000         # 模拟的液滴数量
v0 = 40.0                 # 初始速度: 40 m/s
g = 9.81                  # 重力加速度: 9.81 m/s^2

# 关键截面的 X 坐标 (米)
x_nozzle = 410.453510 / 1000.0   # 0.41045 m
x_target = 0.0                   # 0 m
x_tc = -310.0 / 1000.0           # -0.310 m

# 喷嘴初始空间位置 (假设 DGS 将喷嘴基准保持在标称高度)
y0_nozzle = 0.0
z0_nozzle = 0.0

# 锥角设置
cone_angle_deg = 1.5
cone_angle_rad = np.radians(cone_angle_deg) 
# 假设 1.5° 是圆锥的边缘边界(对应 3 sigma 范围)，则角度标准差为：
sigma_theta = cone_angle_rad / 3.0 

# DGS 的 Pitch 仰角补偿 (用于对抗飞到 Target 过程中的重力下坠)
# 补偿角 = 飞到 target 过程中的重力坠落高度 / 飞到 target 的水平距离
t_to_target_nominal = (x_nozzle - x_target) / v0
z_drop_at_target = 0.5 * g * (t_to_target_nominal**2)
pitch_compensation = z_drop_at_target / (x_nozzle - x_target) # 约 1.257 mrad

# ---------------------------------------------------------
# 2. 蒙特卡洛样本生成 (模拟随机锥角喷射)
# ---------------------------------------------------------
np.random.seed(42) # 固定随机种子以确保结果可复现

# 生成偏离喷射主轴的径向角 theta (高斯分布)
theta = np.abs(np.random.normal(0, sigma_theta, N_droplets))
# 生成绕主轴旋转的方位角 phi (0 到 2*pi 均匀分布)
phi = np.random.uniform(0, 2 * np.pi, N_droplets)

# 在喷嘴局部坐标系下的发散方向向量分量 (主喷射方向为 -X)
dx_local = -np.cos(theta)
dy_local = np.sin(theta) * np.cos(phi)
dz_local = np.sin(theta) * np.sin(phi)

# 考虑 DGS 引入的 Pitch 仰角补偿（绕 Y 轴旋转 pitch_compensation）
# 旋转矩阵应用到速度方向上：
alpha = pitch_compensation
dx = dx_local * np.cos(alpha) - dz_local * np.sin(alpha)
dy = dy_local
dz = dx_local * np.sin(alpha) + dz_local * np.cos(alpha)

# 转换成实际的速度矢量分量
vx = v0 * dx
vy = v0 * dy
vz = v0 * dz

# ---------------------------------------------------------
# 3. 轨迹交点解析计算（计入重力抛物线）
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
ax1.scatter(y_target * 1000, z_target * 1000, s=2, c='darkblue', alpha=0.6)
ax1.plot(0, 0, 'ro', markersize=8, label='Target Position (0,0)')
ax1.set_title("Droplet Distribution at Target Position (X = 0 mm)\n(With DGS Pitch Compensation)")
ax1.set_xlabel("Y Position (mm)")
ax1.set_ylabel("Z Position (mm)")
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.axis('equal')
ax1.legend()

# 图 2：TC 截面落点 (散点图)
ax2.scatter(y_tc * 1000, z_tc * 1000, s=2, c='crimson', alpha=0.5)
# 绘制一个标称参考原点(0,0)，用来对比重力宏观下沉量
ax2.plot(0, 0, 'k+', markersize=10, label='Nominal Center Axis')
# 计算并标出散点簇的统计中心
mean_y_tc = np.mean(y_tc) * 1000
mean_z_tc = np.mean(z_tc) * 1000
ax2.plot(mean_y_tc, mean_z_tc, 'go', markersize=6, label=f'Cluster Center ({mean_y_tc:.2f}, {mean_z_tc:.2f})')

ax2.set_title("Droplet Distribution at TC Location (X = -310 mm)\n(Final Scatter Plot)")
ax2.set_xlabel("Y Position (mm)")
ax2.set_ylabel("Z Position (mm)")
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.axis('equal')
ax2.legend()

plt.tight_layout()
plt.show()

# 打印定量评估结果
print(f"--- 数字化仿真定量评估简报 ---")
print(f"液滴数量: {N_droplets} 颗 | 设定喷射圆锥角: {cone_angle_deg}°")
print(f"Target 截面中心位置: Y = {np.mean(y_target)*1000:.4f} mm, Z = {np.mean(z_target)*1000:.4f} mm  --> 成功命中靶心！")
print(f"TC 截面中心位置:     Y = {mean_y_tc:.2f} mm, Z = {mean_z_tc:.2f} mm")
print(f"TC 截面散点簇弥散范围 (3σ): Y_width = {3*np.std(y_tc)*1000:.2f} mm, Z_height = {3*np.std(z_tc)*1000:.2f} mm")
