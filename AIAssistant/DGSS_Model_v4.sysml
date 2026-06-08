package DGSS_Model {

    doc /* Version: 4.1
         * Author: Simon
         * Date: 2026-05-28
         * Change history:
         *   v1.0 (2026-05-28): 初始版本
         *   v2.0 (2026-05-28): 增加 FDS 组件、DGSS 状态机、液滴稳定性定义、
         *                      SafetyConstraint, DGS 调节分辨率/承重/承压/重量需求.
         *   v3.0 (2026-05-28): 集成附录坐标变换 (三系初版, 部分约定为推测).
         *   v4.0 (2026-05-28): 依据完整附录矩阵 (PF/机器/DGS 四系往返) 重写坐标变换,
         *                      并经 verify_transforms.py 数值验证:
         *                      - 显式建模四坐标系 DGSLCS/DGSCS/MCS/PFCS
         *                      - 五环节链: DGSLCS-(命令6DOF)->DGSCS-(β=78°)->MCS-(α=28°)->PFCS
         *                        实时控制运动在 DGSLCS→DGSCS; α/β 为 DGSCS 之上的固定安装关系
         *                      - 用附录精确矩阵替换推测变换 (DGSCS↔MCS, MCS↔PFCS,
         *                        DGSCS↔PFCS 及其逆), 全部互逆/复合一致
         *                      - 修正标定约束为附录 T 分量公式 (L_op/L_os 在机器系,
         *                        L_ps 在 PFCS, 经 Rx(α) 联系而非简单相减)
         *   v4.1 (2026-05-28): 机器坐标系命名由 OriginFrame 改为 MCS (Machine Coordinate
         *                      System); 同步更新枚举值、变换函数名 (DGSCS_to_MCS /
         *                      MCS_to_PFCS) 及全部注释.
         */

    private import SI::*;
    private import ISQ::*;
    private import ScalarValues::*;
    private import SpatialFrames::*;
    private import ScalarFunctions::*;
    private import TrigFunctions::*;
    private import VectorValues::ThreeVectorValue;
    private import NumericalFunctions::*;

    // ============================================================
    //                       基础属性定义
    // ============================================================
    attribute def Position3D {
        attribute x : LengthValue;
        attribute y : LengthValue;
        attribute z : LengthValue;
    }

    attribute def Direction3D {
        attribute dx : Real;
        attribute dy : Real;
        attribute dz : Real;
    }

    attribute def Pose6DOF {
        attribute position : Position3D;
        attribute <Rx> roll  : ISQ::PlaneAngleValue;
        attribute <Ry> pitch : ISQ::PlaneAngleValue;
        attribute <Rz> yaw   : ISQ::PlaneAngleValue;
    }

    attribute def AngularDeviation {
        attribute thetaX : ISQ::PlaneAngleValue;
        attribute thetaY : ISQ::PlaneAngleValue;
        attribute thetaZ : ISQ::PlaneAngleValue;
    }

    // ============================================================
    // ★ v4 重写: 四坐标系 (Reference Frames) 与齐次变换链
    //
    // 经数值验证 (见 verify_transforms.py) 的五环节变换链:
    //
    //   DGSLCS --(命令6DOF)--> DGSCS --(β=78°固定)--> MCS --(α=28°固定)--> PFCS
    //   └ 实时控制运动 ──────┘ └────────── 固定安装关系 (标定) ──────────────┘
    //
    // 四个坐标系:
    //   DGSLCS      — DGS Local Coord. System (下标 l)
    //                 与 DGS 命令 6-DOF 位姿绑定; 命令归零时与 DGSCS 重合.
    //   DGSCS       — DGS Coord. System (机器坐标系中的 "s" 系, 下标 s)
    //                 与 DGS 机械结构刚性绑定.
    //   MCS — 机器坐标系 (Machine Coordinate System, 下标 o)
    //                 DGS 与 PF 共同的安装基准系; L_op, L_os 在此系下表达.
    //   PFCS        — Primary Focus Coord. System (下标 p)
    //                 以 PF 为原点的工艺参考系; 几何/需求评估在此系下表达.
    //
    // 固定安装关系 (来自附录, 已验证):
    //   α = 28°  连接 PFCS ↔ MCS   (绕 X 轴, Rx(α))
    //   β = 78°  连接 MCS ↔ DGSCS  (绕 Y 类轴, 含轴重映射)
    //   L_ps = (916.328834, 43.549721, 81.905114) mm 即 DGSCS→PFCS 复合变换的平移列.
    // ============================================================

    enum def FrameKind {
        DGSLCS;        // DGS Local
        DGSCS;         // DGS (机器系中的 s 系)
        MCS;   // 机器坐标系 (Machine Coordinate System, o 系)
        PFCS;          // Primary Focus
    }

    attribute def Position3DInFrame :> Position3D {
        doc /* 带坐标系标签的 3D 位置, 明确区分同一数值在不同坐标系下的含义. */
        attribute frame : FrameKind;
    }

    attribute def Vector4 {
        doc /* 齐次坐标 (x, y, z, 1), 用于 4×4 矩阵乘法. */
        attribute x : LengthValue;
        attribute y : LengthValue;
        attribute z : LengthValue;
        attribute w : Real = 1.0;
    }

    attribute def HomogeneousTransform {
        doc /* 4×4 齐次变换矩阵, 行优先存储. 底行恒为 (0,0,0,1).
             *   旋转 3×3 子矩阵元素无量纲, 平移列分量为长度. */
        attribute r11 : Real;        attribute r12 : Real;        attribute r13 : Real;        attribute tx : LengthValue;
        attribute r21 : Real;        attribute r22 : Real;        attribute r23 : Real;        attribute ty : LengthValue;
        attribute r31 : Real;        attribute r32 : Real;        attribute r33 : Real;        attribute tz : LengthValue;
    }

    // ---- 通用算子 ----
    calc def applyTransform {
        in T : HomogeneousTransform;
        in p : Position3D;
        return : Position3D;
        // out.x = r11*p.x + r12*p.y + r13*p.z + tx  (其余两行类推)
    }

    calc def composeTransforms {
        in T1 : HomogeneousTransform;     // outer (左乘)
        in T2 : HomogeneousTransform;     // inner (右乘)
        return : HomogeneousTransform;    // = T1 · T2
    }

    calc def invertTransform {
        in T : HomogeneousTransform;
        return : HomogeneousTransform;    // 刚体逆: R^T, -R^T·t
    }

    // ============================================================
    // 实时段: 变换 0 — DGSLCS → DGSCS  (由 DGS 命令 6-DOF 决定)
    //
    //   旋转约定: α_cmd = Rz(yaw), β_cmd = Ry(pitch), γ_cmd = Rx(roll)
    //   R = Rz(α_cmd)·Ry(β_cmd)·Rx(γ_cmd), 平移 = 命令 (tx, ty, tz).
    //   注意: 此处的 β_cmd 是 *实时可变控制量*, 与固定安装角 β=78° 无关.
    // ============================================================
    calc def buildTransform_DGSLCS_to_DGSCS {
        in alpha_Rz : PlaneAngleValue;     // yaw  (实时命令)
        in beta_Ry  : PlaneAngleValue;     // pitch(实时命令)
        in gamma_Rx : PlaneAngleValue;     // roll (实时命令)
        in tx       : LengthValue;
        in ty       : LengthValue;
        in tz       : LengthValue;
        return : HomogeneousTransform;
        // let ca=cos(alpha_Rz), sa=sin(alpha_Rz);
        // let cb=cos(beta_Ry),  sb=sin(beta_Ry);
        // let cg=cos(gamma_Rx), sg=sin(gamma_Rx);
        // r11= ca*cb; r12= ca*sb*sg - sa*cg; r13= ca*sb*cg + sa*sg; tx=tx;
        // r21= sa*cb; r22=-sa*sb*sg + ca*cg; r23=-sa*sb*cg - ca*sg; ty=ty;
        // r31=-sb;    r32= cb*sg;            r33= cb*cg;            tz=tz;
    }

    // ============================================================
    // 固定段①: 变换 1 — DGSCS(s) → MCS(o)   [附录图 "DGS坐标系转换到机器坐标系"]
    //   R = [  0     cos β   sin β ] , 平移 = (Losx, Losy, Losz)
    //       [ -1     0       0     ]
    //       [  0    -sin β   cos β ]
    //   β = 78° (固定安装角)
    // ============================================================
    calc def buildTransform_DGSCS_to_MCS {
        in beta_install : PlaneAngleValue;
        in Losx : LengthValue;
        in Losy : LengthValue;
        in Losz : LengthValue;
        return : HomogeneousTransform;
        // let cb=cos(beta_install), sb=sin(beta_install);
        // r11= 0;  r12= cb; r13= sb; tx=Losx;
        // r21=-1;  r22= 0;  r23= 0;  ty=Losy;
        // r31= 0;  r32=-sb; r33= cb; tz=Losz;
    }

    // ============================================================
    // 固定段②: 变换 2 — MCS(o) → PFCS(p)   [附录图 "机器坐标系转换到PF点坐标系"]
    //   R = [ 1     0        0     ] , 平移 = (-Lopx,
    //       [ 0   cos α    sin α ]              -Lopy·cosα - Lopz·sinα,
    //       [ 0  -sin α    cos α ]               Lopy·sinα - Lopz·cosα )
    //   α = 28° (固定安装角)
    // ============================================================
    calc def buildTransform_MCS_to_PFCS {
        in alpha_install : PlaneAngleValue;
        in Lopx : LengthValue;
        in Lopy : LengthValue;
        in Lopz : LengthValue;
        return : HomogeneousTransform;
        // let ca=cos(alpha_install), sa=sin(alpha_install);
        // r11=1; r12= 0;  r13=0;  tx= -Lopx;
        // r21=0; r22= ca; r23=sa; ty= -Lopy*ca - Lopz*sa;
        // r31=0; r32=-sa; r33=ca; tz=  Lopy*sa - Lopz*ca;
    }

    // ============================================================
    // 固定段 (复合): 变换 3 — DGSCS(s) → PFCS(p)   [附录图 "DGS坐标系转换到PF点坐标系"]
    //   = MCS_to_PFCS · DGSCS_to_MCS (已验证两种算法一致)
    //   R = [   0          cos β         sin β   ] , 平移 = T = (Tx, Ty, Tz)
    //       [ -cos α    -sin α·sin β   sin α·cos β]
    //       [  sin α    -cos α·sin β   cos α·cos β]
    //   Tx = Losx - Lopx
    //   Ty = (Losy - Lopy)·cosα + (Losz - Lopz)·sinα
    //   Tz = (Lopy - Losy)·sinα + (Losz - Lopz)·cosα
    //   验证: T 精确等于 L_ps = (916.328834, 43.549721, 81.905114) mm.
    // ============================================================
    calc def buildTransform_DGSCS_to_PFCS {
        in alpha_install : PlaneAngleValue;
        in beta_install  : PlaneAngleValue;
        in Lopx : LengthValue;   in Lopy : LengthValue;   in Lopz : LengthValue;
        in Losx : LengthValue;   in Losy : LengthValue;   in Losz : LengthValue;
        return : HomogeneousTransform;
        // let ca=cos(alpha_install), sa=sin(alpha_install);
        // let cb=cos(beta_install),  sb=sin(beta_install);
        // r11= 0;   r12= cb;     r13= sb;     tx= Losx - Lopx;
        // r21=-ca;  r22=-sa*sb;  r23= sa*cb;  ty= (Losy-Lopy)*ca + (Losz-Lopz)*sa;
        // r31= sa;  r32=-ca*sb;  r33= ca*cb;  tz= (Lopy-Losy)*sa + (Losz-Lopz)*ca;
    }

    // ============================================================
    // 固定段 (逆): 变换 4 — PFCS(p) → DGSCS(s)   [附录图 "PF点坐标系转换到DGS坐标系"]
    //   = invertTransform(DGSCS_to_PFCS) (已验证 ⑤·⑥ = I)
    //   R = [   0       -cos α        sin α    ] , 平移 = T' = (Tx', Ty', Tz')
    //       [ cos β   -sin α·sin β  -cos α·sin β]
    //       [ sin β    sin α·cos β   cos α·cos β]
    //   Tx' = Losy - Lopy
    //   Ty' = (Lopx - Losx)·cosβ - (Lopz - Losz)·sinβ
    //   Tz' = (Lopx - Losx)·sinβ + (Lopz - Losz)·cosβ
    // ============================================================
    calc def buildTransform_PFCS_to_DGSCS {
        in alpha_install : PlaneAngleValue;
        in beta_install  : PlaneAngleValue;
        in Lopx : LengthValue;   in Lopy : LengthValue;   in Lopz : LengthValue;
        in Losx : LengthValue;   in Losy : LengthValue;   in Losz : LengthValue;
        return : HomogeneousTransform;
        // let ca=cos(alpha_install), sa=sin(alpha_install);
        // let cb=cos(beta_install),  sb=sin(beta_install);
        // r11= 0;   r12=-ca;     r13= sa;      tx= Losy - Lopy;
        // r21= cb;  r22=-sa*sb;  r23=-ca*sb;   ty= (Lopx-Losx)*cb - (Lopz-Losz)*sb;
        // r31= sb;  r32= sa*cb;  r33= ca*cb;   tz= (Lopx-Losx)*sb + (Lopz-Losz)*cb;
    }

    // ============================================================
    // 安装标定 (Calibration) — 全局只读常量, 系统集成阶段标定填入.
    //   L_op, L_os 在 MCS(机器系) 下; L_ps 在 PFCS 下.
    // ============================================================
    attribute def DGSInstallationCalibration {
        doc /* DGS 安装姿态与平移标定. α 连接 PFCS↔MCS, β 连接 MCS↔DGSCS.
             *   均为系统集成期一次性标定的固定值, 运行期视为常量. */
        attribute alpha_install : PlaneAngleValue = 28 ['°'];   // PFCS ↔ MCS
        attribute beta_install  : PlaneAngleValue = 78 ['°'];   // MCS ↔ DGSCS

        // PF 点在机器系 (MCS) 下的位置  L_op
        attribute Lopx : LengthValue =    0.000000 [mm];
        attribute Lopy : LengthValue =  -41.524524 [mm];
        attribute Lopz : LengthValue = 1248.672688 [mm];

        // DGS 原点在机器系 (MCS) 下的位置  L_os
        attribute Losx : LengthValue =  916.328834 [mm];
        attribute Losy : LengthValue =  -41.524524 [mm];
        attribute Losz : LengthValue = 1341.435967 [mm];

        // DGS 原点相对 PF 在 PFCS 下的位置  L_ps (= DGSCS→PFCS 复合平移列 T)
        attribute Lpsx : LengthValue =  916.328834 [mm];
        attribute Lpsy : LengthValue =   43.549721 [mm];
        attribute Lpsz : LengthValue =   81.905114 [mm];

        // ---- 已验证的标定一致性约束 (来自附录 T 分量公式) ----
        //   L_op, L_os 在机器系; L_ps 在 PFCS. 经 Rx(α) 旋转联系, 非简单相减.
        assert constraint translationConsistency {
            Lpsx == (Losx - Lopx) and
            Lpsy == (Losy - Lopy) * cos(alpha_install) + (Losz - Lopz) * sin(alpha_install) and
            Lpsz == (Lopy - Losy) * sin(alpha_install) + (Losz - Lopz) * cos(alpha_install)
        }
    }

    item dgsInstallCalibration : DGSInstallationCalibration;

    // ============================================================
    // 高层封装: DGS 命令位姿 + 安装标定 → 完整 DGSLCS → PFCS 变换 (五环节复合)
    //   T_pl = T_po · T_os · T_sl
    //        = MCS_to_PFCS · DGSCS_to_MCS · DGSLCS_to_DGSCS
    //   或等价地 = DGSCS_to_PFCS · DGSLCS_to_DGSCS  (固定段先合并)
    // ============================================================
    calc def buildTransform_DGSLCS_to_PFCS {
        in alpha_Rz : PlaneAngleValue;   // 实时命令 yaw
        in beta_Ry  : PlaneAngleValue;   // 实时命令 pitch
        in gamma_Rx : PlaneAngleValue;   // 实时命令 roll
        in tx : LengthValue;   in ty : LengthValue;   in tz : LengthValue;
        in calib : DGSInstallationCalibration;
        return : HomogeneousTransform;
        // let T_sl = buildTransform_DGSLCS_to_DGSCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz);
        // let T_ps = buildTransform_DGSCS_to_PFCS(
        //              calib.alpha_install, calib.beta_install,
        //              calib.Lopx, calib.Lopy, calib.Lopz,
        //              calib.Losx, calib.Losy, calib.Losz);
        // return composeTransforms(T_ps, T_sl);
    }

    // ---- 便捷点变换: 把 DGSLCS 下的点直接表达到 PFCS ----
    calc def transformPoint_DGSLCS_to_PFCS {
        in pointInDGSLCS : Position3D;
        in alpha_Rz : PlaneAngleValue;
        in beta_Ry  : PlaneAngleValue;
        in gamma_Rx : PlaneAngleValue;
        in tx : LengthValue;   in ty : LengthValue;   in tz : LengthValue;
        in calib : DGSInstallationCalibration;
        return : Position3D;
        // let T = buildTransform_DGSLCS_to_PFCS(alpha_Rz,beta_Ry,gamma_Rx,tx,ty,tz,calib);
        // return applyTransform(T, pointInDGSLCS);
    }

    attribute def EnvironmentCondition {
        doc /* 飞行环境条件，影响液滴轨迹 */
        attribute gravity : AccelerationValue = 9.80665 ['m/s²'];
        attribute gravityDirection : Direction3D {
            :>> dx =  0.0;
            :>> dy =  0.0;
            :>> dz = -1.0;
        }
        attribute ambientPressure        : PressureValue;
        attribute meanGasVelocity        : ThreeVectorValue;
        attribute gasVelocityRMS         : SpeedValue;
        attribute gasVelocityCorrTime    : TimeValue;
    }

    item ambientEnvironment : EnvironmentCondition {
        :>> ambientPressure    = 1.0e-3 [Pa];   // 典型 EUV 舱体真空度
    }

    // ============================================================
    //                       空间参考点定义
    // ============================================================
    item def SpatialReferencePoint {
        attribute position : Position3D;
        attribute frame    : FrameKind = FrameKind::PFCS;   // 默认在 PFCS 下定义
    }

    item <PF> 'Primary Focus' : SpatialReferencePoint {
        :>> position {
            :>> x = 0 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> frame = FrameKind::PFCS;   // PF 是 PFCS 的原点
    }

    item monitorPosition : SpatialReferencePoint {
        :>> position {
            :>> x = 5 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> frame = FrameKind::PFCS;
    }

    item nominalNozzlePosition : SpatialReferencePoint {
        :>> position {
            :>> x = 410.453510 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> frame = FrameKind::PFCS;
    }

    // ============================================================
    //                  ★ 新增: 液滴稳定性定义
    // ============================================================
    attribute def DropletStability {
        doc /* 液滴稳定性定义:
             *   - 评估窗口: 任意连续 30 s
             *   - 评估指标: 在指定参考平面 (默认: FDS 监测平面 x = monitorPosition.x)
             *               处, 液滴穿越点在 Y / Z / X 方向位置时间序列的 3σ 值.
             *   - 合格判据: 各方向 3σ 不稳定性 < 15 μm.
             *   - X 方向 (纵向): 反映液滴生成时序/速度的稳定性, 可由液滴到达
             *                    监测平面的时间抖动 (Δt) 乘以名义速度 |V0| 估算.
             *   - Y / Z 方向 (横向): 反映液滴落点空间散布. */

        attribute evaluationWindow      : TimeValue   = 30 [s];
        attribute sigmaMultiplier       : Real        = 3.0;
        attribute stabilityThreshold    : LengthValue = 15 [um];

        // 当前测得的 3σ 不稳定性 (由 FDS 实时计算)
        attribute longitudinalSigma3_X  : LengthValue;   // 纵向 (X 方向)
        attribute transverseSigma3_Y    : LengthValue;   // 横向 Y
        attribute transverseSigma3_Z    : LengthValue;   // 横向 Z

        // 派生布尔标志: 用于状态机转移条件
        attribute isLongitudinalStable  : Boolean;       // X 方向稳定性达标
        attribute isTransverseStable    : Boolean;       // Y 与 Z 方向均达标
        attribute isFullyStable         : Boolean;       // X / Y / Z 全部达标

        // ---- 派生标志计算约束 ----
        assert constraint deriveLongitudinalStable {
            isLongitudinalStable == (longitudinalSigma3_X < stabilityThreshold)
        }
        assert constraint deriveTransverseStable {
            isTransverseStable == (
                transverseSigma3_Y < stabilityThreshold and
                transverseSigma3_Z < stabilityThreshold
            )
        }
        assert constraint deriveFullyStable {
            isFullyStable == (isLongitudinalStable and isTransverseStable)
        }
    }

    // ============================================================
    //                          Nozzle
    // ============================================================
    part def Nozzle {
        attribute currentPose : Pose6DOF;
        attribute nominalDirection : Direction3D {
            :>> dx = -1.0;
            :>> dy =  0.0;
            :>> dz =  0.0;
        }
        attribute actualDirection     : Direction3D;
        attribute directionDeviation  : AngularDeviation;
        attribute positionDeviation   : Position3D;

        // ---- ★ v3 新增: 喷嘴在 DGSLCS 下的标称几何 ----
        //   (相对于 DGS 本地系的固定安装关系, 与 DGS 命令无关)
        attribute outletPositionInDGSLCS  : Position3D;
        attribute outletDirectionInDGSLCS : Direction3D;

        // 喷嘴出口在 PFCS 下的实时位置 (由 DGS 当前变换映射)
        attribute outletPositionInPFCS    : Position3D;

        port outletPort   : DropletOutletPort;
        port mountingPort : MechanicalMountPort;
    }

    // ============================================================
    //                          DG / DGS
    // ============================================================
    part def <DG> 'Droplet Generator'{
        port nozzleMount : MechanicalMountPort;
        port stageMount  : MechanicalMountPort;
        attribute weight :>> ISQ::mass;
    }

    part def <DGS> 'Droplet Generator Stage' {
        attribute commandedTx     : LengthValue;
        attribute commandedTy     : LengthValue;
        attribute commandedTz     : LengthValue;
        attribute <commandedRx> commandedRoll  : PlaneAngleValue;
        attribute <commandedRy> commandedPitch : PlaneAngleValue;
        attribute <commandedRz> commandedYaw   : PlaneAngleValue;

        attribute currentPose : Pose6DOF;

        attribute txMin : LengthValue     = -9   [mm];   attribute txMax : LengthValue     = +9   [mm];
        attribute tyMin : LengthValue     = -10  [mm];   attribute tyMax : LengthValue     = +10  [mm];
        attribute tzMin : LengthValue     = -4.5 [mm];   attribute tzMax : LengthValue     = +4.5 [mm];
        attribute rxMin : PlaneAngleValue = -2.1 ['°'];  attribute rxMax : PlaneAngleValue = +2.1 ['°'];
        attribute ryMin : PlaneAngleValue = -2.4 ['°'];  attribute ryMax : PlaneAngleValue = +2.4 ['°'];
        attribute rzMin : PlaneAngleValue = -4.3 ['°'];  attribute rzMax : PlaneAngleValue = +4.3 ['°'];

        attribute axesAreCoupled : Boolean = true;

        attribute requiredAdjustableRegionCenter   : String      = "PF";
        attribute requiredAdjustableRegionPlane    : String      = "YZ-plane through PF";
        attribute requiredAdjustableRegionShape    : String      = "circle";
        attribute requiredAdjustableRegionDiameter : LengthValue = 15 [mm];

        attribute achievableAdjustableRegionCenter   : String      = "PF";
        attribute achievableAdjustableRegionPlane    : String      = "YZ-plane through PF";
        attribute achievableAdjustableRegionShape    : String      = "circle";
        attribute achievableAdjustableRegionDiameter : LengthValue;

        // ---- 闭环带宽 + 测试规范 ----
        attribute closedLoopBandwidth          : FrequencyValue;
        attribute bandwidthDefinitionCriterion : String = "-3dB magnitude attenuation";
        attribute bandwidthTestMethod          : String = "Frequency-response identification via swept-sine on each DOF";
        attribute bandwidthTestExcitationType  : String = "swept sine | chirp | white noise";
        attribute bandwidthTestAmplitude       : LengthValue;
        attribute bandwidthTestAmplitudeAng    : PlaneAngleValue;
        attribute bandwidthTestFreqRange_low   : FrequencyValue = 0.1 [Hz];
        attribute bandwidthTestFreqRange_high  : FrequencyValue = 100 [Hz];
        attribute bandwidthTestSensor          : String = "Internal capacitive sensor or external interferometer";

        // ---- ★ 新增: 调节分辨率属性 ----
        attribute translationResolution        : LengthValue;        // 平移轴最小可分辨步长
        attribute rotationResolution           : PlaneAngleValue;    // 旋转轴最小可分辨步长

        // ---- ★ 新增: 承重 / 承压 / 自身重量属性 ----
        attribute loadCapacity                 : MassValue;          // 可承载负载质量 (kg)
        attribute workingPressureCapacity      : PressureValue;      // 工作中可承受的环境压力上限 (Pa)
        attribute stageOwnWeight               : MassValue;          // DGS 自身质量 (kg)

        // ---- ★ v3 新增: 坐标系标定与当前变换 ----
        ref item installCalib = dgsInstallCalibration;

        // 当前 DGSLCS → DGSCS 变换 (由命令位姿派生)
        attribute currentTransform_DGSLCS_to_DGSCS : HomogeneousTransform;
        // 当前 DGSLCS → PFCS  变换 (复合)
        attribute currentTransform_DGSLCS_to_PFCS  : HomogeneousTransform;

        // 变换一致性约束: 当前变换必须与命令位姿吻合
        assert constraint transformConsistency_sl {
            currentTransform_DGSLCS_to_DGSCS ==
                buildTransform_DGSLCS_to_DGSCS(
                    commandedYaw, commandedPitch, commandedRoll,
                    commandedTx, commandedTy, commandedTz
                )
        }
        assert constraint transformConsistency_pl {
            currentTransform_DGSLCS_to_PFCS ==
                buildTransform_DGSLCS_to_PFCS(
                    commandedYaw, commandedPitch, commandedRoll,
                    commandedTx, commandedTy, commandedTz,
                    installCalib
                )
        }

        port controlPort : ControlSignalPort;
        port loadPort    : MechanicalMountPort;

        // ============================================================
        // 运动学可达性占位声明
        // ============================================================
        calc def reachableWorkspace {
            in tx : LengthValue;
            in ty : LengthValue;
            in tz : LengthValue;
            in rx : PlaneAngleValue;
            in ry : PlaneAngleValue;
            in rz : PlaneAngleValue;
            return : Boolean;
        }

        // ---- 单轴行程极限 ----
        assert constraint singleAxisTravelLimits {
            commandedTx    >= txMin and commandedTx    <= txMax and
            commandedTy    >= tyMin and commandedTy    <= tyMax and
            commandedTz    >= tzMin and commandedTz    <= tzMax and
            commandedRoll  >= rxMin and commandedRoll  <= rxMax and
            commandedPitch >= ryMin and commandedPitch <= ryMax and
            commandedYaw   >= rzMin and commandedYaw   <= rzMax
        }

        // ---- 耦合工作空间可达性 ----
        assert constraint coupledWorkspaceReachability {
            reachableWorkspace(
                commandedTx, commandedTy, commandedTz,
                commandedRoll, commandedPitch, commandedYaw
            ) == true
        }

        // ============================================================
        // Req-1: 可调区域覆盖
        // ============================================================
        requirement <'Req-1'> dropletAdjustableRangeRequirement {
            doc /* 在过 PF 点的 YZ 平面上, 以 PF 点为中心, 液滴落点应在
                 * 直径不小于 15 mm 的圆内任意可调. */
            assume constraint regionGeometryMatch {
                achievableAdjustableRegionCenter == requiredAdjustableRegionCenter and
                achievableAdjustableRegionPlane  == requiredAdjustableRegionPlane  and
                achievableAdjustableRegionShape  == requiredAdjustableRegionShape
            }
            require constraint diameterCoverage {
                achievableAdjustableRegionDiameter >= requiredAdjustableRegionDiameter
            }
        }

        // ============================================================
        // Req-2: 闭环带宽
        // ============================================================
        requirement <'Req-2'> closedLoopBandwidthRequirement {
            doc /* DGS 各受控自由度闭环控制带宽 ≥ 5 Hz, 以 -3dB 幅值衰减为判定准则. */
            require constraint bandwidthLowerBound {
                closedLoopBandwidth >= 5 [Hz]
            }
        }

        // ============================================================
        // ★ Req-3: DGS 调节分辨率需求
        // ============================================================
        requirement <'Req-3'> adjustmentResolutionRequirement {
            doc /* DGS 平移自由度 (Tx, Ty, Tz) 的调节分辨率应不大于 1 μm,
                 *   即最小可控/可分辨步长 ≤ 1 μm.
                 *   测试方法: 在闭环控制下连续发送小步长指令, 用外部高分辨率
                 *   位移传感器 (干涉仪) 监测实际位移, 验证步长可重复分辨. */
            attribute requiredTranslationResolution : LengthValue = 1 [um];
            require constraint translationResolutionBound {
                translationResolution <= requiredTranslationResolution
            }
        }

        // ============================================================
        // ★ Req-4: DGS 承重能力需求
        // ============================================================
        requirement <'Req-4'> loadCapacityRequirement {
            doc /* DGS 承重能力应不小于 25 kg, 用于支撑 Droplet Generator 本体
                 *   及附属管路/电缆. */
            attribute requiredLoadCapacity : MassValue = 25 [kg];
            require constraint loadCapacityBound {
                loadCapacity >= requiredLoadCapacity
            }
        }

        // ============================================================
        // ★ Req-5: DGS 工作中承压能力需求
        // ============================================================
        requirement <'Req-5'> workingPressureRequirement {
            doc /* DGS 工作中的环境承压能力应不小于 1×10^5 Pa
                 *   (约一个大气压, 覆盖从真空到大气环境的全工况). */
            attribute requiredWorkingPressure : PressureValue = 1.0e5 [Pa];
            require constraint pressureBound {
                workingPressureCapacity >= requiredWorkingPressure
            }
        }

        // ============================================================
        // ★ Req-6: DGS 自身重量约束
        // ============================================================
        requirement <'Req-6'> stageWeightRequirement {
            doc /* DGS 自身重量应小于 40 kg, 以满足上级结构承载与装配限制. */
            attribute maxStageWeight : MassValue = 40 [kg];
            require constraint stageWeightBound {
                stageOwnWeight < maxStageWeight
            }
        }
    }

    // ============================================================
    //                  闭环带宽测试用例
    // ============================================================
    use case def MeasureClosedLoopBandwidth {
        subject dgs : DGS;
        objective {
            doc /* 在每个受控自由度上, 通过阶跃响应法测量上升时间,
                 * 并按一阶系统近似公式 BW ≈ 0.35 / tr 估算闭环带宽. */
        }
        first start;
        then action setupSensor;
        then action setOperatingPoint;
        then action injectStepCommand;
        then action recordStepResponse;
        then action extractRiseTime;
        then action estimateBandwidth;
        then action verifyAgainstRequirement;
    }

    // ============================================================
    //                  弹道轨迹模型
    // ============================================================
    attribute def BallisticTrajectory {
        doc /* 考虑重力与气流扰动的弹道轨迹 */
        attribute startPoint     : Position3D;
        attribute initialVelocity: ThreeVectorValue;
        attribute env            : EnvironmentCondition;

        attribute nominalCrossingAtPF      : Position3D;
        attribute nominalCrossingAtMonitor : Position3D;
        attribute nominalCrossingAtShroud  : Position3D;
        attribute nominalCrossingAtHole    : Position3D;
        attribute nominalCrossingAtTC      : Position3D;

        attribute crossingDispersionAtTC_y : LengthValue;
        attribute crossingDispersionAtTC_z : LengthValue;
    }

    item def Droplet {
        attribute generationPosition : Position3D;
        attribute velocityDirection  : Direction3D;
        attribute initialSpeed       : SpeedValue;
        attribute mass               : MassValue;
        attribute diameter           : LengthValue;
        attribute trajectory         : BallisticTrajectory;
        attribute stability          : DropletStability;   // ★ 新增: 液滴稳定性
    }

    // ============================================================
    //                  通孔 / TC 开口几何
    // ============================================================
    attribute def ApertureGeometry {
        attribute shape          : String;
        attribute centerPosition : Position3D;
        attribute width          : LengthValue;
        attribute height         : LengthValue;
    }

    attribute def TCApertureGeometry {
        doc /* TC 开口: 上半圆 (R=33mm) + 下矩形 (66 × 20.8873 mm) */
        attribute openingPlane          : String = "YZ-plane at TC";
        attribute centerPosition        : Position3D{
            :>> x = -310.0930 [mm];
            :>> y = 0 [mm];
            :>> z = 5.9873 [mm];
        }
        attribute upperSemicircleRadius : LengthValue = 33      [mm];
        attribute lowerRectangleWidth   : LengthValue = 66      [mm];
        attribute lowerRectangleHeight  : LengthValue = 20.8873 [mm];

        attribute idealHitOffsetFromBottom : LengthValue = 14.9 [mm];

        calc def getIdealHitLocalZ {
            in geo : TCApertureGeometry;
            return : LengthValue;
        }
    }

    calc def pointInsideTCAperture {
        in yp  : LengthValue;
        in zp  : LengthValue;
        in geo : TCApertureGeometry;
        return : Boolean;
    }

    calc def pointInsideAperture {
        in yp  : LengthValue;
        in zp  : LengthValue;
        in geo : ApertureGeometry;
        return : Boolean;
    }

    // ============================================================
    //                  Shroud / InletHole / TC
    // ============================================================
    part def Shroud {
        attribute mountingPose : Pose6DOF;
        attribute aperture     : ApertureGeometry {
            :>> shape = "circle";
        }
        port mountingPort : MechanicalMountPort;
    }

    part def InletHole {
        attribute aperture : ApertureGeometry {
            :>> shape = "circle";
        }
    }

    part def <TC> 'Tin Capture' {
        attribute mountingPose           : Pose6DOF;
        attribute aperture               : TCApertureGeometry;
        attribute idealCollectionPosition: Position3D;
        port inletPort                   : DropletInletPort;
        attribute hasCollected           : Boolean = false;
    }

    // ============================================================
    //                  ★ 新增: FDS (小视场组件)
    // ============================================================
    part def <FDS> 'Fine Detection System' {
        doc /* 小视场探测系统, 用于实时监测液滴在 monitorPosition 附近的
             *   位置/到达时序, 并为 DGS 闭环锁定提供反馈. */

        // ---- 子组件 ----
        part def IlluminationComponent {
            doc /* 照明子组件: 提供照亮液滴所需的光源 (例如脉冲激光),
                 *   与探测器同步触发, 形成短曝光成像. */
            attribute wavelength      : LengthValue;
            attribute pulseDuration   : TimeValue;
            attribute repetitionRate  : FrequencyValue;
            attribute opticalPower    : PowerValue;
            port lightOutPort : MechanicalMountPort;
        }

        part def DetectionComponent {
            doc /* 探测子组件: 高速成像/光电传感器, 测量液滴穿越参考平面
                 *   时的横向位置 (Y, Z) 与到达时间 (用于 X 方向稳定性估计). */
            attribute frameRate          : FrequencyValue;
            attribute spatialResolution  : LengthValue;
            attribute fieldOfView_Y      : LengthValue;
            attribute fieldOfView_Z      : LengthValue;
            attribute timingResolution   : TimeValue;
            port imagingInPort  : MechanicalMountPort;
            port dataOutPort    : MechanicalMountPort;
        }

        part illumination : IlluminationComponent;
        part detection    : DetectionComponent;

        // FDS 整体参考监测位置
        ref item monitorRef = monitorPosition;

        // 由 FDS 实时输出的液滴稳定性测量结果
        attribute measuredStability : DropletStability;

        port feedbackPort : ControlSignalPort;

        // ---- FDS 性能需求 ----
        requirement <'Req-FDS-1'> stabilityMeasurementAccuracy {
            doc /* FDS 在 Y / Z 方向的空间测量分辨率应不大于 1 μm,
                 *   以保证可分辨 15 μm/3σ 的稳定性阈值. */
            require constraint spatialResolutionBound {
                detection.spatialResolution <= 1 [um]
            }
        }
    }

    // ============================================================
    //                          端口定义
    // ============================================================
    port def DropletOutletPort {
        out dropletFlow : Droplet;
    }

    port def DropletInletPort {
        in dropletFlow : Droplet;
    }

    port def MechanicalMountPort { }

    port def ControlSignalPort {
        in command : DGSCommand;
    }

    item def DGSCommand {
        attribute tx    : LengthValue;
        attribute ty    : LengthValue;
        attribute tz    : LengthValue;
        attribute roll  : ScalarValues::Real;
        attribute pitch : ScalarValues::Real;
        attribute yaw   : ScalarValues::Real;
    }

    // ============================================================
    //                          DGSS 系统
    // ============================================================
    part def DGSS {
        part dgs    : DGS;
        part body   : DG;
        part nozzle : Nozzle;

        part shroud : Shroud;
        part inlet  : InletHole;

        part tc : TC {
            :>> idealCollectionPosition {
                attribute :>> x = -310.0930 [mm];
                attribute :>> y = 0 [mm];
                attribute :>> z = 0 [mm];
            }
        }

        // ★ 新增: FDS 作为 DGSS 的组成部件
        part fds : FDS;

        // ★ v3 新增: 喷嘴出口在 PFCS 下的实时位置 = T_pl · 喷嘴出口在 DGSLCS 下的位置
        assert constraint nozzleOutletInPFCSConsistency {
            nozzle.outletPositionInPFCS ==
                applyTransform(
                    dgs.currentTransform_DGSLCS_to_PFCS,
                    nozzle.outletPositionInDGSLCS
                )
        }

        ref item targetRef  = PF;
        ref item monitorRef = monitorPosition;
        ref item envRef     = ambientEnvironment;

        // ---- 当前运行状态 (供外部观察) ----
        attribute currentState : String;        // "Initialization" | "Idle" | "ClosedLoopLocking" | "Fault"
        attribute faultFlag    : Boolean = false;

        // ============================================================
        // ★ 新增: DGSS 状态机
        //   States: Initialization → Idle → ClosedLoopLocking
        //   Any state → Fault (异常发生)
        //   Fault → Initialization (复位后重新初始化)
        // ============================================================
        state DGSSStateMachine {
            entry; then Initialization;

            // ---------------- Initialization ----------------
            state Initialization {
                entry action homeDGS;
                    doc /* DGS 回零: 所有平移/旋转轴执行 homing 序列,
                         *   确认编码器/限位开关状态. */
                then action moveDropletToFDSRef;
                    doc /* 控制 DGS 使液滴落点对准 FDS 参考位置 (monitorPosition),
                         *   完成粗对准. */
                then action verifyInitComplete;
                    doc /* 校验所有轴位置在期望区间内, 无故障标志. */
            }
            // 初始化完成 → Idle
            transition initToIdle
                first Initialization
                accept after (1 [s])      // 占位: 实际由 verifyInitComplete 完成事件触发
                if (not faultFlag)
                then Idle;
            // 初始化失败 → Fault
            transition initToFault
                first Initialization
                if (faultFlag)
                then Fault;

            // ---------------- Idle ----------------
            state Idle {
                doc /* 初始化完成或从 ClosedLoopLocking 退出后的待命状态.
                     *   液滴可能正在喷射, 但 DGS 控制环未参与闭环锁定. */
                entry action openControlLoop;
                    doc /* 释放闭环锁定 (若曾激活), DGS 维持当前姿态. */
                do action monitorStability;
                    doc /* 持续从 FDS 读取液滴稳定性测量结果, 但不闭环修正. */
            }
            // Idle → ClosedLoopLocking : 仅当 Y/Z 横向稳定性达标
            transition idleToLocking
                first Idle
                if (fds.measuredStability.isTransverseStable and not faultFlag)
                then ClosedLoopLocking;
            // Idle → Fault
            transition idleToFault
                first Idle
                if (faultFlag)
                then Fault;

            // ---------------- ClosedLoopLocking ----------------
            state ClosedLoopLocking {
                doc /* 利用 FDS 实时反馈, 闭环调整 DGS, 使液滴稳定锁定在
                     *   FDS 参考位置 (monitorPosition). */
                entry action engageClosedLoop;
                    doc /* 接入 FDS → DGS 闭环路径, 开启位置闭环控制器. */
                do action trackReference;
                    doc /* 实时计算液滴位置误差, 输出 DGS 修正指令. */
            }
            // ClosedLoopLocking → Idle : 横向稳定性失守 或 收到退出指令
            transition lockingToIdle
                first ClosedLoopLocking
                if (
                    (not fds.measuredStability.isTransverseStable) and
                    (not faultFlag)
                )
                then Idle;
            // ClosedLoopLocking → Fault
            transition lockingToFault
                first ClosedLoopLocking
                if (faultFlag)
                then Fault;

            // ---------------- Fault ----------------
            state Fault {
                doc /* 故障状态: 检测到不可恢复异常 (运动学超限/传感器丢失/
                     *   液滴丢失/TC 收集失败超阈值/承压超限 等).
                     *   触发安全响应序列 (见 SafetyConstraint). */
                entry action raiseAlarm;
                entry action executeSafetyResponse;
                    doc /* 见 dropletCannotBeCollectedSafetyResponse 定义的应对措施. */
            }
            // Fault → Initialization : 外部复位
            transition faultToInit
                first Fault
                accept after (1 [s])      // 占位: 实际由复位指令触发
                if (not faultFlag)
                then Initialization;
        }
    }

    // ============================================================
    //                  弹道计算函数
    // ============================================================
    calc def computeBallisticTrajectory {
        in startPoint      : Position3D;
        in initialVelocity : ThreeVectorValue;
        in env             : EnvironmentCondition;
        return             : BallisticTrajectory;
    }

    calc def crossingPointAtX {
        in trajectory : BallisticTrajectory;
        in xPlane     : LengthValue;
        return        : Position3D;
    }

    // ============================================================
    //                          动作定义
    // ============================================================
    action def GenerateDroplet {
        in nozzle : Nozzle;
        out droplet : Droplet {
            attribute :>> generationPosition = nozzle.currentPose.position;
            attribute :>> velocityDirection  = nozzle.actualDirection;
        }
    }

    action def PropagateDroplet {
        in  droplet  : Droplet;
        in  target   : SpatialReferencePoint;
        in  monitor  : SpatialReferencePoint;
        out passesTarget  : Boolean;
        out passesMonitor : Boolean;
    }

    action def CollectDroplet {
        in  droplet   : Droplet;
        in  collector : TC;
        out collected : Boolean;
    }

    calc def distancePointToLine {
        in point : Position3D;
        in line  : BallisticTrajectory;
        return   : LengthValue;
    }

    // ============================================================
    //                  轨迹通过性需求
    // ============================================================
    requirement def DropletPassesTarget {
        subject d : Droplet;
        ref item target = PF;
        require constraint {
            distancePointToLine(target.position, d.trajectory) <= 0.001[mm]
        }
    }

    requirement def DropletPassesMonitor {
        subject d : Droplet;
        ref item monitor = monitorPosition;
        require constraint {
            distancePointToLine(monitor.position, d.trajectory) <= 0.001[mm]
        }
    }

    requirement def DropletCollectedByTC {
        subject d  : Droplet;
        ref part c : TC;
        require constraint {
            pointInsideTCAperture(
                crossingPointAtX(d.trajectory, c.aperture.centerPosition.x).y
                    - c.aperture.centerPosition.y,
                crossingPointAtX(d.trajectory, c.aperture.centerPosition.x).z
                    - c.aperture.centerPosition.z,
                c.aperture
            ) == true
        }
    }

    // ============================================================
    //                  ★ 新增: 液滴稳定性需求
    // ============================================================
    requirement def DropletStabilityRequirement {
        doc /* 液滴在喷射过程中, 任意连续 30 s 窗口内,
             *   X / Y / Z 三方向的 3σ 不稳定性均应 < 15 μm.
             *   测量方法: 由 FDS 在 monitorPosition 平面处实时统计. */
        subject d : Droplet;
        require constraint stabilityBound {
            d.stability.longitudinalSigma3_X < d.stability.stabilityThreshold and
            d.stability.transverseSigma3_Y  < d.stability.stabilityThreshold and
            d.stability.transverseSigma3_Z  < d.stability.stabilityThreshold
        }
    }

    // ============================================================
    //                  目标层次结构
    // ============================================================
    requirement def IdealGoal {
        subject d : Droplet;
        requirement passMonitor : DropletPassesMonitor;
        requirement passTarget  : DropletPassesTarget;
        requirement collected   : DropletCollectedByTC;
        requirement stable      : DropletStabilityRequirement;  // ★ 新增

        bind passMonitor.d = d;
        bind passTarget.d  = d;
        bind collected.d   = d;
        bind stable.d      = d;
    }

    requirement def FallbackGoal {
        subject d : Droplet;
        requirement passTarget : DropletPassesTarget;
        requirement collected  : DropletCollectedByTC;

        bind passTarget.d  = d;
        bind collected.d   = d;
    }

    requirement goalHierarchy {
        requirement ideal : IdealGoal;
    }

    // ============================================================
    //                  弹道预算需求
    // ============================================================
    requirement def DropletBallisticBudget {
        subject d  : Droplet;
        ref part tcRef : TC;

        require constraint nominalDropBudget {
            NumericalFunctions::abs(d.trajectory.nominalCrossingAtTC.z - tcRef.aperture.centerPosition.z)
                <= tcRef.aperture.lowerRectangleHeight / 2
        }
        require constraint dispersionBudget {
            3 * d.trajectory.crossingDispersionAtTC_y
                <= tcRef.aperture.lowerRectangleWidth / 6  and
            3 * d.trajectory.crossingDispersionAtTC_z
                <= tcRef.aperture.lowerRectangleHeight / 6
        }
    }

    // ============================================================
    // ★ 新增: SafetyConstraint - 液滴无法被 TC 收集时的应对措施
    // ============================================================
    requirement def SafetyConstraint {
        doc /* 安全约束基类: 定义不安全场景的检测与应对响应. */
        attribute hazardDescription   : String;
        attribute detectionMethod     : String;
        attribute responseAction      : String;
        attribute severityLevel       : String;     // "Low" | "Medium" | "High" | "Critical"
    }

    requirement <'Safety-1'> dropletCannotBeCollectedSafetyResponse :> SafetyConstraint {
        doc /* 液滴无法被 TC 收集时的安全约束.
             *
             * 触发条件 (任一成立):
             *   (a) 液滴名义穿越点落在 TC 开口外 (pointInsideTCAperture == false)
             *   (b) 连续 N 个液滴 (建议 N ≥ 10) 检测到 TC 收集失败
             *   (c) FDS 丢失液滴信号且重试无效
             *   (d) DGS 控制误差超出可恢复范围
             *
             * 应对措施 (按优先级顺序执行):
             *   1) 立即关闭液滴生成器 (停止液滴喷射), 防止未被收集的液滴
             *      污染收集镜 (Collector Mirror) 与光学元件;
             *   2) 切换 DGSS 状态机至 Fault 状态, 释放 ClosedLoopLocking;
             *   3) 关闭/隔离驱动激光 (Drive Laser) 触发, 避免无效辐射;
             *   4) 触发 Shroud / Inlet Hole 区域的局部抽气, 抽除残留锡蒸汽;
             *   5) 向上位机/操作员报警, 报告故障类型与最近一次液滴轨迹偏差;
             *   6) 记录故障时刻完整的遥测数据 (DGS 姿态、FDS 测量、环境压力)
             *      供后续根因分析;
             *   7) 在确认 TC 入口与收集镜内舱清洁后, 由操作员手动触发复位,
             *      进入 Initialization 重新初始化. */
        :>> hazardDescription =
            "Droplet trajectory does not enter TC aperture, risking tin contamination of collector optics.";
        :>> detectionMethod   =
            "Geometric check via pointInsideTCAperture + statistical TC-miss counter + FDS loss-of-droplet detection.";
        :>> responseAction    =
            "Stop droplet generator → enter Fault → inhibit drive laser → local pump-down → alarm → log telemetry → manual reset to Initialization.";
        :>> severityLevel     = "Critical";

        subject system : DGSS;

        // 安全约束: 在检测到不可收集场景时, DGSS 必须进入 Fault 状态
        require constraint enterFaultOnUncollectable {
            system.currentState == "Fault"  or
            system.faultFlag    == false
            // 语义: 若 faultFlag 已经被置位 (检测到不可收集),
            //       则系统必须处于 Fault 状态.
        }
    }

    // ============================================================
    //                  用例: 调整喷嘴命中目标
    // ============================================================
    use case def AdjustNozzleToHitTarget {
        subject system : DGSS;
        objective {
            doc /* 通过控制 DGS, 使液滴满足 IdealGoal; 否则至少满足 FallbackGoal. */
        }
        first start;
        then action measureDeviation;
        then action computeDGSCommand;
        then action applyDGSCommand;
        then action emitDroplet;
        then action evaluateGoal;
        then done;
    }

    // ============================================================
    // ★ v3 新增: 坐标变换用例 — DGS 命令位姿求解
    //
    // 给定 PF 系下的目标落点 (例如 monitorPosition 或 PF 本身),
    // 求解使液滴出口位于该目标处所需的 DGS 命令位姿.
    // 该用例描述系统集成期的标定流程与运行期的反解流程.
    // ============================================================
    use case def SolveDGSCommandFromPFTarget {
        subject system : DGSS;

        objective {
            doc /* 由 PFCS 中给定目标点反解 DGS 6-DOF 命令.
                 *
                 *   正向链:
                 *     P_p = T_ps(α=28°, β=78°, L_ps) · T_sl(α_cmd, β_cmd, γ_cmd; tx, ty, tz) · P_l
                 *
                 *   反解流程 (利用 invertTransform):
                 *     1) 已知 P_p (PF 系目标) 与 P_l (喷嘴出口在 DGSLCS 下的固定位置)
                 *     2) 求 T_sl 使得  T_ps · T_sl · P_l = P_p
                 *     3) 从 T_sl 提取 (α, β, γ, tx, ty, tz) 作为 DGS 命令
                 *     4) 验证落在 DGS 行程极限内 (singleAxisTravelLimits) 且
                 *        在耦合工作空间内 (coupledWorkspaceReachability)
                 */
        }

        first start;

        then action loadCalibration;
            doc /* 加载 dgsInstallCalibration: α=28°, β=78°,
                 *   Lpsx=916.328834 mm, Lpsy=43.549721 mm, Lpsz=81.905114 mm. */

        then action specifyTargetInPFCS;
            doc /* 指定 PF 系下的目标点 P_p. 典型取值: PF (0,0,0) 或 monitorPosition (5,0,0). */

        then action computeT_ps;
            doc /* 由标定参数构建 DGSCS → PFCS 固定变换 T_ps. */

        then action computeRequiredT_sl;
            doc /* 求解 T_sl: 利用 invertTransform(T_ps) 得 T_sp,
                 *   并由 T_sp · P_p 给出喷嘴出口在 DGSCS 下应到达的位置 P_s.
                 *   再结合 P_l (喷嘴出口在 DGSLCS 下的标称位置) 解 T_sl. */

        then action extractEulerAndTranslation;
            doc /* 从 T_sl 的旋转 3×3 子矩阵按 ZYX 反解 (α, β, γ),
                 *   并从第 4 列得 (tx, ty, tz). */

        then action verifyWithinTravelLimits;
            doc /* 校验 (tx, ty, tz, α, β, γ) 满足 DGS 的 singleAxisTravelLimits
                 *   与 coupledWorkspaceReachability. 若不满足则报告反解失败. */

        then action issueDGSCommand;
            doc /* 把反解结果打包为 DGSCommand 通过 controlPort 下发. */

        then done;
    }
}
