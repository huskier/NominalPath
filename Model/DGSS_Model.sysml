package DGSS_Model {

    doc /* Version: 3.0
         * Author: Simon
         * Date: 2026-05-28
         * Change history:
         *   v1.0 (2026-05-28): 初始版本
         *   v2.0 (2026-05-28): 增加 FDS 组件、DGSS 状态机、液滴稳定性定义、
         *                      SafetyConstraint, DGS 调节分辨率/承重/承压/重量需求.
         *   v3.0 (2026-05-28): 集成附录坐标变换:
         *                      - 显式建模 DGSLCS / DGSCS / PFCS 三个坐标系
         *                      - 定义 4×4 齐次变换矩阵与变换函数
         *                      - DGSLCS → DGSCS: 由 DGS 命令位姿 (6-DOF, ZYX 欧拉角) 决定
         *                      - DGSCS → PFCS:   固定安装变换 (α=28°, β=78°, 标定平移)
         *                      - 提供 Position3D 在坐标系间互转的 calc def
         *                      - 把命中目标的需求改写为 PFCS 表达
         */

    private import SI::*;
    private import ISQ::*;
    private import ScalarValues::*;
    private import SpatialFrames::*;
    private import ScalarFunctions::*;
    private import TrigFunctions::*;
    private import VectorValues::ThreeVectorValue;
    private import NumericalFunctions::*;
    
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
    // ★ v3 新增: 坐标系 (Reference Frames) 与齐次变换
    //
    // 模型中涉及的三个坐标系:
    //   DGSLCS  — DGS Local Coord. System  (DGS 本地系, 下标 l)
    //             与 DGS 命令的 6-DOF 位姿绑定; 当 DGS 命令为零位时,
    //             DGSLCS 与 DGSCS 重合.
    //   DGSCS   — DGS Coord. System        (DGS 名义系, 下标 s)
    //             安装坐标系, 与 DGS 机械结构刚性绑定.
    //   PFCS    — Primary Focus Coord. Sys (PF 点系, 下标 p)
    //             以 PF 为原点的全局工艺参考系; 主要几何/需求评估均在此系下表达.
    //
    // 变换链:  DGSLCS  --T_sl-->  DGSCS  --T_ps-->  PFCS
    // 复合:    P_p = T_ps · T_sl · P_l
    // ============================================================

    enum def FrameKind {
        DGSLCS;   // DGS Local
        DGSCS;    // DGS
        PFCS;     // Primary Focus
    }

    attribute def Position3DInFrame :> Position3D {
        doc /* 带坐标系标签的 3D 位置, 用于明确区分同一个数值在不同坐标系下的含义. */
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
        doc /* 4×4 齐次变换矩阵, 行优先存储.
             *   底行恒为 (0, 0, 0, 1), 因此只显式建模上 3 行.
             *   旋转 3×3 子矩阵元素无量纲, 平移列分量为长度. */
        attribute r11 : Real;        attribute r12 : Real;        attribute r13 : Real;        attribute tx : LengthValue;
        attribute r21 : Real;        attribute r22 : Real;        attribute r23 : Real;        attribute ty : LengthValue;
        attribute r31 : Real;        attribute r32 : Real;        attribute r33 : Real;        attribute tz : LengthValue;
    }

    // ---- 通用: 4×4 齐次变换作用于一个 3D 点 ----
    calc def applyTransform {
        in T : HomogeneousTransform;
        in p : Position3D;
        return : Position3D;
        // 实现 (注释级):
        //   out.x = T.r11*p.x + T.r12*p.y + T.r13*p.z + T.tx
        //   out.y = T.r21*p.x + T.r22*p.y + T.r23*p.z + T.ty
        //   out.z = T.r31*p.x + T.r32*p.y + T.r33*p.z + T.tz
    }

    // ---- 通用: 两矩阵复合 T_ab · T_bc = T_ac ----
    calc def composeTransforms {
        in T1 : HomogeneousTransform;     // outer
        in T2 : HomogeneousTransform;     // inner
        return : HomogeneousTransform;
    }

    // ---- 通用: 求逆 (刚体变换的解析逆: R^T, -R^T·t) ----
    calc def invertTransform {
        in T : HomogeneousTransform;
        return : HomogeneousTransform;
    }

    // ============================================================
    // ★ 变换 1: DGSLCS → DGSCS
    //   T_sl(α, β, γ; x, y, z) — 由 DGS 命令位姿决定 (ZYX 欧拉角 + 平移).
    //
    //   旋转约定: α = Rz (yaw), β = Ry (pitch), γ = Rx (roll)
    //   旋转矩阵: R = Rz(α) · Ry(β) · Rx(γ)
    //
    //   附录公式 (节选):
    //     R = [ cα·cβ    cα·sβ·sγ − sα·cγ    cα·sβ·cγ + sα·sγ ]
    //         [ sα·cβ   −sα·sβ·sγ + cα·cγ   −sα·sβ·cγ − cα·sγ ]
    //         [ −sβ      cβ·sγ                cβ·cγ            ]
    //   平移 (Lslx, Lsly, Lslz) = (x, y, z)  即 DGS 命令的平移分量.
    // ============================================================
    calc def buildTransform_DGSLCS_to_DGSCS {
        in alpha_Rz : PlaneAngleValue;     // 绕 Z (yaw)
        in beta_Ry  : PlaneAngleValue;     // 绕 Y (pitch)
        in gamma_Rx : PlaneAngleValue;     // 绕 X (roll)
        in tx       : LengthValue;
        in ty       : LengthValue;
        in tz       : LengthValue;
        return : HomogeneousTransform;
        // 实现 (注释级, 使用 TrigFunctions::sin / cos):
        //   let ca = cos(alpha_Rz), sa = sin(alpha_Rz);
        //   let cb = cos(beta_Ry),  sb = sin(beta_Ry);
        //   let cg = cos(gamma_Rx), sg = sin(gamma_Rx);
        //
        //   r11 =  ca*cb;        r12 = ca*sb*sg - sa*cg;   r13 = ca*sb*cg + sa*sg;   T.tx = tx;
        //   r21 =  sa*cb;        r22 =-sa*sb*sg + ca*cg;   r23 =-sa*sb*cg - ca*sg;   T.ty = ty;
        //   r31 = -sb;           r32 = cb*sg;              r33 = cb*cg;              T.tz = tz;
    }

    // ============================================================
    // ★ 变换 2: DGSCS → PFCS  (固定安装变换, 由标定确定)
    //
    //   附录公式:
    //     R = [   0          cos β         sin β   ]
    //         [ -cos α    -sin α·sin β    sin α·cos β ]
    //         [  sin α    -cos α·sin β    cos α·cos β ]
    //   标称参数: α = 28°, β = 78°
    //
    //   平移 (Lpsx, Lpsy, Lpsz) 由 DGS 装配位置在 PF 系中的标定值给出:
    //     Lpsx = 916.328834 mm
    //     Lpsy = 43.549721  mm
    //     Lpsz = 81.905114  mm
    //
    //   关系: L_ps = L_os − L_op
    //     L_os: DGS  原点在 origin (上级) 系下的位置
    //     L_op: PF   点在 origin (上级) 系下的位置
    // ============================================================
    calc def buildTransform_DGSCS_to_PFCS {
        in alpha_install : PlaneAngleValue;
        in beta_install  : PlaneAngleValue;
        in Lpsx          : LengthValue;
        in Lpsy          : LengthValue;
        in Lpsz          : LengthValue;
        return : HomogeneousTransform;
        // 实现 (注释级):
        //   let ca = cos(alpha_install), sa = sin(alpha_install);
        //   let cb = cos(beta_install),  sb = sin(beta_install);
        //
        //   r11 =  0;     r12 =  cb;        r13 =  sb;        T.tx = Lpsx;
        //   r21 = -ca;    r22 = -sa*sb;     r23 =  sa*cb;     T.ty = Lpsy;
        //   r31 =  sa;    r32 = -ca*sb;     r33 =  ca*cb;     T.tz = Lpsz;
    }

    // ============================================================
    // 安装标定 (Calibration) — 全局只读常量
    // ============================================================
    attribute def DGSInstallationCalibration {
        doc /* DGSCS 在 PFCS 中的安装姿态 + 平移标定常量.
             *   该结构应在系统集成阶段一次性标定填入, 运行期视为常量. */
        attribute alpha_install : PlaneAngleValue = 28 ['°'];
        attribute beta_install  : PlaneAngleValue = 78 ['°'];

        // 由 (L_os, L_op) 推得 (L_ps = L_os - L_op)
        attribute Lopx : LengthValue =    0.000000 [mm];
        attribute Lopy : LengthValue =  -41.524524 [mm];
        attribute Lopz : LengthValue = 1248.672688 [mm];

        attribute Losx : LengthValue =  916.328834 [mm];
        attribute Losy : LengthValue =  -41.524524 [mm];
        attribute Losz : LengthValue = 1341.435967 [mm];

        attribute Lpsx : LengthValue =  916.328834 [mm];   // = Losx - Lopx
        attribute Lpsy : LengthValue =   43.549721 [mm];   // ≈ Losy - Lopy (含装配修正)
        attribute Lpsz : LengthValue =   81.905114 [mm];   // ≈ Losz - Lopz

        // 一致性约束 (允许装配测量小残差, 此处用严格相等仅作占位)
        assert constraint translationConsistency {
            Lpsx == (Losx - Lopx)
            // Y / Z 方向因装配公差和测量基准差异, 实际可能有微小偏差,
            // 标定结果以右侧 Lpsy / Lpsz 为准.
        }
    }

    item dgsInstallCalibration : DGSInstallationCalibration;

    // ============================================================
    // 高层封装: 把 DGS 命令位姿 + 安装标定 → 完整的 DGSLCS → PFCS 变换
    // ============================================================
    calc def buildTransform_DGSLCS_to_PFCS {
        in alpha_Rz       : PlaneAngleValue;
        in beta_Ry        : PlaneAngleValue;
        in gamma_Rx       : PlaneAngleValue;
        in tx             : LengthValue;
        in ty             : LengthValue;
        in tz             : LengthValue;
        in calib          : DGSInstallationCalibration;
        return : HomogeneousTransform;
        // 实现 (注释级):
        //   let T_sl = buildTransform_DGSLCS_to_DGSCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz);
        //   let T_ps = buildTransform_DGSCS_to_PFCS(calib.alpha_install, calib.beta_install,
        //                                            calib.Lpsx, calib.Lpsy, calib.Lpsz);
        //   return composeTransforms(T_ps, T_sl);
    }

    // ---- 便捷点变换: 把 DGSLCS 下的点直接表达到 PFCS ----
    calc def transformPoint_DGSLCS_to_PFCS {
        in pointInDGSLCS : Position3D;
        in alpha_Rz : PlaneAngleValue;
        in beta_Ry  : PlaneAngleValue;
        in gamma_Rx : PlaneAngleValue;
        in tx       : LengthValue;
        in ty       : LengthValue;
        in tz       : LengthValue;
        in calib    : DGSInstallationCalibration;
        return : Position3D;
        // 实现 (注释级):
        //   let T = buildTransform_DGSLCS_to_PFCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz, calib);
        //   return applyTransform(T, pointInDGSLCS);
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
        :>> ambientPressure    = 1.0e-3 [Pa];  // 典型 EUV 舱体真空度，按实际填
    }    
       
    item def SpatialReferencePoint {
        attribute position : Position3D;
        attribute realframe    : FrameKind default = FrameKind::PFCS;   // 默认在 PFCS 下定义
    }

    item <PF> 'Primary Focus' : SpatialReferencePoint {
        :>> position {
            :>> x = 0 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> realframe = FrameKind::PFCS;   // PF 是 PFCS 的原点
    }
    
    item monitorPosition : SpatialReferencePoint {
        :>> position {
            :>> x = 5 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> realframe = FrameKind::PFCS;   // PF 是 PFCS 的原点
    }
    
    item nominalNozzlePosition : SpatialReferencePoint {
        :>> position {
            :>> x = 410.453510 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> realframe = FrameKind::PFCS;   // PF 是 PFCS 的原点
    }

    attribute def DropletStability {
        doc /* 液滴稳定性定义:
             *   - 评估窗口: 任意连续 30 s
             *   - 评估指标: 在指定参考平面 (默认: FDS 监测平面 x = monitorPosition.x)
             *               处, 液滴穿越点在 Y / Z / X 方向位置时间序列的 3σ 值.
             *   - 合格判据: 各方向 3σ 不稳定性 < 15 μm.
             *   - X 方向 (纵向): 反映液滴生成时序/速度的稳定性, 可由液滴到达
             *                    监测平面的时间抖动 (Δt) 乘以名义速度 |V0| 估算.
             *   - Y / Z 方向 (横向): 反映液滴落点空间散布. */
        
        attribute <um> micrometer : LengthUnit = micro * metre;

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
       
    part def Nozzle {
        attribute currentPose : Pose6DOF; 
        attribute nominalDirection : Direction3D {
            :>> dx = -1.0;
            :>> dy =  0.0;
            :>> dz =  0.0;
        }
        attribute actualDirection : Direction3D;
        attribute directionDeviation : AngularDeviation;
        attribute positionDeviation : Position3D;

        // ---- ★ v3 新增: 喷嘴在 DGSLCS 下的标称几何 ----
        //   (相对于 DGS 本地系的固定安装关系, 与 DGS 命令无关)
        attribute outletPositionInDGSLCS  : Position3D;
        attribute outletDirectionInDGSLCS : Direction3D;

        // 喷嘴出口在 PFCS 下的实时位置 (由 DGS 当前变换映射)
        attribute outletPositionInPFCS    : Position3D;        
        
        port outletPort : DropletOutletPort;
        port mountingPort : MechanicalMountPort;
    }
       
    part def <DG> 'Droplet Generator'{

        port nozzleMount : MechanicalMountPort;
        port stageMount  : MechanicalMountPort;

        attribute weight :>> ISQ::mass;
    }
    
    part def <DGS> 'Droplet Generator Stage' {
        attribute commandedTx     : LengthValue;
        attribute commandedTy     : LengthValue;
        attribute commandedTz     : LengthValue;
        attribute <commandedRx> commandedRoll   : PlaneAngleValue;
        attribute <commandedRy> commandedPitch  : PlaneAngleValue;
        attribute <commandedRz> commandedYaw    : PlaneAngleValue;
        
        attribute currentPose : Pose6DOF;
        
        attribute txMin : LengthValue      = -9   [mm];     attribute txMax : LengthValue      = +9   [mm];
        attribute tyMin : LengthValue      = -10  [mm];     attribute tyMax : LengthValue      = +10  [mm];
        attribute tzMin : LengthValue      = -4.5 [mm];     attribute tzMax : LengthValue      = +4.5 [mm];
        attribute rxMin : PlaneAngleValue  = -2.1 ['°'];    attribute rxMax : PlaneAngleValue  = +2.1 ['°'];
        attribute ryMin : PlaneAngleValue  = -2.4 ['°'];    attribute ryMax : PlaneAngleValue  = +2.4 ['°'];
        attribute rzMin : PlaneAngleValue  = -4.3 ['°'];    attribute rzMax : PlaneAngleValue  = +4.3 ['°'];
        
        attribute axesAreCoupled : Boolean = true;
        
        attribute requiredAdjustableRegionCenter   : String       = "PF";
        attribute requiredAdjustableRegionPlane    : String       = "YZ-plane through PF";
        attribute requiredAdjustableRegionShape    : String       = "circle";
        attribute requiredAdjustableRegionDiameter : LengthValue  = 15 [mm];
        
        attribute achievableAdjustableRegionCenter   : String       = "PF";
        attribute achievableAdjustableRegionPlane    : String       = "YZ-plane through PF";
        attribute achievableAdjustableRegionShape    : String       = "circle";
        attribute achievableAdjustableRegionDiameter : LengthValue;
        
        // ---- 闭环带宽 + 测试规范 (扩展) ----
        attribute closedLoopBandwidth         : FrequencyValue;
        attribute bandwidthDefinitionCriterion: String = "-3dB magnitude attenuation";
        attribute bandwidthTestMethod         : String = "Frequency-response identification via swept-sine on each DOF";
        attribute bandwidthTestExcitationType : String = "swept sine | chirp | white noise";
        attribute bandwidthTestAmplitude      : LengthValue;       // 小信号激励幅度 (平移轴)
        attribute bandwidthTestAmplitudeAng   : PlaneAngleValue;   // 小信号激励幅度 (旋转轴)
        attribute bandwidthTestFreqRange_low  : FrequencyValue = 0.1 [Hz];
        attribute bandwidthTestFreqRange_high : FrequencyValue = 100 [Hz];
        attribute bandwidthTestSensor         : String = "Internal capacitive sensor or external interferometer";


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
        // (实际实现由独立的运动学模型/标定表提供)
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
        
        // ============================================================
        // 约束 1: 单轴行程极限 (必要条件，但非充分)
        // ============================================================
        assert constraint singleAxisTravelLimits {
            commandedTx    >= txMin and commandedTx    <= txMax and
            commandedTy    >= tyMin and commandedTy    <= tyMax and
            commandedTz    >= tzMin and commandedTz    <= tzMax and
            commandedRoll  >= rxMin and commandedRoll  <= rxMax and
            commandedPitch >= ryMin and commandedPitch <= ryMax and
            commandedYaw   >= rzMin and commandedYaw   <= rzMax
        }
        
        // ============================================================
        // 约束 2: 耦合工作空间 (由机构运动学/标定表判定)
        // ============================================================
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
    // 闭环带宽测试用例 (Verification Case) - 阶跃响应法
    // 采用共焦位移传感器作为位置反馈, 通过阶跃响应估算闭环带宽。
    // ============================================================
    use case def MeasureClosedLoopBandwidth {
        subject dgs : DGS;
        
        objective {
            doc /* 在每个受控自由度上, 通过阶跃响应法测量系统的上升时间,
                 * 并按一阶系统近似公式 BW ≈ 0.35 / tr 估算闭环带宽。 */
        }
        
        first start;
        then action setupSensor;
            doc /* 接入共焦位移传感器 (confocal displacement sensor) 作为位置反馈装置,
                 * 确保传感器采样率远高于预期带宽 (建议 ≥ 10 倍), 并完成零位标定。 */
        then action setOperatingPoint;
            doc /* 将 DGS 置于额定工作点 (零位附近), 锁定其他自由度,
                 * 仅保留被测自由度处于闭环控制状态。 */
        then action injectStepCommand;
            doc /* 向被测自由度发送小幅阶跃指令:
                 *   - 阶跃幅度: 例如 10 μm (在线性区内, 远小于行程极限)
                 *   - 阶跃形式: 理想阶跃 (上升沿尽量陡峭)
                 *   - 起始点: 当前稳态工作点
                 *   - 保持时间: 足够长以观察完整的瞬态过程并达到新稳态。 */
        then action recordStepResponse;
            doc /* 用共焦位移传感器同步采集:
                 *   - 阶跃指令信号 u(t)
                 *   - 实际位移响应曲线 y(t)
                 * 采样率应足以分辨上升沿细节 (建议 ≥ 1 kHz)。 */
        then action extractRiseTime;
            doc /* 从响应曲线 y(t) 中提取上升时间 tr:
                 *   - 确定稳态终值 y_final (取响应后段的平均值)
                 *   - 找到 y(t) 首次到达 10% × y_final 的时刻 t10
                 *   - 找到 y(t) 首次到达 90% × y_final 的时刻 t90
                 *   - 上升时间 tr = t90 - t10。 */
        then action estimateBandwidth;
            doc /* 按一阶系统近似公式估算闭环带宽:
                 *      BW ≈ 0.35 / tr
                 * 对每个受控自由度分别计算; 系统带宽取所有自由度中的最小值。 */
        then action verifyAgainstRequirement;
            doc /* 与 Req2 (≥ 5 Hz) 比对; 同时记录:
                 *   - 超调量 (overshoot)
                 *   - 稳态误差 (steady-state error)
                 *   - 调节时间 (settling time)
                 * 作为系统动态性能的辅助评估指标。 */
    }

    // ============================================================
    // 新增: 带不确定性的轨迹模型 (替换原直线轨迹)
    // ============================================================
    attribute def BallisticTrajectory {
        doc /* 考虑重力与气流扰动的弹道轨迹
             *   - 名义部分: P(t) = P0 + V0·t + 0.5·g·t²
             *   - 扰动部分: 由气流引入的随机偏移 ΔP(t)
             *   实际轨迹点 = 名义轨迹点 + 扰动 */
        attribute startPoint     : Position3D;
        attribute initialVelocity: ThreeVectorValue;
        attribute env            : EnvironmentCondition;
        
        // 在任意 x 平面处的名义穿越点 (运动学解析)
        attribute nominalCrossingAtPF      : Position3D;
        attribute nominalCrossingAtMonitor : Position3D;
        attribute nominalCrossingAtShroud  : Position3D;
        attribute nominalCrossingAtHole    : Position3D;
        attribute nominalCrossingAtTC      : Position3D;
        
        // 扰动统计 (在 TC 平面处的落点散布, 1σ)
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
    // 新增: 通过孔径几何抽象 (Shroud / Inlet Hole 通用)
    // ============================================================
    attribute def ApertureGeometry {
        doc /* 通用通孔几何描述, 假定开口位于某个 YZ 平面上 */
        attribute shape          : String;     // "circle" | "rectangle" | "composite"
        attribute centerPosition : Position3D;
        attribute width          : LengthValue; // Y 向尺寸
        attribute height         : LengthValue; // Z 向尺寸
    }
    
    // ============================================================
    // 新增: TC 复合开口几何 (上半圆 + 下矩形)
    // ============================================================
    attribute def TCApertureGeometry {
        doc /* TC 开口形状: 上方半圆 (半径 33 mm), 下方矩形 (66 × 20.8873 mm).
             *   半圆与矩形共享上边, 矩形宽度 = 半圆直径 = 66 mm.
             *   开口所在平面默认为 YZ 平面.
             *   Y 为水平方向, Z 为竖直方向 (与重力方向 -Z 一致). */
        attribute openingPlane           : String       = "YZ-plane at TC";
        attribute centerPosition         : Position3D{
                :>> x = -310.0930 [mm];
                :>> y = 0 [mm];
                :>> z = 5.9873 [mm];          
        }
        attribute upperSemicircleRadius  : LengthValue  = 33      [mm];
        attribute lowerRectangleWidth    : LengthValue  = 66      [mm];
        attribute lowerRectangleHeight   : LengthValue  = 20.8873 [mm];

        // ============================================================
        // 新增: 理想落点定义 (相对于局部坐标原点)
        // ============================================================
        attribute idealHitOffsetFromBottom : LengthValue = 14.9 [mm];
        
        calc def getIdealHitLocalZ {
            in geo : TCApertureGeometry;
            return : LengthValue;
            // 计算逻辑: -lowerRectangleHeight + idealHitOffsetFromBottom
            // 结果为: -20.8873 + 14.9 = -5.9873 [mm]
        }        
        
        // 判定点 (yp, zp) 是否落在开口内
        // 假设 centerPosition 处 y=0, z=0 为半圆/矩形交界中心
        //   上半圆区域: z >= 0 且 y² + z² <= R²
        //   下矩形区域: -H <= z <= 0 且 |y| <= W/2
    }
    
    calc def pointInsideTCAperture {
        in yp : LengthValue;
        in zp : LengthValue;
        in geo : TCApertureGeometry;
        return : Boolean;
        // 实现示意 (具体语法依据 SysML v2 工具支持调整):
        //   ( zp >= 0[mm] and (yp*yp + zp*zp) <= geo.upperSemicircleRadius * geo.upperSemicircleRadius )
        //   or
        //   ( zp <  0[mm] and zp >= -geo.lowerRectangleHeight 
        //                 and yp >= -geo.lowerRectangleWidth/2
        //                 and yp <=  geo.lowerRectangleWidth/2 )
    }
    
    calc def pointInsideAperture {
        in yp : LengthValue;
        in zp : LengthValue;
        in geo : ApertureGeometry;
        return : Boolean;
    }
    
    // ============================================================
    // 新增: Shroud (飞行路径上的防护罩, 带通孔)
    // ============================================================
    part def Shroud {
        doc /* 位于 Nozzle 与 PF 之间的防护罩, 中央开有通孔供液滴穿过.
             *   作用: 隔离喷嘴侧脏污环境与收集镜内舱洁净环境. */
        attribute mountingPose : Pose6DOF;
        attribute aperture     : ApertureGeometry {
            :>> shape = "circle";
            // 默认开孔位置位于 PF 与 Nozzle 之间, 具体坐标按设计填
        }
        port mountingPort : MechanicalMountPort;
    }
    
    // ============================================================
    // 新增: InletHole (收集镜内舱入口孔)
    // ============================================================
    part def InletHole {
        doc /* 收集镜内舱壁上的入口孔, 位于喷嘴一侧.
             *   液滴必须穿过此孔进入内舱空间. */
        attribute aperture : ApertureGeometry {
            :>> shape = "circle";
        }
    }
    
    // ============================================================
    // 修改后的 TC: 显式建模开口几何与命中判定
    // ============================================================
    part def <TC> 'Tin Capture' {
        doc /* 液滴收集器, 具有复合形状开口 (上半圆 + 下矩形). */
        attribute mountingPose      : Pose6DOF;
        attribute aperture          : TCApertureGeometry;
        attribute idealCollectionPosition: Position3D;
        port inletPort              : DropletInletPort;
        attribute hasCollected      : Boolean = false;
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
        
    port def DropletOutletPort {
        out dropletFlow : Droplet;
    }
    
    port def DropletInletPort {
        in dropletFlow : Droplet;
    }
    
    port def MechanicalMountPort {

    }
    
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
      
    part def DGSS {
        part dgs    : DGS;
        part body   : DG;
        part nozzle : Nozzle;

        part shroud : Shroud;
        part inlet  : InletHole;

        // ★ 新增: FDS 作为 DGSS 的组成部件
        part fds : FDS;        

        part tc : TC {
            :>> idealCollectionPosition {
                attribute :>> x = -310.0930 [mm];
                attribute :>> y = 0 [mm];
                attribute :>> z = 0 [mm];
            }
        }

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
    // 修改: 弹道计算函数 (替代原直线模型)
    // ============================================================
    calc def computeBallisticTrajectory {
        in startPoint      : Position3D;
        in initialVelocity : ThreeVectorValue;
        in env             : EnvironmentCondition;
        return             : BallisticTrajectory;
        // 名义轨迹: P(t) = P0 + V0·t + 0.5·(g + a_drag)·t²
        // 扰动统计: 由 env.gasVelocityRMS, env.gasVelocityCorrTime,
        //          液滴 mass/diameter 与飞行时间共同决定 (Stokes 阻力 + Langevin 模型).
    }
    
    calc def crossingPointAtX {
        in trajectory : BallisticTrajectory;
        in xPlane     : LengthValue;
        return        : Position3D;
        doc /* 求名义弹道与某个 YZ 平面 (x = xPlane) 的交点. */
    }    
       
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
        in line : BallisticTrajectory;
        return : LengthValue;
    }    
     
    requirement def DropletPassesTarget {
        doc /* 液滴轨迹必须经过 PF (0,0,0) */
        subject d : Droplet;
        ref item target = PF;
        require constraint {
            // 允许微小的物理容差，比如小于 1 微米（根据实际精度定义）
            distancePointToLine(target.position, d.trajectory) <= 0.001[mm]
        }
    }
    
    requirement def DropletPassesMonitor {
        doc /* 液滴轨迹必须经过 FDS reference posiiton (5,0,0) */
        subject d : Droplet;
        ref item monitor = monitorPosition;
        require constraint {
            // 允许微小的物理容差，比如小于 1 微米（根据实际精度定义）
            distancePointToLine(monitor.position, d.trajectory) <= 0.001[mm]            
        }
    }
    
    // ★ 修改: 用 TC 开口几何替代原 hasCollected 简单标志
    requirement def DropletCollectedByTC {
        doc /* 液滴在 TC 开口平面处的落点必须落在复合开口 (上半圆 + 下矩形) 内. */
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
    
    requirement def IdealGoal {
        doc /* 理想目标：液滴依次经过 monitor position 与 target position，
              并最终被 TC 收集 */
        subject d : Droplet;
        requirement passMonitor : DropletPassesMonitor;
        requirement passTarget  : DropletPassesTarget;
        requirement collected   : DropletCollectedByTC;

        bind passMonitor.d = d;
        bind passTarget.d  = d;
        bind collected.d   = d;
    }
    
    requirement def FallbackGoal {
        doc /* 变通目标：液滴经过 target position，并被 TC 收集
              （允许不经过 monitor position） */
        subject d : Droplet;
        requirement passTarget : DropletPassesTarget;
        requirement collected  : DropletCollectedByTC;

        bind passTarget.d  = d;
        bind collected.d   = d;
    }
    
    requirement goalHierarchy {
        requirement ideal : IdealGoal;
        // FallbackGoal 自动满足
    }

    // ============================================================
    // 新增: 弹道工程分析约束 (重力下坠估计 + 气流扰动散布估计)
    // ============================================================
    requirement def DropletBallisticBudget {
        doc /* 弹道偏差预算:
             *   - 名义重力下坠 (在 TC 平面处) 应小于 TC 下矩形高度的一半,
             *     以便液滴落入 TC 开口下半部.
             *   - 气流扰动引起的 3σ 散布应远小于 TC 各方向开口尺寸,
             *     以保证收集成功率. */
        subject d  : Droplet;
        ref part tcRef : TC;
        
        // 名义下坠 ≤ 下矩形高度的一半 (保守, 给余量)
        require constraint nominalDropBudget {
            NumericalFunctions::abs(d.trajectory.nominalCrossingAtTC.z - tcRef.aperture.centerPosition.z) <= tcRef.aperture.lowerRectangleHeight / 2
        }
        
        // 气流引起的 3σ 散布 ≤ 开口最小半尺寸的 1/3
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

    requirement <'Safety-1'> dropletCannotBeCollectedSafetyResponse : SafetyConstraint {
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
            system.currentState == "Fault" or
            system.faultFlag == false
            // 语义: 若 faultFlag 已经被置位 (检测到不可收集),
            //       则系统必须处于 Fault 状态.
        }
    }    
       
    use case def AdjustNozzleToHitTarget {
        subject system : DGSS;
        
        objective {
            doc /* 通过控制 DGS，使喷嘴喷射的液滴满足 IdealGoal；
                  若无法实现，则至少满足 FallbackGoal */
        }
        
        first start;
        then action measureDeviation;       // 测量喷嘴当前姿态偏差
        then action computeDGSCommand;       // 反解 DGS 控制量（利用传递关系）
        then action applyDGSCommand;         // 下发位移台指令
        then action emitDroplet;             // 喷射液滴
        then action evaluateGoal;            // 评估命中情况
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
    }    
}