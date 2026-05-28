package DGSS_Model {

    doc /* Version: 2.0
         * Author: Simon
         * Date: 2026-05-28
         * Change history:
         *   v1.0 (2026-05-28): 初始版本
         *   v2.0 (2026-05-28): 增加 FDS 组件、DGSS 状态机、液滴稳定性定义、
         *                      SafetyConstraint, DGS 调节分辨率/承重/承压/重量需求.
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
    }

    item <PF> 'Primary Focus' : SpatialReferencePoint {
        :>> position {
            :>> x = 0 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
    }

    item monitorPosition : SpatialReferencePoint {
        :>> position {
            :>> x = 5 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
    }

    item nominalNozzlePosition : SpatialReferencePoint {
        :>> position {
            :>> x = 410.453510 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
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
}
