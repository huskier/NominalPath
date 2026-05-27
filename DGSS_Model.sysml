package DGSS_Model {

    private import SI::*;
    private import ISQ::*;
    private import ScalarValues::*;
    private import SpatialFrames::*;
    private import ScalarFunctions::*;
    private import TrigFunctions::*;
    private import VectorValues::ThreeVectorValue;
    
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
    // 新增: 物理环境参数 (重力 + 气流)
    // ============================================================
    attribute def EnvironmentCondition {
        doc /* 飞行环境条件，影响液滴轨迹 */
        attribute gravity : AccelerationValue = 9.80665 [m/s²];
        attribute gravityDirection : Direction3D {
            :>> dx =  0.0;
            :>> dy =  0.0;
            :>> dz = -1.0;
        }
        attribute ambientPressure        : PressureValue;     // 舱内残余气压
        attribute meanGasVelocity        : ThreeVectorValue;  // 平均气流速度
        attribute gasVelocityRMS         : SpeedValue;        // 气流脉动 RMS
        attribute gasVelocityCorrTime    : TimeValue;         // 气流相关时间
    }
    
    item ambientEnvironment : EnvironmentCondition {
        :>> ambientPressure    = 1.0e-3 [Pa];  // 典型 EUV 舱体真空度，按实际填
    }    
       
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
        
        
        // ============================================================
        // 端口
        // ============================================================
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
        
        requirement <Req1> dropletAdjustableRangeRequirement {
            doc /* 在过PF点的YZ平面上，以PF点为中心，
                 * 液滴落点应在直径不小于 15 mm 的圆内任意可调。 */
            
            assume constraint regionGeometryMatch {
                achievableAdjustableRegionCenter == requiredAdjustableRegionCenter and
                achievableAdjustableRegionPlane  == requiredAdjustableRegionPlane  and
                achievableAdjustableRegionShape  == requiredAdjustableRegionShape
            }
            
            require constraint diameterCoverage {
                achievableAdjustableRegionDiameter >= requiredAdjustableRegionDiameter
            }
        }
        
        requirement <Req2> closedLoopBandwidthRequirement {
            doc /* DGS 各受控自由度闭环控制带宽 ≥ 5 Hz, 以 -3dB 幅值衰减为判定准则.
                 *   测试方法见 bandwidthTestMethod / Verification case BandwidthTest. */
            require constraint bandwidthLowerBound {
                closedLoopBandwidth >= 5 [Hz]
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
        attribute collectionPosition: Position3D; // 开口中心位置 (与 aperture.centerPosition 对齐)
        port inletPort              : DropletInletPort;
        attribute hasCollected      : Boolean = false;
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

        part tc : TC {
            :>> collectionPosition {
                attribute :>> x = -310.0930 [mm];
                attribute :>> y = 0 [mm];
                attribute :>> z = 0 [mm];
            }
        }
        
        ref item targetRef  = PF;
        ref item monitorRef = monitorPosition;
        ref item envRef     = ambientEnvironment;
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
            abs(d.trajectory.nominalCrossingAtTC.z
                - tcRef.aperture.centerPosition.z)
                <= tcRef.aperture.lowerRectangleHeight / 2
        }
        
        // 气流引起的 3σ 散布 ≤ 开口最小半尺寸的 1/3
        require constraint dispersionBudget {
            3 * d.trajectory.crossingDispersionAtTC_y
                <= tcRef.aperture.lowerRectangleWidth / 6  and
            3 * d.trajectory.crossingDispersionAtTC_z
                <= tcRef.aperture.lowerRectangleHeight / 6
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
}

