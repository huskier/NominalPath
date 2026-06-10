package DGSS_Model {

    doc /* Version: 5.2
         * Author: Simon  (变换层优化 by 协作迭代)
         * Date: 2026-06-10
         * Change history (见 ChangeLog.txt). v5.2 修订摘要:
         *   FIX④ [H] 运行期可变量由绑定 `=` 改为初值 `:=` (uncollectableDetected /
         *            lastDropletMissedTC / resetCommand / hasCollected),
         *            消除 "绑定恒假 → faultFlag 恒假 → Fault 不可达" 的回归.
         *   FIX⑤ [H] 仅含 if 守卫的迁移改为 accept when 变化触发
         *            (源状态含无界 do 动作时, 守卫式完成迁移永无评估时机).
         *   FIX⑥ [M] DGSS 状态机改为 exhibit state, 确保由 DGSS 实例履行.
         *   FIX⑦ [M] DG.weight 由非法重定义 `:>> ISQ::mass` 改为 `: MassValue`.
         *   FIX⑧ [M] do 动作内的故障锁存由单次 assign 改为 loop + accept when
         *            变化事件驱动, 实现真正的"持续锁存"语义.
         *   FIX⑨ [M] Req-1..6 / Req-FDS-1 提升至包级, 约束统一经 subject 引用,
         *            并在 DGSS 中以 satisfy 建立需求-主体绑定与追溯.
         *   FIX⑩ [M] translationConsistency 由精确 `==` 改为 0.1 µm 容差形式
         *            (标定字面量第 6 位小数存在 ~0.6 nm 舍入, 精确相等会判 false).
         *   其余: pointInsideTCAperture / getIdealHitLocalZ 落地为可求值表达式并
         *         提升至包级; 液滴/指令交互由裸 connect 升级为 flow; use case 中
         *         doc 移入对应 action 体内 (原挂错宿主); 删除冗余 `== true`;
         *         micrometer → micrometre 拼写统一; 修正两处复制错误注释.
         */

    private import SI::*;
    private import ISQ::*;
    private import ScalarValues::*;
    private import MeasurementReferences::*;
    private import ScalarFunctions::*;
    private import TrigFunctions::*;
    private import VectorValues::ThreeVectorValue;
    private import NumericalFunctions::*;

    // [M] 若导入的 SI 库已声明 micrometre/<um>, 本地声明将按名称遮蔽之 (合法);
    //     请在工具中确认无歧义告警.
    attribute <um> micrometre : LengthUnit = micro * metre;
  
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
    // 四坐标系与齐次变换链 (v5.0 链条)
    //
    //   树 (根在 MCS):
    //     MCS (root)
    //      └─ PFCS   : transformation = PFCS_to_MCS  (R=Rx(α), T=L_op),  source=MCS
    //          └─ DGSCS : transformation = DGSCS_to_PFCS (复合, T=L_ps),  source=PFCS
    //              └─ DGSLCS : transformation = DGSLCS_to_DGSCS (命令6DOF), source=DGSCS
    //
    //   数据流 (子→父): DGSLCS → DGSCS → PFCS → MCS
    //   规则: [p]_父 = M · [p]_子；M = 子系.transformation；source=父，target=that=子。
    //
    //   固定安装关系 (已验证):
    //     α = 28°  Rx(α) 连接 PFCS ↔ MCS   (PFCS_to_MCS 平移 = L_op)
    //     β = 78°  连接 DGSCS ↔ (经 MCS) → 并入 DGSCS_to_PFCS 复合
    //     L_ps = (916.328834, 43.549721, 81.905114) mm = DGSCS→PFCS 平移列.
    //     L_op = (0, -41.524524, 1248.672688) mm       = PFCS→MCS  平移列.
    // ============================================================

    enum def FrameKind {
        DGSLCS_FRAME;         // DGS Local (DGS局部坐标系，l 系)
        DGSCS_FRAME;          // DGS (机器系中的 s 系)
        MCS_FRAME;            // 机器坐标系 (Machine Coordinate System, o 系)
        PFCS_FRAME;           // Primary Focus
    }

    enum def DGSSOperatingState {
        initialization;
        idle;
        closedLoopLocking;
        fault;
    }   

    attribute def Position3DInFrame :> Position3D {
        doc /* 带坐标系标签的 3D 位置, 明确区分同一数值在不同坐标系下的含义. */
        attribute coordinateFrame : FrameKind;
    }

    // ============================================================
    // 位移向量基础定义 (附录 "位置/位移向量定义")
    //
    //   设点 P, O 在某坐标系中的位置向量分别为 p, o, 则从 O 指向 P 的
    //   位移向量定义为:
    //                       L_op = p - o
    //   其分量 L_opx, L_opy, L_opz 即 L_op 在 x/y/z 轴上的标量分量.
    //
    //   本模型中所有 L_** 量 (L_op, L_os, L_ps, …) 均遵循此定义:
    //   下标 (起点 终点) → 向量 = 终点位置 − 起点位置.
    //   注意各 L_** 量须配合其所在坐标系解读 (见 coordinateFrame 标签).
    // ============================================================
    attribute def DisplacementVector {
        doc /* 从起点 O 指向终点 P 的位移向量 L_op = p - o.
             *   x/y/z 为其三个标量分量; coordinateFrame 标明该向量在哪个坐标系下表达. */
        attribute x               : LengthValue;
        attribute y               : LengthValue;
        attribute z               : LengthValue;
        attribute coordinateFrame : FrameKind;
    }

    // ---- 构造: 由起点 O 与终点 P 的位置求位移向量 L = P - O ----
    calc def displacementVector {
        in originPoint : Position3D;     // O
        in targetPoint : Position3D;     // P
        in coordinateFrame       : FrameKind;      // 两点所在坐标系 (须一致)
        return         : DisplacementVector;
        // 实现 (注释级):
        //   out.x = targetPoint.x - originPoint.x
        //   out.y = targetPoint.y - originPoint.y
        //   out.z = targetPoint.z - originPoint.z
        //   out.coordinateFrame = coordinateFrame
    }

    // ============================================================
    // 角度正方向约定 (附录 "坐标系旋转方向规定")
    //   按右手定则: 右手大拇指指向某轴 (X / Y / Z) 正方向时,
    //   其余四指弯曲方向即为该轴旋转角度的正方向.
    //   本模型中所有 PlaneAngleValue 旋转量 (α, β, γ, roll/pitch/yaw,
    //   命令角与安装角) 均遵循此右手定则正方向约定.
    // ============================================================
    attribute def RotationSignConvention {
        doc /* 右手定则角度正方向约定. 该结构作为全局约定声明, 供旋转量引用.
             *   isRightHanded = true 表示遵循右手定则. */
        attribute isRightHanded : Boolean = true;
        attribute description   : String  =
            "Right-hand rule: thumb along +axis, curled fingers give positive angle direction.";
    }

    attribute globalRotationConvention : RotationSignConvention;    

    // =================================================================
    //  坐标系类型: 三轴单位均为 mm 的 3D 正交坐标系
    //  [M] millimetre 的确切元素名请核 SI 库 (可能是 millimetre 或 mm 别名)
    // =================================================================
    attribute def MillimetreFrame3d :> ThreeDCoordinateFrame {
        attribute :>> mRefs = (millimetre, millimetre, millimetre);
        attribute :>> isOrthogonal = true;
        attribute kind : FrameKind;
    }    

    // ---- 应用变换: [p]_source = M · [p]_target ----  [H]
    calc def applyTransform {
        in M : AffineTransformationMatrix3d;
        in p : Position3D;

        private attribute px : Real = p.x / (1 [mm]);
        private attribute py : Real = p.y / (1 [mm]);
        private attribute pz : Real = p.z / (1 [mm]);

        private attribute r11 : Real = M.rotationMatrix.elements#(1);
        private attribute r12 : Real = M.rotationMatrix.elements#(2);
        private attribute r13 : Real = M.rotationMatrix.elements#(3);
        private attribute r21 : Real = M.rotationMatrix.elements#(4);
        private attribute r22 : Real = M.rotationMatrix.elements#(5);
        private attribute r23 : Real = M.rotationMatrix.elements#(6);
        private attribute r31 : Real = M.rotationMatrix.elements#(7);
        private attribute r32 : Real = M.rotationMatrix.elements#(8);
        private attribute r33 : Real = M.rotationMatrix.elements#(9);
        private attribute tx  : Real = M.translationVector.elements#(1);
        private attribute ty  : Real = M.translationVector.elements#(2);
        private attribute tz  : Real = M.translationVector.elements#(3);

        return result : Position3D {
            attribute :>> x = (r11*px + r12*py + r13*pz + tx) * (1 [mm]);
            attribute :>> y = (r21*px + r22*py + r23*pz + ty) * (1 [mm]);
            attribute :>> z = (r31*px + r32*py + r33*pz + tz) * (1 [mm]);
        }
    }

    // ---- 复合: composeTransforms(T1, T2) = T1 · T2 (T1 外层) ----  [M] 复核索引
    calc def composeTransforms {
        in T1 : AffineTransformationMatrix3d;
        in T2 : AffineTransformationMatrix3d;

        private attribute a11 : Real = T1.rotationMatrix.elements#(1);
        private attribute a12 : Real = T1.rotationMatrix.elements#(2);
        private attribute a13 : Real = T1.rotationMatrix.elements#(3);
        private attribute a21 : Real = T1.rotationMatrix.elements#(4);
        private attribute a22 : Real = T1.rotationMatrix.elements#(5);
        private attribute a23 : Real = T1.rotationMatrix.elements#(6);
        private attribute a31 : Real = T1.rotationMatrix.elements#(7);
        private attribute a32 : Real = T1.rotationMatrix.elements#(8);
        private attribute a33 : Real = T1.rotationMatrix.elements#(9);
        private attribute t1x : Real = T1.translationVector.elements#(1);
        private attribute t1y : Real = T1.translationVector.elements#(2);
        private attribute t1z : Real = T1.translationVector.elements#(3);

        private attribute b11 : Real = T2.rotationMatrix.elements#(1);
        private attribute b12 : Real = T2.rotationMatrix.elements#(2);
        private attribute b13 : Real = T2.rotationMatrix.elements#(3);
        private attribute b21 : Real = T2.rotationMatrix.elements#(4);
        private attribute b22 : Real = T2.rotationMatrix.elements#(5);
        private attribute b23 : Real = T2.rotationMatrix.elements#(6);
        private attribute b31 : Real = T2.rotationMatrix.elements#(7);
        private attribute b32 : Real = T2.rotationMatrix.elements#(8);
        private attribute b33 : Real = T2.rotationMatrix.elements#(9);
        private attribute t2x : Real = T2.translationVector.elements#(1);
        private attribute t2y : Real = T2.translationVector.elements#(2);
        private attribute t2z : Real = T2.translationVector.elements#(3);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    a11*b11+a12*b21+a13*b31,  a11*b12+a12*b22+a13*b32,  a11*b13+a12*b23+a13*b33,
                    a21*b11+a22*b21+a23*b31,  a21*b12+a22*b22+a23*b32,  a21*b13+a22*b23+a23*b33,
                    a31*b11+a32*b21+a33*b31,  a31*b12+a32*b22+a33*b32,  a31*b13+a32*b23+a33*b33 );
            }
            attribute :>> translationVector {
                attribute :>> elements = (
                    a11*t2x + a12*t2y + a13*t2z + t1x,
                    a21*t2x + a22*t2y + a23*t2z + t1y,
                    a31*t2x + a32*t2y + a33*t2z + t1z );
            }
        }
    }

    // ---- 刚体逆: R'=Rᵀ, t'=-Rᵀ·t ----  [M] 复核
    calc def invertTransform {
        in T : AffineTransformationMatrix3d;

        private attribute r11 : Real = T.rotationMatrix.elements#(1);
        private attribute r12 : Real = T.rotationMatrix.elements#(2);
        private attribute r13 : Real = T.rotationMatrix.elements#(3);
        private attribute r21 : Real = T.rotationMatrix.elements#(4);
        private attribute r22 : Real = T.rotationMatrix.elements#(5);
        private attribute r23 : Real = T.rotationMatrix.elements#(6);
        private attribute r31 : Real = T.rotationMatrix.elements#(7);
        private attribute r32 : Real = T.rotationMatrix.elements#(8);
        private attribute r33 : Real = T.rotationMatrix.elements#(9);
        private attribute tx  : Real = T.translationVector.elements#(1);
        private attribute ty  : Real = T.translationVector.elements#(2);
        private attribute tz  : Real = T.translationVector.elements#(3);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    r11, r21, r31,
                    r12, r22, r32,
                    r13, r23, r33 );
            }
            attribute :>> translationVector {
                attribute :>> elements = (
                    -(r11*tx + r21*ty + r31*tz),
                    -(r12*tx + r22*ty + r32*tz),
                    -(r13*tx + r23*ty + r33*tz) );
            }
        }
    }

    // =================================================================
    //  变换 0 (实时): DGSLCS → DGSCS，命令 6-DOF 决定
    //   R = Rz(α)·Ry(β)·Rx(γ)，正确符号；[p]_DGSCS = M·[p]_DGSLCS
    // =================================================================  [H]
    calc def buildTransform_DGSLCS_to_DGSCS {
        in alpha_Rz : PlaneAngleValue;   // yaw  (Rz)
        in beta_Ry  : PlaneAngleValue;   // pitch(Ry)
        in gamma_Rx : PlaneAngleValue;   // roll (Rx)
        in tx : LengthValue;  in ty : LengthValue;  in tz : LengthValue;

        private attribute a : Real = alpha_Rz / (1 [radian]);
        private attribute b : Real = beta_Ry  / (1 [radian]);
        private attribute g : Real = gamma_Rx / (1 [radian]);
        private attribute ca : Real = cos(a);  private attribute sa : Real = sin(a);
        private attribute cb : Real = cos(b);  private attribute sb : Real = sin(b);
        private attribute cg : Real = cos(g);  private attribute sg : Real = sin(g);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    ca*cb,  ca*sb*sg - sa*cg,  ca*sb*cg + sa*sg,
                    sa*cb,  sa*sb*sg + ca*cg,  sa*sb*cg - ca*sg,
                     -sb,   cb*sg,             cb*cg );
            }
            attribute :>> translationVector {
                attribute :>> elements = ( tx/(1 [mm]), ty/(1 [mm]), tz/(1 [mm]) );
            }
        }
    }

    // =================================================================
    //  PFCS → MCS  (新树边, α=28°)   [H] 用户提供, 已验证 = invert(MCS_to_PFCS)
    //   R = Rx(α) = [1,0,0; 0,cα,-sα; 0,sα,cα] ,  T = L_op = (Lopx,Lopy,Lopz)
    //   [p]_MCS = M · [p]_PFCS
    // =================================================================
    calc def buildTransform_PFCS_to_MCS {
        in alpha_install : PlaneAngleValue;
        in Lopx : LengthValue;  in Lopy : LengthValue;  in Lopz : LengthValue;

        private attribute a : Real = alpha_install / (1 [radian]);
        private attribute ca : Real = cos(a);  private attribute sa : Real = sin(a);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    1,    0,    0,
                    0,   ca,  -sa,
                    0,   sa,   ca );
            }
            attribute :>> translationVector {
                attribute :>> elements = ( Lopx/(1 [mm]), Lopy/(1 [mm]), Lopz/(1 [mm]) );
            }
        }
    }

    // =================================================================
    //  MCS → PFCS  (= invert(PFCS_to_MCS), α=28°)   [M] 转写自 v4.3 验证注释
    // =================================================================
    calc def buildTransform_MCS_to_PFCS {
        in alpha_install : PlaneAngleValue;
        in Lopx : LengthValue;  in Lopy : LengthValue;  in Lopz : LengthValue;

        private attribute a : Real = alpha_install / (1 [radian]);
        private attribute ca : Real = cos(a);  private attribute sa : Real = sin(a);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    1,    0,    0,
                    0,   ca,   sa,
                    0,  -sa,   ca );
            }
            attribute :>> translationVector {
                attribute :>> elements = (
                    (0 [mm] - Lopx) / (1 [mm]),
                    (0 [mm] - Lopy*ca - Lopz*sa) / (1 [mm]),
                    (Lopy*sa - Lopz*ca) / (1 [mm]) );
            }
        }
    }

    // =================================================================
    //  DGSCS → MCS  (β=78°, 单段)   [M]
    // =================================================================
    calc def buildTransform_DGSCS_to_MCS {
        in beta_install : PlaneAngleValue;
        in Losx : LengthValue;  in Losy : LengthValue;  in Losz : LengthValue;

        private attribute b : Real = beta_install / (1 [radian]);
        private attribute cb : Real = cos(b);  private attribute sb : Real = sin(b);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    0,   cb,   sb,
                   -1,    0,    0,
                    0,  -sb,   cb );
            }
            attribute :>> translationVector {
                attribute :>> elements = ( Losx/(1 [mm]), Losy/(1 [mm]), Losz/(1 [mm]) );
            }
        }
    }

    // =================================================================
    //  DGSCS → PFCS  (复合, 树边)  [M]  平移列 = L_ps (附录已验证)
    // =================================================================
    calc def buildTransform_DGSCS_to_PFCS {
        in alpha_install : PlaneAngleValue;
        in beta_install  : PlaneAngleValue;
        in Lopx : LengthValue;  in Lopy : LengthValue;  in Lopz : LengthValue;
        in Losx : LengthValue;  in Losy : LengthValue;  in Losz : LengthValue;

        private attribute a : Real = alpha_install / (1 [radian]);
        private attribute b : Real = beta_install  / (1 [radian]);
        private attribute ca : Real = cos(a);  private attribute sa : Real = sin(a);
        private attribute cb : Real = cos(b);  private attribute sb : Real = sin(b);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                     0,    cb,      sb,
                   -ca,  -sa*sb,  sa*cb,
                    sa,  -ca*sb,  ca*cb );
            }
            attribute :>> translationVector {
                attribute :>> elements = (
                    (Losx - Lopx) / (1 [mm]),
                    ((Losy-Lopy)*ca + (Losz-Lopz)*sa) / (1 [mm]),
                    ((Lopy-Losy)*sa + (Losz-Lopz)*ca) / (1 [mm]) );
            }
        }
    }

    // =================================================================
    //  PFCS → DGSCS  (= invert(DGSCS_to_PFCS))  [M]
    // =================================================================
    calc def buildTransform_PFCS_to_DGSCS {
        in alpha_install : PlaneAngleValue;
        in beta_install  : PlaneAngleValue;
        in Lopx : LengthValue;  in Lopy : LengthValue;  in Lopz : LengthValue;
        in Losx : LengthValue;  in Losy : LengthValue;  in Losz : LengthValue;

        private attribute a : Real = alpha_install / (1 [radian]);
        private attribute b : Real = beta_install  / (1 [radian]);
        private attribute ca : Real = cos(a);  private attribute sa : Real = sin(a);
        private attribute cb : Real = cos(b);  private attribute sb : Real = sin(b);

        return result : AffineTransformationMatrix3d {
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    0,    -ca,     sa,
                   cb,  -sa*sb,  -ca*sb,
                   sb,   sa*cb,   ca*cb );
            }
            attribute :>> translationVector {
                attribute :>> elements = (
                    (Losy - Lopy) / (1 [mm]),
                    ((Lopx-Losx)*cb - (Lopz-Losz)*sb) / (1 [mm]),
                    ((Lopx-Losx)*sb + (Lopz-Losz)*cb) / (1 [mm]) );
            }
        }
    }


    // ============================================================
    // 安装标定 (Calibration) — 全局只读常量, 系统集成阶段标定填入.
    //   L_op, L_os 在 MCS(机器系) 下; L_ps 在 PFCS 下.
    //   依据位移向量定义 L_ab = (终点 b 位置) − (起点 a 位置):
    //     L_op = (PF 点位置 p) − (MCS 原点 o)   [在 MCS 下]
    //     L_os = (DGS 原点 s 位置) − (MCS 原点 o) [在 MCS 下]
    //     L_ps = (DGS 原点 s 位置) − (PF 点 p 位置) [在 PFCS 下]
    // ============================================================
    attribute def DGSInstallationCalibration {
        doc /* DGS 安装姿态与平移标定. α 连接 PFCS↔MCS, β 连接 MCS↔DGSCS.
             *   均为系统集成期一次性标定的固定值, 运行期视为常量.
             *   旋转角遵循右手定则 (见 RotationSignConvention). */
        attribute alpha_install : PlaneAngleValue = 28 ['°'];   // PFCS ↔ MCS
        attribute beta_install  : PlaneAngleValue = 78 ['°'];   // MCS ↔ DGSCS

        // ---- 标定原始量: PF 点与 DGS 原点在 MCS 下的位置 ----
        //   (位移向量 L_op / L_os 的终点; MCS 原点 o = (0,0,0) 为起点)
        attribute P_in_MCS : Position3D {        // PF 点在 MCS 下
            :>> x =    0.000000 [mm];
            :>> y =  -41.524524 [mm];
            :>> z = 1248.672688 [mm];
        }
        attribute S_in_MCS : Position3D {        // DGS 原点在 MCS 下
            :>> x =  916.328834 [mm];
            :>> y =  -41.524524 [mm];
            :>> z = 1341.435967 [mm];
        }

        // ---- 位移向量 L_op = P - o (在 MCS, o=原点) ----
        attribute L_op : DisplacementVector {
            :>> x =    0.000000 [mm];
            :>> y =  -41.524524 [mm];
            :>> z = 1248.672688 [mm];
            :>> coordinateFrame = FrameKind::MCS_FRAME;
        }
        // ---- 位移向量 L_os = S - o (在 MCS, o=原点) ----
        attribute L_os : DisplacementVector {
            :>> x =  916.328834 [mm];
            :>> y =  -41.524524 [mm];
            :>> z = 1341.435967 [mm];
            :>> coordinateFrame = FrameKind::MCS_FRAME;
        }
        // ---- 位移向量 L_ps = S - P (在 PFCS) = DGSCS→PFCS 复合平移列 T ----
        attribute L_ps : DisplacementVector {
            :>> x =  916.328834 [mm];
            :>> y =   43.549721 [mm];
            :>> z =   81.905114 [mm];
            :>> coordinateFrame = FrameKind::PFCS_FRAME;
        }

        // 为向后兼容保留的标量别名 (供既有约束/函数引用)
        attribute Lopx : LengthValue = L_op.x;   attribute Lopy : LengthValue = L_op.y;   attribute Lopz : LengthValue = L_op.z;
        attribute Losx : LengthValue = L_os.x;   attribute Losy : LengthValue = L_os.y;   attribute Losz : LengthValue = L_os.z;
        attribute Lpsx : LengthValue = L_ps.x;   attribute Lpsy : LengthValue = L_ps.y;   attribute Lpsz : LengthValue = L_ps.z;

        // ---- 位移向量定义自洽: L_op / L_os 由 MCS 位置点按 L = P - o 派生 ----
        //   (MCS 原点 o = (0,0,0), 故 L_op = P_in_MCS, L_os = S_in_MCS)
        assert constraint displacementDefinitionConsistency {
            L_op.x == P_in_MCS.x and L_op.y == P_in_MCS.y and L_op.z == P_in_MCS.z and
            L_os.x == S_in_MCS.x and L_os.y == S_in_MCS.y and L_os.z == S_in_MCS.z
        }

        // ---- 标定一致性 (cos/sin 消量纲) ----  [H]
        // ★ FIX⑩: 标定字面量按 6 位小数手填, L_ps.y 精确值 ≈ 43.549722 (差 ~0.6 nm),
        //   精确 `==` 会因舍入判 false. 改为 0.1 µm 容差 (远高于舍入、远低于物理意义).
        assert constraint translationConsistency {
            NumericalFunctions::abs( Lpsx - (Losx - Lopx) ) <= 1.0e-4 [mm] and
            NumericalFunctions::abs( Lpsy - ((Losy - Lopy) * cos(alpha_install/(1 [radian])) + (Losz - Lopz) * sin(alpha_install/(1 [radian]))) ) <= 1.0e-4 [mm] and
            NumericalFunctions::abs( Lpsz - ((Lopy - Losy) * sin(alpha_install/(1 [radian])) + (Losz - Lopz) * cos(alpha_install/(1 [radian]))) ) <= 1.0e-4 [mm]
        }
    }

    attribute dgsInstallCalibration : DGSInstallationCalibration;

    // =================================================================
    //  四坐标系实例 + 安置关系 (新链: MCS 根)
    //  [M] frame↔矩阵的等值绑定语法请在工具中复核；若 == 整体比较因
    //      source/target 参与而失败，改为只绑 rotationMatrix/translationVector。
    // =================================================================
    attribute mcs : MillimetreFrame3d {                         // 根参考系，无 transformation
        attribute :>> kind = FrameKind::MCS_FRAME;
    }

    attribute pfcs : MillimetreFrame3d {
        attribute :>> kind = FrameKind::PFCS_FRAME;
        attribute :>> transformation : AffineTransformationMatrix3d {
            attribute :>> source = mcs;                       // [p]_MCS = M·[p]_PFCS
        }
    }
    attribute dgscs : MillimetreFrame3d {
        attribute :>> kind = FrameKind::DGSCS_FRAME;
        attribute :>> transformation : AffineTransformationMatrix3d {
            attribute :>> source = pfcs;                      // [p]_PFCS = M·[p]_DGSCS
        }
    }
    attribute dgslcs : MillimetreFrame3d {
        attribute :>> kind = FrameKind::DGSLCS_FRAME;
        attribute :>> transformation : AffineTransformationMatrix3d {
            attribute :>> source = dgscs;                     // [p]_DGSCS = M·[p]_DGSLCS
        }
    }

    // ★ D6: 预构建矩阵 (只算一次), 供按元素比较——避免整体 == 把 source/target 卷入
    private attribute pfcsXformBuilt : AffineTransformationMatrix3d =
        buildTransform_PFCS_to_MCS(
            dgsInstallCalibration.alpha_install,
            dgsInstallCalibration.Lopx, dgsInstallCalibration.Lopy, dgsInstallCalibration.Lopz);
    private attribute dgscsXformBuilt : AffineTransformationMatrix3d =
        buildTransform_DGSCS_to_PFCS(
            dgsInstallCalibration.alpha_install, dgsInstallCalibration.beta_install,
            dgsInstallCalibration.Lopx, dgsInstallCalibration.Lopy, dgsInstallCalibration.Lopz,
            dgsInstallCalibration.Losx, dgsInstallCalibration.Losy, dgsInstallCalibration.Losz);

    
    assert constraint wire_pfcs_transformation {     // ★ D6: 只比较 R/T 元素
        pfcs.transformation.rotationMatrix.elements    == pfcsXformBuilt.rotationMatrix.elements and
        pfcs.transformation.translationVector.elements == pfcsXformBuilt.translationVector.elements
    }
    assert constraint wire_dgscs_transformation {    // ★ D6
        dgscs.transformation.rotationMatrix.elements    == dgscsXformBuilt.rotationMatrix.elements and
        dgscs.transformation.translationVector.elements == dgscsXformBuilt.translationVector.elements
    }

    // =================================================================
    //  高层: DGSLCS → PFCS = compose(DGSCS_to_PFCS, DGSLCS_to_DGSCS)
    //   (需求在 PFCS 评估；PFCS→MCS 树边仅用于需要机器系坐标时)
    // =================================================================  [H 结构]
    calc def buildTransform_DGSLCS_to_PFCS {
        in alpha_Rz : PlaneAngleValue;
        in beta_Ry  : PlaneAngleValue;
        in gamma_Rx : PlaneAngleValue;
        in tx : LengthValue;  in ty : LengthValue;  in tz : LengthValue;
        in calib : DGSInstallationCalibration;

        private attribute T_sl : AffineTransformationMatrix3d =
            buildTransform_DGSLCS_to_DGSCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz);
        private attribute T_ps : AffineTransformationMatrix3d =
            buildTransform_DGSCS_to_PFCS(
                calib.alpha_install, calib.beta_install,
                calib.Lopx, calib.Lopy, calib.Lopz,
                calib.Losx, calib.Losy, calib.Losz);

        return result : AffineTransformationMatrix3d = composeTransforms(T_ps, T_sl);
    }

    calc def transformPoint_DGSLCS_to_PFCS {
        in pointInDGSLCS : Position3D;
        in alpha_Rz : PlaneAngleValue;  in beta_Ry : PlaneAngleValue;  in gamma_Rx : PlaneAngleValue;
        in tx : LengthValue;  in ty : LengthValue;  in tz : LengthValue;
        in calib : DGSInstallationCalibration;

        private attribute T : AffineTransformationMatrix3d =
            buildTransform_DGSLCS_to_PFCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz, calib);
        return result : Position3D = applyTransform(T, pointInDGSLCS);
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
    
    attribute ambientEnvironment : EnvironmentCondition {
        :>> ambientPressure    = 1.0e-5 [Pa];  // 典型 EUV 舱体真空度，按实际填
    }    
       
    item def SpatialReferencePoint {
        attribute position : Position3D;
        attribute coordinateFrame    : FrameKind default = FrameKind::PFCS_FRAME;   // 默认在 PFCS 下定义
    }

    item <PF> 'Primary Focus' : SpatialReferencePoint {
        :>> position {
            :>> x = 0 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> coordinateFrame = FrameKind::PFCS_FRAME;   // PF 是 PFCS 的原点
    }
    
    item monitorPosition : SpatialReferencePoint {
        :>> position {
            :>> x = 5 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> coordinateFrame = FrameKind::PFCS_FRAME;   // FDS 监测参考点, 在 PFCS 下定义
    }
    
    item nominalNozzlePosition : SpatialReferencePoint {
        :>> position {
            :>> x = 410.453510 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> coordinateFrame = FrameKind::PFCS_FRAME;   // 喷嘴出口标称位置, 在 PFCS 下定义
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

        // ----  喷嘴在 DGSLCS 下的标称几何 ----
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

        // ★ FIX⑦: 原 `:>> ISQ::mass` 为非法重定义 (DG 并不继承 ISQ::mass).
        attribute weight : MassValue;
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

        // ---- 调节分辨率属性 ----
        attribute translationResolution        : LengthValue;        // 平移轴最小可分辨步长
        attribute rotationResolution           : PlaneAngleValue;    // 旋转轴最小可分辨步长

        // ---- 承重 / 承压 / 自身重量属性 ----
        attribute loadCapacity                 : MassValue;          // 可承载负载质量 (kg)

        //   真实载荷是放气/抽真空瞬态下密封结构内外两侧的压差 ΔP, 见 Req-5.
        attribute differentialPressureCapacity : PressureValue;      // 可承受的内外压差上限 ΔP (Pa)
        attribute stageOwnWeight               : MassValue;          // DGS 自身质量 (kg)

        // ---- 坐标系标定与当前变换 ----
        ref attribute installCalib = dgsInstallCalibration;

        // 当前变换 (v5.0: 改用 AffineTransformationMatrix3d)
        attribute currentTransform_DGSLCS_to_DGSCS : AffineTransformationMatrix3d;
        attribute currentTransform_DGSLCS_to_PFCS  : AffineTransformationMatrix3d;

        // ★ D6: 预构建当前命令位姿对应的矩阵 (只算一次)
        private attribute slXformBuilt : AffineTransformationMatrix3d =
            buildTransform_DGSLCS_to_DGSCS(
                commandedYaw, commandedPitch, commandedRoll,
                commandedTx, commandedTy, commandedTz);
        private attribute plXformBuilt : AffineTransformationMatrix3d =
            buildTransform_DGSLCS_to_PFCS(
                commandedYaw, commandedPitch, commandedRoll,
                commandedTx, commandedTy, commandedTz,
                installCalib);
        
        // 变换一致性约束: 当前变换必须与命令位姿吻合
        assert constraint transformConsistency_sl {     // ★ D6
            currentTransform_DGSLCS_to_DGSCS.rotationMatrix.elements    == slXformBuilt.rotationMatrix.elements and
            currentTransform_DGSLCS_to_DGSCS.translationVector.elements == slXformBuilt.translationVector.elements
        }
        assert constraint transformConsistency_pl {     // ★ D6
            currentTransform_DGSLCS_to_PFCS.rotationMatrix.elements    == plXformBuilt.rotationMatrix.elements and
            currentTransform_DGSLCS_to_PFCS.translationVector.elements == plXformBuilt.translationVector.elements
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
            )
        }
        
        // ============================================================
        // ★ FIX⑨: 原嵌于此处的 Req-1 ~ Req-6 已提升至包级 (见 DGS 定义之后),
        //   约束统一经 subject 引用, 并在 DGSS 中以 satisfy 绑定主体.
        // ============================================================
    }    

    // ============================================================
    // ★ FIX⑨: DGS 性能/结构需求 (包级; 约束一律经 subject `dgs.` 引用,
    //   消除 v5.1 中 "subject 形同虚设、约束按词法直引外围特征" 的问题;
    //   主体绑定见 DGSS 中的 satisfy 声明)
    // ============================================================

    // ---- Req-1: 可调区域覆盖 ----
    requirement <'Req-1'> dropletAdjustableRangeRequirement {
        doc /* 在过 PF 点的 YZ 平面上, 以 PF 点为中心, 液滴落点应在
             * 直径不小于 15 mm 的圆内任意可调. */
        subject dgs : DGS;
        assume constraint regionGeometryMatch {
            dgs.achievableAdjustableRegionCenter == dgs.requiredAdjustableRegionCenter and
            dgs.achievableAdjustableRegionPlane  == dgs.requiredAdjustableRegionPlane  and
            dgs.achievableAdjustableRegionShape  == dgs.requiredAdjustableRegionShape
        }
        require constraint diameterCoverage {
            dgs.achievableAdjustableRegionDiameter >= dgs.requiredAdjustableRegionDiameter
        }
    }

    // ---- Req-2: 闭环带宽 ----
    requirement <'Req-2'> closedLoopBandwidthRequirement {
        doc /* DGS 各受控自由度闭环控制带宽 ≥ 5 Hz, 以 -3dB 幅值衰减为判定准则. */
        subject dgs : DGS;
        require constraint bandwidthLowerBound {
            dgs.closedLoopBandwidth >= 5 [Hz]
        }
    }

    // ---- Req-3: DGS 调节分辨率需求 ----
    requirement <'Req-3'> adjustmentResolutionRequirement {
        doc /* DGS 平移自由度 (Tx, Ty, Tz) 的调节分辨率应不大于 1 μm,
             *   即最小可控/可分辨步长 ≤ 1 μm.
             *   测试方法: 在闭环控制下连续发送小步长指令, 用外部高分辨率
             *   位移传感器 (干涉仪) 监测实际位移, 验证步长可重复分辨. */
        subject dgs : DGS;
        attribute requiredTranslationResolution : LengthValue = 1 [um];
        require constraint translationResolutionBound {
            dgs.translationResolution <= requiredTranslationResolution
        }
    }

    // ---- Req-4: DGS 承重能力需求 ----
    requirement <'Req-4'> loadCapacityRequirement {
        doc /* DGS 承重能力应不小于 25 kg, 用于支撑 Droplet Generator 本体
             *   及附属管路/电缆. */
        subject dgs : DGS;
        attribute requiredLoadCapacity : MassValue = 25 [kg];
        require constraint loadCapacityBound {
            dgs.loadCapacity >= requiredLoadCapacity
        }
    }

    // ---- Req-5: DGS 抗压差 (耐放气 / 抽真空) 能力需求 ----
    requirement <'Req-5'> differentialPressureRequirement {
        doc /* DGS 结构可承受的内外压差 ΔP 应不小于 1×10^5 Pa (约一个大气压).
             *
             *   语义澄清: DGS 的名义工作环境是高真空 (见 ambientEnvironment,
             *   约 1×10^-5 Pa). 本需求并非"工作环境绝对压力上限", 而是结构在
             *   以下瞬态工况下两侧出现的压差承受能力:
             *     - 放气 (vent):       舱体由高真空回充至大气 (~1×10^5 Pa);
             *     - 抽真空 (pump-down): 舱体由大气抽降至高真空;
             *   上述过程中, DGS 的密封面 / 夹层 / 波纹管 / 法兰 / 焊缝 两侧可能
             *   承受接近一个大气压的压差, 结构不得发生屈服、永久变形或泄漏.
             *   定义 ΔP = |p_inside − p_outside|, 取全工况包络的最大值.
             *
             *   测试方法: 压差/检漏试验 (如氦质谱检漏 + 充压保持), 在最大 ΔP
             *   下保持规定时长, 验证无永久变形且漏率在允许范围内. */
        subject dgs : DGS;
        attribute requiredDifferentialPressure : PressureValue = 1.0e5 [Pa];
        require constraint differentialPressureBound {
            dgs.differentialPressureCapacity >= requiredDifferentialPressure
        }
    }

    // ---- Req-6: DGS 自身重量约束 ----
    requirement <'Req-6'> stageWeightRequirement {
        doc /* DGS 自身重量应小于 40 kg, 以满足上级结构承载与装配限制. */
        subject dgs : DGS;
        attribute maxStageWeight : MassValue = 40 [kg];
        require constraint stageWeightBound {
            dgs.stageOwnWeight < maxStageWeight
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
        // ★ FIX: 以下 doc 原以独立成员形式出现, 实际注解的是外层 use case 而非
        //   其前的 action; 现移入各 action 体内, 与状态机内的写法保持一致.
        then action setupSensor {
            doc /* 接入共焦位移传感器 (confocal displacement sensor) 作为位置反馈装置,
                 * 确保传感器采样率远高于预期带宽 (建议 ≥ 10 倍), 并完成零位标定。 */
        }
        then action setOperatingPoint {
            doc /* 将 DGS 置于额定工作点 (零位附近), 锁定其他自由度,
                 * 仅保留被测自由度处于闭环控制状态。 */
        }
        then action injectStepCommand {
            doc /* 向被测自由度发送小幅阶跃指令:
                 *   - 阶跃幅度: 例如 10 μm (在线性区内, 远小于行程极限)
                 *   - 阶跃形式: 理想阶跃 (上升沿尽量陡峭)
                 *   - 起始点: 当前稳态工作点
                 *   - 保持时间: 足够长以观察完整的瞬态过程并达到新稳态。 */
        }
        then action recordStepResponse {
            doc /* 用共焦位移传感器同步采集:
                 *   - 阶跃指令信号 u(t)
                 *   - 实际位移响应曲线 y(t)
                 * 采样率应足以分辨上升沿细节 (建议 ≥ 1 kHz)。 */
        }
        then action extractRiseTime {
            doc /* 从响应曲线 y(t) 中提取上升时间 tr:
                 *   - 确定稳态终值 y_final (取响应后段的平均值)
                 *   - 找到 y(t) 首次到达 10% × y_final 的时刻 t10
                 *   - 找到 y(t) 首次到达 90% × y_final 的时刻 t90
                 *   - 上升时间 tr = t90 - t10。 */
        }
        then action estimateBandwidth {
            doc /* 按一阶系统近似公式估算闭环带宽:
                 *      BW ≈ 0.35 / tr
                 * 对每个受控自由度分别计算; 系统带宽取所有自由度中的最小值。 */
        }
        then action verifyAgainstRequirement {
            doc /* 与 Req2 (≥ 5 Hz) 比对; 同时记录:
                 *   - 超调量 (overshoot)
                 *   - 稳态误差 (steady-state error)
                 *   - 调节时间 (settling time)
                 * 作为系统动态性能的辅助评估指标。 */
        }
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
        attribute stability          : DropletStability;
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
        
        // ★ FIX: 原嵌于此处的 calc def getIdealHitLocalZ 已提升至包级
        //   (DataType 内嵌 calc def 的工具兼容性存疑, 且与 pointInsideTCAperture
        //    的放置位置不一致), 见下方包级定义, 并已落地为可求值表达式.
        
        // 判定点 (yp, zp) 是否落在开口内
        // 假设 centerPosition 处 y=0, z=0 为半圆/矩形交界中心
        //   上半圆区域: z >= 0 且 y² + z² <= R²
        //   下矩形区域: -H <= z <= 0 且 |y| <= W/2
    }
    
    // ---- 理想落点局部 Z 坐标: -H + offset = -20.8873 + 14.9 = -5.9873 [mm] ----
    calc def getIdealHitLocalZ {
        in geo : TCApertureGeometry;
        return result : LengthValue =
            geo.idealHitOffsetFromBottom - geo.lowerRectangleHeight;
    }        
    
    // ★ FIX: 落地为可求值表达式 (原仅注释级占位)
    calc def pointInsideTCAperture {
        in yp : LengthValue;
        in zp : LengthValue;
        in geo : TCApertureGeometry;
        return result : Boolean =
            ( zp >= 0 [mm] and
              (yp*yp + zp*zp) <= geo.upperSemicircleRadius * geo.upperSemicircleRadius )
            or
            ( zp <  0 [mm] and zp >= -geo.lowerRectangleHeight
                           and yp >= -geo.lowerRectangleWidth/2
                           and yp <=  geo.lowerRectangleWidth/2 );
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
        attribute hasCollected      : Boolean := false;   // ★ FIX④: 运行期可变状态用初值 `:=`
    }

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

        port feedbackPort : ~ControlSignalPort;
        // [设计注记] ~ControlSignalPort 意味着 FDS 直接输出 DGSCommand, 即控制律
        //   内嵌于 FDS. 若实际架构中控制器位于 DGS 侧, 建议 FDS 改为输出测量量
        //   (DropletStability), 由 DGS 内部 calc 生成命令. 留待 v6 架构决策.

        // ★ FIX⑨: 原嵌于此处的 Req-FDS-1 已提升至包级 (见 FDS 定义之后).
    }    

    // ---- Req-FDS-1: FDS 空间测量分辨率需求 (包级, ★ FIX⑨) ----
    requirement <'Req-FDS-1'> stabilityMeasurementAccuracy {
        doc /* FDS 在 Y / Z 方向的空间测量分辨率应不大于 1 μm,
             *   以保证可分辨 15 μm/3σ 的稳定性阈值. */
        subject theFds : FDS;
        require constraint spatialResolutionBound {
            theFds.detection.spatialResolution <= 1 [um]
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
        attribute roll  : ISQ::PlaneAngleValue;
        attribute pitch : ISQ::PlaneAngleValue;
        attribute yaw   : ISQ::PlaneAngleValue;
    }
      
    part def DGSS {
        part dgs    : DGS;
        part body   : DG;
        part nozzle : Nozzle;

        part shroud : Shroud;
        part inlet  : InletHole;

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
        
        // ============================================================
        // ★ D3 / FIX: 端口互连. 物项交换 (液滴/指令) 升级为 flow,
        //   显式表达传递方向与载荷; 机械安装关系保留为 connect.
        // ============================================================
        flow from nozzle.outletPort.dropletFlow to tc.inletPort.dropletFlow;   // 液滴流: out Droplet → in Droplet
        flow from fds.feedbackPort.command     to dgs.controlPort.command;     // 闭环反馈: FDS(out command) → DGS(in command)
        connect body.nozzleMount  to nozzle.mountingPort;  // 机械: DG 本体 ↔ Nozzle
        connect dgs.loadPort      to body.stageMount;      // 机械: DGS 承载口 ↔ DG 台架口
        // 注: Shroud.mountingPort / InletHole 当前无明确配对端口, 暂未连        

        // ============================================================
        // ★ FIX⑨: 需求满足声明 — 把包级需求绑定到本系统的具体主体
        // ============================================================
        satisfy requirement dropletAdjustableRangeRequirement by dgs;
        satisfy requirement closedLoopBandwidthRequirement    by dgs;
        satisfy requirement adjustmentResolutionRequirement   by dgs;
        satisfy requirement loadCapacityRequirement           by dgs;
        satisfy requirement differentialPressureRequirement   by dgs;
        satisfy requirement stageWeightRequirement            by dgs;
        satisfy requirement stabilityMeasurementAccuracy      by fds;
        
        ref item targetRef   = PF;
        ref item monitorRef  = monitorPosition;
        ref attribute envRef = ambientEnvironment;

        // ---- 当前运行状态 (供外部观察) ----
        attribute currentState : DGSSOperatingState;
        // ★ FIX④ (修正 FIX③ 回归): 以下三个量在运行期被 assign 改写, 故必须用
        //   初值 `:=` 声明. v5.1 误用绑定 `=` —— FeatureValue 绑定是贯穿生命周期
        //   的等值约束, 声明语义下 uncollectableDetected 被钉死为 false, 导致
        //   faultFlag(= uncollectableDetected) 恒假, Fault 再次不可达.
        // 不可收集检测锁存. 初值 false; 运行态由 do 动作内的事件驱动循环置位
        //   (见 monitorStability/trackReference), Fault 复位时由 faultToInit 的
        //   do 动作清零. 置位语义: pointInsideTCAperture(...)==false ∨ TC-miss≥N ∨
        //   FDS-loss ∨ 超行程, 这些条件经 FDS/几何判定汇入 lastDropletMissedTC.
        attribute uncollectableDetected : Boolean := false;
        // FDS/几何判定输出 — 最近一颗液滴 TC 穿越点是否落在开口外 (外部检测输入信号).
        attribute lastDropletMissedTC : Boolean := false;
        // 操作员手动复位指令 (外部输入; 见 Safety-1 应对措施第 7 步). 触发 Fault→Initialization 并清锁存.
        attribute resetCommand : Boolean := false;
        // ★ D4: faultFlag 绑定到检测锁存. 此处保留绑定 `=` —— faultFlag 是纯派生
        //   别名, 自身从不被 assign, 其值随 uncollectableDetected 变化.
        attribute faultFlag    : Boolean = uncollectableDetected;

        // ============================================================
        // ★ DGSS 状态机
        //   States: Initialization → Idle → ClosedLoopLocking
        //   Any state → Fault (异常发生)
        //   Fault → Initialization (复位后重新初始化)
        // ★ FIX⑥: 改为 exhibit state, 确保该状态行为由 DGSS 实例履行
        //   (裸 state 成员仅声明特征, 不保证被宿主表现).
        // ★ FIX⑤: 仅含 if 守卫的迁移改为 accept when 变化触发. 原因: 无触发器的
        //   迁移仅在源状态"完成"时被使能, 而 Idle/ClosedLoopLocking 含无界 do
        //   动作永不完成 → 守卫永无评估时机; Fault 虽可完成, 但完成瞬间若
        //   resetCommand 为假, 迁移错过后再无重评机制. 改用变化触发后, 条件
        //   变为真的时刻即触发迁移.
        // ============================================================
        exhibit state DGSSStateMachine {
            entry; then Initialization;

            // ---------------- Initialization ----------------
            state Initialization {
                // ★ FIX①②: 统一为单个复合 entry 动作 (与 Idle/Locking 风格一致), doc 内置于各 action
                entry action initSequence {
                    first start;
                    then action homeDGS {
                        doc /* DGS 回零: 所有平移/旋转轴执行 homing 序列, 确认编码器/限位开关状态. (行为占位) */
                        assign currentState := DGSSOperatingState::initialization;     // ★ D5
                    }
                    then action moveDropletToFDSRef {
                        doc /* 控制 DGS 使液滴落点对准 FDS 参考位置 (monitorPosition), 完成粗对准. */
                    }
                    then action verifyInitComplete {
                        doc /* 校验所有轴位置在期望区间内, 无故障标志. */
                    }
                    then done;
                }
            }
            // 初始化完成 → Idle
            transition initToIdle
                first Initialization
                accept after (1 [s])      // 占位: 实际由 verifyInitComplete 完成事件触发
                if (not faultFlag)
                then Idle;
            // 初始化失败 → Fault   (★ FIX⑤: 改为变化触发)
            transition initToFault
                first Initialization
                accept when faultFlag
                then Fault;

            // ---------------- Idle ----------------
            state Idle {
                doc /* 初始化完成或从 ClosedLoopLocking 退出后的待命状态.
                     *   液滴可能正在喷射, 但 DGS 控制环未参与闭环锁定. */
                entry action openControlLoop {
                    doc /* 释放闭环锁定 (若曾激活), DGS 维持当前姿态. */
                    assign currentState := DGSSOperatingState::idle;               // ★ D5
                }
                do action monitorStability {
                    doc /* 持续从 FDS 读取液滴稳定性测量结果, 但不闭环修正. */
                    // ★ FIX⑧: 锁存改为变化事件驱动的循环. v5.1 的单次 assign 仅在
                    //   进入状态的瞬间采样一次, 不具备"持续锁存"语义.
                    loop action latchOnMissIdle {
                        accept when lastDropletMissedTC;
                        then assign uncollectableDetected := true;
                    }
                }
            }
            // Idle → ClosedLoopLocking : 仅当 Y/Z 横向稳定性达标   (★ FIX⑤)
            transition idleToLocking
                first Idle
                accept when (fds.measuredStability.isTransverseStable and not faultFlag)
                then ClosedLoopLocking;
            // Idle → Fault   (★ FIX⑤)
            transition idleToFault
                first Idle
                accept when faultFlag
                then Fault;

            // ---------------- ClosedLoopLocking ----------------
            state ClosedLoopLocking {
                doc /* 利用 FDS 实时反馈, 闭环调整 DGS, 使液滴稳定锁定在
                     *   FDS 参考位置 (monitorPosition). */
                entry action engageClosedLoop {
                    doc /* 接入 FDS → DGS 闭环路径, 开启位置闭环控制器. */
                    assign currentState := DGSSOperatingState::closedLoopLocking;  // ★ D5
                }                    
                do action trackReference {
                    doc /* 实时计算液滴位置误差, 输出 DGS 修正指令. */
                    // ★ FIX⑧: 同 monitorStability, 事件驱动锁存
                    loop action latchOnMissLocking {
                        accept when lastDropletMissedTC;
                        then assign uncollectableDetected := true;
                    }
                }
            }
            // ClosedLoopLocking → Idle : 横向稳定性失守 (★ FIX⑤)
            transition lockingToIdle
                first ClosedLoopLocking
                accept when ((not fds.measuredStability.isTransverseStable) and (not faultFlag))
                then Idle;
            // ClosedLoopLocking → Fault   (★ FIX⑤)
            transition lockingToFault
                first ClosedLoopLocking
                accept when faultFlag
                then Fault;

            // ---------------- Fault ----------------
            state Fault {
                doc /* 故障状态: 检测到不可恢复异常 (运动学超限/传感器丢失/
                     *   液滴丢失/TC 收集失败超阈值/压差超限 等).
                     *   触发安全响应序列 (见 SafetyConstraint). */
                entry action raiseAlarm {                                          // ★ D5
                    assign currentState := DGSSOperatingState::fault;
                }                
                //entry action executeSafetyResponse;
                //    doc /* 见 dropletCannotBeCollectedSafetyResponse 定义的应对措施. */
            }
            // Fault → Initialization : 外部复位 (★ FIX⑤: 复位指令到来即触发)
            transition faultToInit
                first Fault
                accept when resetCommand
                do action clearFaultLatch {
                    assign uncollectableDetected := false;     // 清除故障锁存
                    assign resetCommand          := false;     // 复位指令自清除 (单次触发)
                }
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

    // 注: distancePointToLine 适用于直线轨迹; 自 v4.3 起轨迹为抛物线
    //   (BallisticTrajectory), 穿越判据已统一改用 crossingPointAtX. 该占位函数
    //   保留以兼容历史引用, 不再被穿越类需求使用.
    calc def distancePointToLine {
        in point : Position3D;
        in line : BallisticTrajectory;
        return : LengthValue;
    }   
     
    requirement def DropletPassesTarget {
        doc /* 液滴轨迹必须经过 PF (0,0,0).
             *   [意见2] 改为基于抛物线弹道在 x = target.x 平面处的穿越点偏差判据,
             *   与 BallisticTrajectory 自洽 (不再用"点到直线距离"). */
        subject d : Droplet;
        ref item target = PF;
        require constraint {
            // 允许微小的物理容差 (< 1 微米, 按实际精度定义)
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, target.position.x).y - target.position.y
            ) <= 0.001 [mm] and
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, target.position.x).z - target.position.z
            ) <= 0.001 [mm]
        }
    }
    
    requirement def DropletPassesMonitor {
        doc /* 液滴轨迹必须经过 FDS reference position (5,0,0).
             *   [意见2] 同 DropletPassesTarget, 改用 crossingPointAtX 穿越点偏差. */
        subject d : Droplet;
        ref item monitor = monitorPosition;
        require constraint {
            // 允许微小的物理容差 (< 1 微米, 按实际精度定义)
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, monitor.position.x).y - monitor.position.y
            ) <= 0.001 [mm] and
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, monitor.position.x).z - monitor.position.z
            ) <= 0.001 [mm]
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
            )
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
            system.currentState == DGSSOperatingState::fault or
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
    // 坐标变换用例 — DGS 命令位姿求解
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

        // ★ FIX: doc 移入各 action 体内 (原以独立成员形式注解了外层 use case)
        then action loadCalibration {
            doc /* 加载 dgsInstallCalibration: α=28°, β=78°,
                 *   Lpsx=916.328834 mm, Lpsy=43.549721 mm, Lpsz=81.905114 mm. */
        }
        then action specifyTargetInPFCS {
            doc /* 指定 PF 系下的目标点 P_p. 典型取值: PF (0,0,0) 或 monitorPosition (5,0,0). */
        }
        then action computeT_ps {
            doc /* 由标定参数构建 DGSCS → PFCS 固定变换 T_ps. */
        }
        then action computeRequiredT_sl {
            doc /* 求解 T_sl: 利用 invertTransform(T_ps) 得 T_sp,
                 *   并由 T_sp · P_p 给出喷嘴出口在 DGSCS 下应到达的位置 P_s.
                 *   再结合 P_l (喷嘴出口在 DGSLCS 下的标称位置) 解 T_sl. */
        }
        then action extractEulerAndTranslation {
            doc /* 从 T_sl 的旋转 3×3 子矩阵按 ZYX 反解 (α, β, γ),
                 *   并从第 4 列得 (tx, ty, tz). */
        }
        then action verifyWithinTravelLimits {
            doc /* 校验 (tx, ty, tz, α, β, γ) 满足 DGS 的 singleAxisTravelLimits
                 *   与 coupledWorkspaceReachability. 若不满足则报告反解失败. */
        }
        then action issueDGSCommand {
            doc /* 把反解结果打包为 DGSCommand 通过 controlPort 下发. */
        }
    }    
}