package DGSS_Model {

    doc /* Version: 5.0  (基于 4.3)
         * Author: Simon  (变换层优化 by 协作迭代)
         * Date: 2026-06-08
         * Change history (相对 v4.3):
         *   1. 删除手搓 HomogeneousTransform / Vector4，统一改用标准库
         *      MeasurementReferences::AffineTransformationMatrix3d
         *      (自动补底行 (0,0,0,1)；分解 rotationMatrix:Real[9] + translationVector)。
         *   2. 旋转矩阵 row-major (库 elements 派生式钉死)；语义 [p]_source = M·[p]_target。
         *   3. buildTransform_DGSLCS_to_DGSCS：角度先 /(1[radian]) 消量纲再喂 sin/cos；
         *      R = Rz(α)·Ry(β)·Rx(γ) 采用验证过的正确符号。
         *   4. 原为注释桩的所有 buildTransform_* / applyTransform / compose / invert 全部实现。
         *   5. 四坐标系升级为 ThreeDCoordinateFrame，按链 MCS(根)⊃PFCS⊃DGSCS⊃DGSLCS
         *      用 transformation.source 固化安置关系 (target=that 由库保证)。
         *   6. 新增 buildTransform_PFCS_to_MCS (R=Rx(α), T=L_op)，即 PFCS→MCS 树边。
         *   7. 平移元素统一 /(1[mm])，与 MillimetreFrame3d 的 mRefs 对齐。
         *   8. cos/sin 作用于 PlaneAngleValue 处统一加 /(1[radian]) 消量纲。
         *
         * 复核标注：[H] 高；[M] 需在 Syside/SymPy 复核 (4×4 展开/单位名/绑定语法)。
         */

    private import SI::*;
    private import ISQ::*;
    private import ScalarValues::*;
    private import MeasurementReferences::*;   // 替换 v4.3 的 SpatialFrames::*
    private import ScalarFunctions::*;
    private import TrigFunctions::*;
    private import VectorValues::ThreeVectorValue;
    private import NumericalFunctions::*;

    attribute <um> micrometer : LengthUnit = micro * metre;

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
        DGSLCS_FRAME;
        DGSCS_FRAME;
        MCS_FRAME;
        PFCS_FRAME;
    }

    enum def DGSSOperatingState {
        initialization;
        idle;
        closedLoopLocking;
        fault;
    }

    attribute def Position3DInFrame :> Position3D {
        doc /* 带坐标系标签的 3D 位置. */
        attribute coordinateFrame : FrameKind;
    }

    attribute def DisplacementVector {
        doc /* 位移向量 L_op = p - o；coordinateFrame 标明所在坐标系. */
        attribute x               : LengthValue;
        attribute y               : LengthValue;
        attribute z               : LengthValue;
        attribute coordinateFrame : FrameKind;
    }

    calc def displacementVector {
        in originPoint : Position3D;
        in targetPoint : Position3D;
        in coordinateFrame : FrameKind;
        return : DisplacementVector;
    }

    attribute def RotationSignConvention {
        doc /* 右手定则角度正方向约定. */
        attribute isRightHanded : Boolean = true;
        attribute description   : String  =
            "Right-hand rule: thumb along +axis, curled fingers give positive angle direction.";
    }

    item globalRotationConvention : RotationSignConvention;

    // =================================================================
    //  坐标系类型: 三轴单位均为 mm 的 3D 正交坐标系
    //  [M] millimetre 的确切元素名请核 SI 库 (可能是 millimetre 或 mm 别名)
    // =================================================================
    attribute def MillimetreFrame3d :> ThreeDCoordinateFrame {
        attribute :>> mRefs = (millimetre, millimetre, millimetre);
        attribute :>> isOrthogonal = true;
    }

    // =================================================================
    //  通用算子
    // =================================================================

    // ---- 应用变换: [p]_source = M · [p]_target ----  [H]
    calc def applyTransform {
        in M : AffineTransformationMatrix3d;
        in p : Position3D;
        return : Position3D;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
        return : AffineTransformationMatrix3d;

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
    // 安装标定 — 全局只读常量
    //   L_op, L_os 在 MCS; L_ps 在 PFCS.
    // ============================================================
    attribute def DGSInstallationCalibration {
        doc /* DGS 安装姿态与平移标定. α 连接 PFCS↔MCS, β 连接 DGSCS↔MCS.
             *   旋转角遵循右手定则. */
        attribute alpha_install : PlaneAngleValue = 28 ['°'];   // PFCS ↔ MCS
        attribute beta_install  : PlaneAngleValue = 78 ['°'];   // DGSCS ↔ MCS

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

        attribute L_op : DisplacementVector {    // L_op = P - o (MCS)
            :>> x =    0.000000 [mm];
            :>> y =  -41.524524 [mm];
            :>> z = 1248.672688 [mm];
            :>> coordinateFrame = FrameKind::MCS_FRAME;
        }
        attribute L_os : DisplacementVector {    // L_os = S - o (MCS)
            :>> x =  916.328834 [mm];
            :>> y =  -41.524524 [mm];
            :>> z = 1341.435967 [mm];
            :>> coordinateFrame = FrameKind::MCS_FRAME;
        }
        attribute L_ps : DisplacementVector {    // L_ps = S - P (PFCS) = DGSCS→PFCS 平移列
            :>> x =  916.328834 [mm];
            :>> y =   43.549721 [mm];
            :>> z =   81.905114 [mm];
            :>> coordinateFrame = FrameKind::PFCS_FRAME;
        }

        attribute Lopx : LengthValue = L_op.x;   attribute Lopy : LengthValue = L_op.y;   attribute Lopz : LengthValue = L_op.z;
        attribute Losx : LengthValue = L_os.x;   attribute Losy : LengthValue = L_os.y;   attribute Losz : LengthValue = L_os.z;
        attribute Lpsx : LengthValue = L_ps.x;   attribute Lpsy : LengthValue = L_ps.y;   attribute Lpsz : LengthValue = L_ps.z;

        assert constraint displacementDefinitionConsistency {
            L_op.x == P_in_MCS.x and L_op.y == P_in_MCS.y and L_op.z == P_in_MCS.z and
            L_os.x == S_in_MCS.x and L_os.y == S_in_MCS.y and L_os.z == S_in_MCS.z
        }

        // ---- 标定一致性 (cos/sin 消量纲) ----  [H]
        assert constraint translationConsistency {
            Lpsx == (Losx - Lopx) and
            Lpsy == (Losy - Lopy) * cos(alpha_install/(1 [radian])) + (Losz - Lopz) * sin(alpha_install/(1 [radian])) and
            Lpsz == (Lopy - Losy) * sin(alpha_install/(1 [radian])) + (Losz - Lopz) * cos(alpha_install/(1 [radian]))
        }
    }

    item dgsInstallCalibration : DGSInstallationCalibration;

    // =================================================================
    //  四坐标系实例 + 安置关系 (新链: MCS 根)
    //  [M] frame↔矩阵的等值绑定语法请在工具中复核；若 == 整体比较因
    //      source/target 参与而失败，改为只绑 rotationMatrix/translationVector。
    // =================================================================
    attribute mcs : MillimetreFrame3d;                        // 根参考系，无 transformation

    attribute pfcs : MillimetreFrame3d {
        attribute :>> transformation : AffineTransformationMatrix3d {
            attribute :>> source = mcs;                       // [p]_MCS = M·[p]_PFCS
        }
    }
    attribute dgscs : MillimetreFrame3d {
        attribute :>> transformation : AffineTransformationMatrix3d {
            attribute :>> source = pfcs;                      // [p]_PFCS = M·[p]_DGSCS
        }
    }
    attribute dgslcs : MillimetreFrame3d {
        attribute :>> transformation : AffineTransformationMatrix3d {
            attribute :>> source = dgscs;                     // [p]_DGSCS = M·[p]_DGSLCS
        }
    }

    // 静态树边绑定
    assert constraint wire_pfcs_transformation {
        pfcs.transformation == buildTransform_PFCS_to_MCS(
            dgsInstallCalibration.alpha_install,
            dgsInstallCalibration.Lopx, dgsInstallCalibration.Lopy, dgsInstallCalibration.Lopz)
    }
    assert constraint wire_dgscs_transformation {
        dgscs.transformation == buildTransform_DGSCS_to_PFCS(
            dgsInstallCalibration.alpha_install, dgsInstallCalibration.beta_install,
            dgsInstallCalibration.Lopx, dgsInstallCalibration.Lopy, dgsInstallCalibration.Lopz,
            dgsInstallCalibration.Losx, dgsInstallCalibration.Losy, dgsInstallCalibration.Losz)
    }
    //  dgslcs.transformation 为实时量；由 DGS part 的 currentTransform_DGSLCS_to_DGSCS 提供
    //  (DGS 内 transformConsistency_sl 约束已将其绑定到命令位姿)。

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
        return : AffineTransformationMatrix3d;

        private attribute T_sl : AffineTransformationMatrix3d =
            buildTransform_DGSLCS_to_DGSCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz);
        private attribute T_ps : AffineTransformationMatrix3d =
            buildTransform_DGSCS_to_PFCS(
                calib.alpha_install, calib.beta_install,
                calib.Lopx, calib.Lopy, calib.Lopz,
                calib.Losx, calib.Losy, calib.Losz);

        return composeTransforms(T_ps, T_sl);
    }

    calc def transformPoint_DGSLCS_to_PFCS {
        in pointInDGSLCS : Position3D;
        in alpha_Rz : PlaneAngleValue;  in beta_Ry : PlaneAngleValue;  in gamma_Rx : PlaneAngleValue;
        in tx : LengthValue;  in ty : LengthValue;  in tz : LengthValue;
        in calib : DGSInstallationCalibration;
        return : Position3D;

        private attribute T : AffineTransformationMatrix3d =
            buildTransform_DGSLCS_to_PFCS(alpha_Rz, beta_Ry, gamma_Rx, tx, ty, tz, calib);
        return applyTransform(T, pointInDGSLCS);
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
        :>> ambientPressure    = 1.0e-5 [Pa];
    }

    item def SpatialReferencePoint {
        attribute position : Position3D;
        attribute coordinateFrame    : FrameKind default = FrameKind::PFCS_FRAME;
    }

    item <PF> 'Primary Focus' : SpatialReferencePoint {
        :>> position {
            :>> x = 0 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> coordinateFrame = FrameKind::PFCS_FRAME;
    }

    item monitorPosition : SpatialReferencePoint {
        :>> position {
            :>> x = 5 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> coordinateFrame = FrameKind::PFCS_FRAME;
    }

    item nominalNozzlePosition : SpatialReferencePoint {
        :>> position {
            :>> x = 410.453510 [mm];
            :>> y = 0 [mm];
            :>> z = 0 [mm];
        }
        :>> coordinateFrame = FrameKind::PFCS_FRAME;
    }

    attribute def DropletStability {
        doc /* 液滴稳定性: 任意连续 30 s 窗口内, X/Y/Z 三方向 3σ < 15 μm. */
        attribute evaluationWindow      : TimeValue   = 30 [s];
        attribute sigmaMultiplier       : Real        = 3.0;
        attribute stabilityThreshold    : LengthValue = 15 [um];

        attribute longitudinalSigma3_X  : LengthValue;
        attribute transverseSigma3_Y    : LengthValue;
        attribute transverseSigma3_Z    : LengthValue;

        attribute isLongitudinalStable  : Boolean;
        attribute isTransverseStable    : Boolean;
        attribute isFullyStable         : Boolean;

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

        attribute outletPositionInDGSLCS  : Position3D;
        attribute outletDirectionInDGSLCS : Direction3D;
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

        attribute closedLoopBandwidth         : FrequencyValue;
        attribute translationResolution        : LengthValue;
        attribute rotationResolution           : PlaneAngleValue;
        attribute loadCapacity                 : MassValue;
        attribute differentialPressureCapacity : PressureValue;
        attribute stageOwnWeight               : MassValue;

        ref item installCalib = dgsInstallCalibration;

        // 当前变换 (v5.0: 改用 AffineTransformationMatrix3d)
        attribute currentTransform_DGSLCS_to_DGSCS : AffineTransformationMatrix3d;
        attribute currentTransform_DGSLCS_to_PFCS  : AffineTransformationMatrix3d;

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

        calc def reachableWorkspace {
            in tx : LengthValue;
            in ty : LengthValue;
            in tz : LengthValue;
            in rx : PlaneAngleValue;
            in ry : PlaneAngleValue;
            in rz : PlaneAngleValue;
            return : Boolean;
        }

        assert constraint singleAxisTravelLimits {
            commandedTx    >= txMin and commandedTx    <= txMax and
            commandedTy    >= tyMin and commandedTy    <= tyMax and
            commandedTz    >= tzMin and commandedTz    <= tzMax and
            commandedRoll  >= rxMin and commandedRoll  <= rxMax and
            commandedPitch >= ryMin and commandedPitch <= ryMax and
            commandedYaw   >= rzMin and commandedYaw   <= rzMax
        }

        assert constraint coupledWorkspaceReachability {
            reachableWorkspace(
                commandedTx, commandedTy, commandedTz,
                commandedRoll, commandedPitch, commandedYaw
            ) == true
        }

        requirement <'Req-1'> dropletAdjustableRangeRequirement {
            doc /* 在过 PF 的 YZ 平面上, 以 PF 为中心, 落点应在直径 ≥ 15 mm 的圆内任意可调. */
            subject dgs : DGS;
            assume constraint regionGeometryMatch {
                achievableAdjustableRegionCenter == requiredAdjustableRegionCenter and
                achievableAdjustableRegionPlane  == requiredAdjustableRegionPlane  and
                achievableAdjustableRegionShape  == requiredAdjustableRegionShape
            }
            require constraint diameterCoverage {
                achievableAdjustableRegionDiameter >= requiredAdjustableRegionDiameter
            }
        }

        requirement <'Req-2'> closedLoopBandwidthRequirement {
            doc /* 各受控自由度闭环带宽 ≥ 5 Hz (-3dB 判据). */
            subject dgs : DGS;
            require constraint bandwidthLowerBound {
                closedLoopBandwidth >= 5 [Hz]
            }
        }

        requirement <'Req-3'> adjustmentResolutionRequirement {
            doc /* 平移自由度调节分辨率 ≤ 1 μm. */
            subject dgs : DGS;
            attribute requiredTranslationResolution : LengthValue = 1 [um];
            require constraint translationResolutionBound {
                translationResolution <= requiredTranslationResolution
            }
        }

        requirement <'Req-4'> loadCapacityRequirement {
            doc /* 承重能力 ≥ 25 kg. */
            subject dgs : DGS;
            attribute requiredLoadCapacity : MassValue = 25 [kg];
            require constraint loadCapacityBound {
                loadCapacity >= requiredLoadCapacity
            }
        }

        requirement <'Req-5'> differentialPressureRequirement {
            doc /* 可承受内外压差 ΔP ≥ 1×10^5 Pa (放气/抽真空瞬态). */
            subject dgs : DGS;
            attribute requiredDifferentialPressure : PressureValue = 1.0e5 [Pa];
            require constraint differentialPressureBound {
                differentialPressureCapacity >= requiredDifferentialPressure
            }
        }

        requirement <'Req-6'> stageWeightRequirement {
            doc /* DGS 自重 < 40 kg. */
            subject dgs : DGS;
            attribute maxStageWeight : MassValue = 40 [kg];
            require constraint stageWeightBound {
                stageOwnWeight < maxStageWeight
            }
        }
    }

    use case def MeasureClosedLoopBandwidth {
        subject dgs : DGS;
        objective {
            doc /* 阶跃响应法, BW ≈ 0.35 / tr 估算闭环带宽. */
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

    attribute def BallisticTrajectory {
        doc /* 考虑重力与气流扰动的弹道轨迹. */
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
        attribute stability          : DropletStability;
    }

    attribute def ApertureGeometry {
        doc /* 通用通孔几何, 开口位于某 YZ 平面. */
        attribute shape          : String;
        attribute centerPosition : Position3D;
        attribute width          : LengthValue;
        attribute height         : LengthValue;
    }

    attribute def TCApertureGeometry {
        doc /* TC 开口: 上半圆 (R=33mm) + 下矩形 (66×20.8873mm), YZ 平面. */
        attribute openingPlane           : String       = "YZ-plane at TC";
        attribute centerPosition         : Position3D{
                :>> x = -310.0930 [mm];
                :>> y = 0 [mm];
                :>> z = 5.9873 [mm];
        }
        attribute upperSemicircleRadius  : LengthValue  = 33      [mm];
        attribute lowerRectangleWidth    : LengthValue  = 66      [mm];
        attribute lowerRectangleHeight   : LengthValue  = 20.8873 [mm];

        attribute idealHitOffsetFromBottom : LengthValue = 14.9 [mm];

        calc def getIdealHitLocalZ {
            in geo : TCApertureGeometry;
            return : LengthValue;
            // -lowerRectangleHeight + idealHitOffsetFromBottom = -20.8873 + 14.9 = -5.9873 mm
        }
    }

    calc def pointInsideTCAperture {
        in yp : LengthValue;
        in zp : LengthValue;
        in geo : TCApertureGeometry;
        return : Boolean;
    }

    calc def pointInsideAperture {
        in yp : LengthValue;
        in zp : LengthValue;
        in geo : ApertureGeometry;
        return : Boolean;
    }

    part def Shroud {
        doc /* Nozzle 与 PF 之间的防护罩, 中央通孔. */
        attribute mountingPose : Pose6DOF;
        attribute aperture     : ApertureGeometry {
            :>> shape = "circle";
        }
        port mountingPort : MechanicalMountPort;
    }

    part def InletHole {
        doc /* 收集镜内舱入口孔. */
        attribute aperture : ApertureGeometry {
            :>> shape = "circle";
        }
    }

    part def <TC> 'Tin Capture' {
        doc /* 液滴收集器, 复合开口. */
        attribute mountingPose      : Pose6DOF;
        attribute aperture          : TCApertureGeometry;
        attribute idealCollectionPosition: Position3D;
        port inletPort              : DropletInletPort;
        attribute hasCollected      : Boolean = false;
    }

    part def <FDS> 'Fine Detection System' {
        doc /* 小视场探测系统, 实时监测液滴位置/时序, 为闭环提供反馈. */

        part def IlluminationComponent {
            attribute wavelength      : LengthValue;
            attribute pulseDuration   : TimeValue;
            attribute repetitionRate  : FrequencyValue;
            attribute opticalPower    : PowerValue;
            port lightOutPort : MechanicalMountPort;
        }

        part def DetectionComponent {
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

        ref item monitorRef = monitorPosition;

        attribute measuredStability : DropletStability;

        port feedbackPort : ControlSignalPort;

        requirement <'Req-FDS-1'> stabilityMeasurementAccuracy {
            doc /* FDS Y/Z 空间分辨率 ≤ 1 μm. */
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

        attribute currentState : DGSSOperatingState;
        attribute faultFlag    : Boolean = false;

        state DGSSStateMachine {
            entry; then Initialization;

            state Initialization {
                entry action homeDGS;
                then action moveDropletToFDSRef;
                then action verifyInitComplete;
            }
            transition initToIdle
                first Initialization
                accept after (1 [s])
                if (not faultFlag)
                then Idle;
            transition initToFault
                first Initialization
                if (faultFlag)
                then Fault;

            state Idle {
                doc /* 待命: 液滴可能喷射, 但 DGS 未参与闭环锁定. */
                entry action openControlLoop;
                do action monitorStability;
            }
            transition idleToLocking
                first Idle
                if (fds.measuredStability.isTransverseStable and not faultFlag)
                then ClosedLoopLocking;
            transition idleToFault
                first Idle
                if (faultFlag)
                then Fault;

            state ClosedLoopLocking {
                doc /* 闭环锁定液滴于 FDS 参考位置. */
                entry action engageClosedLoop;
                do action trackReference;
            }
            transition lockingToIdle
                first ClosedLoopLocking
                if (
                    (not fds.measuredStability.isTransverseStable) and
                    (not faultFlag)
                )
                then Idle;
            transition lockingToFault
                first ClosedLoopLocking
                if (faultFlag)
                then Fault;

            state Fault {
                doc /* 故障状态, 触发安全响应序列. */
                entry action raiseAlarm;
                entry action executeSafetyResponse;
            }
            transition faultToInit
                first Fault
                accept after (1 [s])
                if (not faultFlag)
                then Initialization;
        }
    }

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
        doc /* 名义弹道与 YZ 平面 (x = xPlane) 的交点. */
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
        doc /* 液滴轨迹经过 PF (0,0,0)，用 crossingPointAtX 穿越点偏差判据. */
        subject d : Droplet;
        ref item target = PF;
        require constraint {
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, target.position.x).y - target.position.y
            ) <= 0.001 [mm] and
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, target.position.x).z - target.position.z
            ) <= 0.001 [mm]
        }
    }

    requirement def DropletPassesMonitor {
        doc /* 液滴轨迹经过 FDS reference position (5,0,0). */
        subject d : Droplet;
        ref item monitor = monitorPosition;
        require constraint {
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, monitor.position.x).y - monitor.position.y
            ) <= 0.001 [mm] and
            NumericalFunctions::abs(
                crossingPointAtX(d.trajectory, monitor.position.x).z - monitor.position.z
            ) <= 0.001 [mm]
        }
    }

    requirement def DropletCollectedByTC {
        doc /* 落点必须落在 TC 复合开口内. */
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

    requirement def DropletStabilityRequirement {
        doc /* 任意连续 30 s, X/Y/Z 三方向 3σ < 15 μm. */
        subject d : Droplet;
        require constraint stabilityBound {
            d.stability.longitudinalSigma3_X < d.stability.stabilityThreshold and
            d.stability.transverseSigma3_Y  < d.stability.stabilityThreshold and
            d.stability.transverseSigma3_Z  < d.stability.stabilityThreshold
        }
    }

    requirement def IdealGoal {
        doc /* 理想目标: 依次经过 monitor 与 target, 并被 TC 收集. */
        subject d : Droplet;
        requirement passMonitor : DropletPassesMonitor;
        requirement passTarget  : DropletPassesTarget;
        requirement collected   : DropletCollectedByTC;

        bind passMonitor.d = d;
        bind passTarget.d  = d;
        bind collected.d   = d;
    }

    requirement def FallbackGoal {
        doc /* 变通目标: 经过 target 并被 TC 收集. */
        subject d : Droplet;
        requirement passTarget : DropletPassesTarget;
        requirement collected  : DropletCollectedByTC;

        bind passTarget.d  = d;
        bind collected.d   = d;
    }

    requirement goalHierarchy {
        requirement ideal : IdealGoal;
    }

    requirement def DropletBallisticBudget {
        doc /* 弹道偏差预算: 名义下坠 + 气流散布. */
        subject d  : Droplet;
        ref part tcRef : TC;

        require constraint nominalDropBudget {
            NumericalFunctions::abs(d.trajectory.nominalCrossingAtTC.z - tcRef.aperture.centerPosition.z) <= tcRef.aperture.lowerRectangleHeight / 2
        }

        require constraint dispersionBudget {
            3 * d.trajectory.crossingDispersionAtTC_y
                <= tcRef.aperture.lowerRectangleWidth / 6  and
            3 * d.trajectory.crossingDispersionAtTC_z
                <= tcRef.aperture.lowerRectangleHeight / 6
        }
    }

    requirement def SafetyConstraint {
        doc /* 安全约束基类. */
        attribute hazardDescription   : String;
        attribute detectionMethod     : String;
        attribute responseAction      : String;
        attribute severityLevel       : String;
    }

    requirement <'Safety-1'> dropletCannotBeCollectedSafetyResponse : SafetyConstraint {
        doc /* 液滴无法被 TC 收集时的安全约束: 关闭液滴生成器 → 进入 Fault →
             *   抑制驱动激光 → 局部抽气 → 报警 → 记录遥测 → 手动复位. */
        :>> hazardDescription =
            "Droplet trajectory does not enter TC aperture, risking tin contamination of collector optics.";
        :>> detectionMethod   =
            "Geometric check via pointInsideTCAperture + statistical TC-miss counter + FDS loss-of-droplet detection.";
        :>> responseAction    =
            "Stop droplet generator -> enter Fault -> inhibit drive laser -> local pump-down -> alarm -> log telemetry -> manual reset to Initialization.";
        :>> severityLevel     = "Critical";

        subject system : DGSS;

        require constraint enterFaultOnUncollectable {
            system.currentState == "Fault" or
            system.faultFlag == false
        }
    }

    use case def AdjustNozzleToHitTarget {
        subject system : DGSS;
        objective {
            doc /* 控制 DGS 使液滴满足 IdealGoal; 否则至少 FallbackGoal. */
        }
        first start;
        then action measureDeviation;
        then action computeDGSCommand;
        then action applyDGSCommand;
        then action emitDroplet;
        then action evaluateGoal;
        then done;
    }

    use case def SolveDGSCommandFromPFTarget {
        subject system : DGSS;
        objective {
            doc /* 由 PFCS 目标点反解 DGS 6-DOF 命令.
                 *   正向: P_p = T_ps(α,β,L_ps) · T_sl(α_cmd,β_cmd,γ_cmd; tx,ty,tz) · P_l
                 *   反解: 用 invertTransform 求 T_sp, 由 T_sp·P_p 得 P_s, 结合 P_l 解 T_sl,
                 *         从 T_sl 按 ZYX 提取 (α,β,γ) 与第4列 (tx,ty,tz),
                 *         校验 singleAxisTravelLimits 与 coupledWorkspaceReachability. */
        }

        first start;
        then action loadCalibration;
        then action specifyTargetInPFCS;
        then action computeT_ps;
        then action computeRequiredT_sl;
        then action extractEulerAndTranslation;
        then action verifyWithinTravelLimits;
        then action issueDGSCommand;
    }
}
