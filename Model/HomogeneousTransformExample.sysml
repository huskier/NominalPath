package HomogeneousTransformExample {
    doc
    /*
     * 简单的齐次（仿射）坐标变换示例。
     *
     * 场景：让一个“物体坐标系 bodyFrame”相对“世界坐标系 worldFrame”：
     *   - 绕 Z 轴右手旋转 90°
     *   - 再沿世界 X 轴平移 2 米
     *
     * 文件给出同一放置的两种等价写法：
     *   Version A —— 齐次(仿射)变换矩阵  AffineTransformationMatrix3d
     *   Version B —— 平移/旋转序列        TranslationRotationSequence
     * 以及 Version C —— 恒等变换 NullTransformation（作对照）。
     *
     * 所用类型均来自标准库 MeasurementReferences / Quantities / SI。
     * 注：本环境无 SysML v2 求解器，无法编译验证；构造均依据库源文件，仅供参考。
     */

    private import ScalarValues::Real;
    private import SI::*;                                   // metre, degree, radian, ...
    private import MeasurementReferences::ThreeDCoordinateFrame;
    private import MeasurementReferences::VectorMeasurementReference;
    private import MeasurementReferences::AffineTransformationMatrix3d;
    private import MeasurementReferences::TranslationRotationSequence;
    private import MeasurementReferences::Translation;
    private import MeasurementReferences::Rotation;
    private import MeasurementReferences::NullTransformation;
    private import Quantities::ThreeDVectorQuantityValue;   // 别名 of '3dVectorQuantityValue'

    /* -----------------------------------------------------------------
     * 0. 一个“带单位”的 3D 坐标系：三个轴的单位都是 metre
     *    （mRefs 继承自 TensorMeasurementReference，逐轴给出测量参考）
     * ----------------------------------------------------------------- */
    attribute def MetricFrame3d :> ThreeDCoordinateFrame {
        attribute :>> mRefs = (metre, metre, metre);       // 轴单位：米
        attribute :>> isOrthogonal = true;
    }

    // 参考(源)坐标系：世界坐标系
    attribute worldFrame : MetricFrame3d;

    /* =================================================================
     * Version A —— 齐次(仿射)变换矩阵
     *   左上 3x3 = 旋转矩阵，右上 3x1 = 平移向量，底行库里自动补 (0,0,0,1)
     * ================================================================= */
    attribute bodyFrameA : MetricFrame3d {

        attribute :>> transformation : AffineTransformationMatrix3d {

            attribute :>> source = worldFrame;             // 在世界坐标系中描述
            // target 由 CoordinateFrame::transformation 绑定为 that(即 bodyFrameA)

            // 绕 Z 轴右手旋转 90°：cos90°=0, sin90°=1
            attribute :>> rotationMatrix {
                attribute :>> elements = (
                    0, -1, 0,
                    1,  0, 0,
                    0,  0, 1 );
            }

            // 沿 X 平移 2（数值，单位由 worldFrame 的 mRefs=metre 解释）
            attribute :>> translationVector {
                attribute :>> elements = (2, 0, 0);
            }
        }
    }

    /* 自动组装出的 4x4 齐次变换矩阵 T：
     *   ( 0, -1, 0, 2,
     *     1,  0, 0, 0,
     *     0,  0, 1, 0,
     *     0,  0, 0, 1 )
     * 语义：x_world = R * x_body + t
     *   body 原点 (0,0,0) -> (2,0,0)
     *   body 的 X 轴 (1,0,0) -> (2,1,0)   （绕 Z 转 90°，body 的 X 指向世界 +Y）
     */

    /* =================================================================
     * Version B —— 等价的“平移/旋转序列”（更贴近人类直觉）
     *
     * 说明：Translation/Rotation 的向量属性类型是 VectorQuantityValue，
     *       而它在库里是抽象的，必须用“具体向量量值类型 + 向量测量参考”实例化。
     *       下面就近定义一个以米为单位的 3D 向量量值类型 Vec3m 供复用。
     * ================================================================= */

    // 3D 向量测量参考（三个轴都是米）
    attribute lengthRef3d : VectorMeasurementReference {
        attribute :>> dimensions = 3;
        attribute :>> mRefs = (metre, metre, metre);
    }

    // 具体的 3D 向量量值类型
    attribute def Vec3m :> ThreeDVectorQuantityValue {
        attribute :>> mRef = lengthRef3d;
        attribute :>> isBound = false;                     // 自由向量，匹配 mRef.isBound=false
    }

    // 序列元素一：绕世界 Z 轴旋转 90°
    attribute rotZ90 : Rotation {
        attribute :>> axisDirection : Vec3m {              // 旋转轴：仅方向有意义
            attribute :>> num = (0, 0, 1);
        }
        attribute :>> angle = 90 [degree];                 // value [unit] 写法
        attribute :>> isIntrinsic = false;                 // 外旋：绕世界(固定)Z 轴
    }

    // 序列元素二：沿世界 X 轴平移 2 米
    attribute trX2 : Translation {
        attribute :>> translationVector : Vec3m {
            attribute :>> num = (2, 0, 0);
        }
    }

    // 用序列描述 bodyFrameB 的放置
    attribute bodyFrameB : MetricFrame3d {
        attribute :>> transformation : TranslationRotationSequence {
            attribute :>> source = worldFrame;
            // 有序序列：先旋转，再平移
            attribute :>> elements = ( rotZ90, trX2 );
        }
    }

    /* 注：A、B 描述的是同一放置（轴绕 Z 转 90°、原点在 (2,0,0)）。
     * 但“序列”形式的精确合成结果还取决于：元素顺序，以及每个 Rotation 的
     * isIntrinsic（内旋/外旋）。如需与 Version A 的 [R|t] 严格逐位相等，
     * 应据库文档(Davenport chained rotations)确认顺序与内/外旋约定。 */

    /* =================================================================
     * Version C —— 恒等变换（位置、朝向都与源坐标系相同）
     * ================================================================= */
    attribute identityFrame : MetricFrame3d {
        attribute :>> transformation : NullTransformation {
            attribute :>> source = worldFrame;
            // rotationMatrix=单位阵、translationVector=0 已由 NullTransformation 给定
        }
    }
}


/*
private import SI::*;
private import ScalarValues::*;
in tx : LengthValue = 50 [mm];

// 将 tx 除以 1 米，系统会自动换算，结果为纯数字 0.05
attribute standard_number : Real = tx / (1 [m]); 

// 如果你想获取以毫米为单位的纯数值，结果就是 50
attribute mm_number : Real = tx / (1 [mm]);
*/


/*

private import SI::*;
private import ISQ::*;
private import ScalarValues::*;

// 修改点：将 := 改为 =，表示这在数学上等同于 1 度（静态绑定）
attribute alpha_Rz : PlaneAngleValue = 1 [degree];  // yaw
attribute a : Real = alpha_Rz / (1 [radian]);

*/


/*
calc def buildTransform_DGSLCS_to_DGSCS {
    in alpha_Rz : PlaneAngleValue;   // yaw
    in beta_Ry  : PlaneAngleValue;   // pitch
    in gamma_Rx : PlaneAngleValue;   // roll
    in tx : LengthValue;  
    in ty : LengthValue;  
    in tz : LengthValue;

    // 必须先除以单位 → 得到 radian 制的纯 Real, 才能喂给 sin/cos
    private attribute a : Real = alpha_Rz / (1 [radian]);
    private attribute b : Real = beta_Ry  / (1 [radian]);
    private attribute g : Real = gamma_Rx / (1 [radian]);

    private attribute ca = cos(a);  
    private attribute sa = sin(a);
    private attribute cb = cos(b);  
    private attribute sb = sin(b);
    private attribute cg = cos(g); 
    private attribute sg = sin(g);

    return result : AffineTransformationMatrix3d {
        // R = Rz(α)·Ry(β)·Rx(γ), row-major, 与你附录公式一致
        attribute :>> rotationMatrix {
            attribute :>> elements = (
                ca*cb,  ca*sb*sg - sa*cg,  ca*sb*cg + sa*sg,
                sa*cb,  sa*sb*sg + ca*cg,  sa*sb*cg - ca*sg,
                 -sb,   cb*sg,             cb*cg );
        }
        // 平移同样要消量纲: 元素是 mRefs(=millimetre) 解释下的纯数值
        attribute :>> translationVector {
            attribute :>> elements = (tx/(1 [mm]), ty/(1 [mm]), tz/(1 [mm]));
        }
    }
}
*/

private import ScalarValues::*;
private import ISQ::*;                    // LengthValue, PlaneAngleValue
private import SI::*;                     // radian, mm
private import TrigFunctions::*;          // sin, cos
private import MeasurementReferences::*;  // AffineTransformationMatrix3d, ThreeDCoordinateFrame

calc def buildTransform_DGSLCS_to_DGSCS {
    // —— 实时 6-DOF 命令 ——
    in alpha_Rz : PlaneAngleValue;        // yaw  (Rz)  ※库若用 AngularMeasureValue, 仅改此 3 处
    in beta_Ry  : PlaneAngleValue;        // pitch(Ry)
    in gamma_Rx : PlaneAngleValue;        // roll (Rx)
    in tx : LengthValue;
    in ty : LengthValue;
    in tz : LengthValue;

    // —— 接线方向(库铁律: [p]_source = M · [p]_target)——
    // 你要 [p]_DGSCS = M · [p]_DGSLCS  ⇒  source = DGSCS, target = DGSLCS
    in dgscs  : ThreeDCoordinateFrame;    // = source
    in dgslcs : ThreeDCoordinateFrame;    // = target

    // —— 角度消量纲 → radian 制纯 Real, 才能喂 sin/cos ——
    private attribute a : Real = alpha_Rz / (1 [radian]);
    private attribute b : Real = beta_Ry  / (1 [radian]);
    private attribute g : Real = gamma_Rx / (1 [radian]);

    private attribute ca : Real = cos(a);   private attribute sa : Real = sin(a);
    private attribute cb : Real = cos(b);   private attribute sb : Real = sin(b);
    private attribute cg : Real = cos(g);   private attribute sg : Real = sin(g);

    return result : AffineTransformationMatrix3d {
        attribute :>> source = dgscs;     // 方向写死, 接不反
        attribute :>> target = dgslcs;

        // R = Rz(α)·Ry(β)·Rx(γ), row-major, 正确符号(非上传图)
        attribute :>> rotationMatrix {
            attribute :>> elements = (
                ca*cb,  ca*sb*sg - sa*cg,  ca*sb*cg + sa*sg,
                sa*cb,  sa*sb*sg + ca*cg,  sa*sb*cg - ca*sg,    // ← (2,2),(2,3) 正确
                 -sb,   cb*sg,             cb*cg );
        }
        // 平移 = DGSLCS 原点在 DGSCS(=source) 坐标下的位置; 除数单位须 == DGSCS 的 mRefs
        attribute :>> translationVector {
            attribute :>> elements = ( tx/(1 [mm]), ty/(1 [mm]), tz/(1 [mm]) );  // mm ⇔ DGSCS mRefs
        }
    }
}



/*

private import TrigFunctions::*;     // rad, sin, cos, pi
private import Geometry::*;          // AffineTransformationMatrix3d 等(包名以库源为准)

calc def buildTransform_DGSLCS_to_DGSCS {
    in alpha_Rz_deg : Real;   // yaw   (度)
    in beta_Ry_deg  : Real;   // pitch (度)
    in gamma_Rx_deg : Real;   // roll  (度)
    in tx : Real;             // 平移数值, 单位由目标 frame 的 mRefs 解释
    in ty : Real;
    in tz : Real;

    // 度 → 弧度: 用标准库, 不再手写 *pi/180, 也不再有 PlaneAngleValue/Real 的类型缺口
    private attribute a : Real = rad(alpha_Rz_deg);
    private attribute b : Real = rad(beta_Ry_deg);
    private attribute g : Real = rad(gamma_Rx_deg);

    private attribute ca : Real = cos(a);  private attribute sa : Real = sin(a);
    private attribute cb : Real = cos(b);  private attribute sb : Real = sin(b);
    private attribute cg : Real = cos(g);  private attribute sg : Real = sin(g);

    return result : AffineTransformationMatrix3d {
        attribute :>> rotationMatrix {        // R = Rz(a)·Ry(b)·Rx(g), row-major
            attribute :>> elements = (
                ca*cb,  ca*sb*sg - sa*cg,  ca*sb*cg + sa*sg,
                sa*cb,  sa*sb*sg + ca*cg,  sa*sb*cg - ca*sg,
                 -sb,   cb*sg,             cb*cg );
        }
        attribute :>> translationVector {
            attribute :>> elements = (tx, ty, tz);
        }
    }
}


*/
