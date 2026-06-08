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