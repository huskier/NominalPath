
private import SI::*;
private import Spatial::CoordinateSystems::*;
private import Spatial::SpatialItems::*;

// 1. 你的纯数学矩阵生成器 (保持你的代码不变)
calc def buildTransform_DGSLCS_to_DGSCS { ... }

// 2. 在具体的系统架构中应用
part def LPP_EUV_Chamber {
    
    // 定义液滴发生器坐标系 (DGSLCS)
    // 在这里，底层标准库的 CartesianCoordinateSystem 会隐式或显式地包含 mRefs 的定义（例如基准单位设为 mm）
    part dgslcs_frame : CartesianCoordinateSystem;
    
    // 定义系统全局坐标系 (DGSCS)
    part dgscs_frame : CartesianCoordinateSystem;

    // 3. 调用你的 calc def 获取纯数学矩阵
    private attribute pure_math_matrix = buildTransform_DGSLCS_to_DGSCS(
        alpha_Rz = 1 [degree], 
        beta_Ry = 0 [degree], 
        gamma_Rx = 0 [degree], 
        tx = 50 [mm], ty = 0 [mm], tz = -10 [mm]
    );

    // 4. 将纯数学矩阵绑定到物理坐标变换上
    // Spatial 库中的 CoordinateTransformation 会将 pure_math_matrix 
    // 解释为连接 dgslcs_frame 和 dgscs_frame 的桥梁，
    // 此时它才会真正与坐标系的 mRefs (如毫米尺度) 发生物理层面的结合。
    part transform_dgslcs_to_dgscs : CoordinateTransformation {
        attribute :>> transformationMatrix = pure_math_matrix;
        ref sourceCoordinateSystem = dgslcs_frame;
        ref targetCoordinateSystem = dgscs_frame;
    }
}


// 1. 定义一个通用的三维点仿射变换函数
calc def applyTransformToPoint {
    in matrix : AffineTransformationMatrix3d;
    in p_source : Real[3]; // 输入在源坐标系下的纯数字坐标 (例如: mm)

    private attribute R = matrix.rotationMatrix.elements;
    private attribute T = matrix.translationVector.elements;

    // 显式展开矩阵乘法和向量加法 
    // (注意: SysML v2 数组索引通常从 1 开始)
    return p_target : Real[3] = (
        R[1]*p_source[1] + R[2]*p_source[2] + R[3]*p_source[3] + T[1],
        R[4]*p_source[1] + R[5]*p_source[2] + R[6]*p_source[3] + T[2],
        R[7]*p_source[1] + R[8]*p_source[2] + R[9]*p_source[3] + T[3]
    ); 
}

// 2. 在系统上下文中调用它
part def LPP_EUV_Chamber_Simulation {
    
    // 获取我们之前生成的纯数学转换矩阵
    private attribute pure_math_matrix = buildTransform_DGSLCS_to_DGSCS(
        alpha_Rz = 1 [degree], 
        beta_Ry = 0 [degree], 
        gamma_Rx = 0 [degree], 
        tx = 50 [mm], ty = 0 [mm], tz = -10 [mm]
    );

    // 假设：我们要追踪的液滴在发生器坐标系 (DGSLCS) 中的初始剥离单位位置
    attribute droplet_pos_dgslcs : Real[3] = (0.0, 0.0, -5.0); 

    // 使用函数进行坐标变换求值
    attribute droplet_pos_dgscs : Real[3] = applyTransformToPoint(
        matrix = pure_math_matrix,
        p_source = droplet_pos_dgslcs
    );
}


import Spatial::CoordinateSystems::*;
import Spatial::SpatialItems::*;

part def ChamberArchitecture {
    part dgslcs_frame : CartesianCoordinateSystem;
    part dgscs_frame : CartesianCoordinateSystem;

    // 定义矩阵和坐标系桥梁
    private attribute math_matrix = buildTransform_DGSLCS_to_DGSCS(...);
    part transform_bridge : CoordinateTransformation {
        attribute :>> transformationMatrix = math_matrix;
        ref sourceCoordinateSystem = dgslcs_frame;
        ref targetCoordinateSystem = dgscs_frame;
    }

    // 声明一个真实的物理对象（液滴）
    part droplet : SpatialItem {
        // 声明它在 DGSLCS 下的位置
        attribute pos_in_local : PositionVector {
            ref coordinateSystem = dgslcs_frame;
            attribute :>> coordinates = (0.0, 0.0, -5.0);
        }
        
        // 声明它在 DGSCS 下的位置（只声明，不手动计算）
        attribute pos_in_global : PositionVector {
            ref coordinateSystem = dgscs_frame;
        }
        
        // 系统底层约束：这两个位置描述的是同一个物理对象的同一个物理位置
        // 当模型对接给高级求解器时，求解器会通过 transform_bridge 自动算出 pos_in_global 的值
    }
}


//范式一：参数化绑定（Parametric Binding）—— 最实用、最适合求值
part droplet : SpatialItem {
    
    // 1. 给出源数据（在发生器坐标系下）
    attribute pos_in_local : PositionVector {
        ref coordinateSystem = dgslcs_frame;
        attribute :>> coordinates = (0.0, 0.0, -5.0);
    }
    
    // 2. 施加底层约束：目标坐标系的数值，永远等于源数值经过矩阵变换后的结果
    attribute pos_in_global : PositionVector {
        ref coordinateSystem = dgscs_frame;
        
        // 【核心约束在此】直接使用等号（=）绑定数学表达式
        attribute :>> coordinates = applyTransformToPoint(
            matrix = transform_bridge.transformationMatrix, 
            p_source = pos_in_local.coordinates
        );
    }
}

//范式二：声明式约束断言（Constraint Assertion）—— 最符合 MBSE 精神
// 1. 定义一条“铁律”（物理法则）：
// 规定任何通过矩阵连接的两个点，必须满足公式 P_target == M * P_source
constraint def PointTransformationConstraint {
    in matrix : AffineTransformationMatrix3d;
    in source_pos : Real[3];
    in target_pos : Real[3];
    
    // 核心约束表达式：使用双等号 (==) 表示必须相等的逻辑关系
    target_pos == applyTransformToPoint(matrix, source_pos)
}

part def ChamberArchitecture {
    // ... 前面的坐标系和矩阵定义 ...

    part droplet : SpatialItem {
        // 假设这两个值是从不同地方获取的
        attribute pos_in_local : PositionVector { ... }
        attribute pos_in_global : PositionVector { ... }
        
        // 2. 施加底层约束：把这条铁律像“紧箍咒”一样套在这三个属性上
        assert constraint physical_equivalence : PointTransformationConstraint {
            in matrix = transform_bridge.transformationMatrix;
            in source_pos = pos_in_local.coordinates;
            in target_pos = pos_in_global.coordinates;
        }
    }
}

//总结建议
//如果你的目的是正向设计（已知源位置，想在模型里推导并使用全局位置），用范式一（参数化绑定）。
//如果你的目的是架构校验（汇集了各方数据，想检查它们在物理上是否矛盾），用范式二（约束断言）。

