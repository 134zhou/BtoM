import pandas as pd
import numpy as np

def csv_to_structured_vtk(csv_path, vtk_path):
    # 1. 读取数据并排序 (确保顺序: Z -> Y -> X)
    df = pd.read_csv(csv_path)
    # 必须严格排序以匹配 VTK 遍历规则
    df = df.sort_values(by=['z', 'y', 'x'])
    
    for col in ['Mx', 'My', 'Mz', 'x', 'y', 'z']:
        df[col] = pd.to_numeric(df[col], errors='coerce') 
    df = df.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # 2. 获取网格维度
    x_coords = sorted(df['x'].unique())
    y_coords = sorted(df['y'].unique())
    z_coords = sorted(df['z'].unique())
    
    nx, ny, nz = len(x_coords), len(y_coords), len(z_coords)
    
    # 获取间距（假设是均匀网格）
    dx = x_coords[1] - x_coords[0] if nx > 1 else 1.0
    dy = y_coords[1] - y_coords[0] if ny > 1 else 1.0
    dz = z_coords[1] - z_coords[0] if nz > 1 else 1.0

    with open(vtk_path, 'w') as f:
        # 写入头部
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Magnetic Field Scan Structured\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        f.write(f"ORIGIN {x_coords[0]} {y_coords[0]} {z_coords[0]}\n")
        f.write(f"SPACING {dx} {dy} {dz}\n")
        
        f.write(f"POINT_DATA {len(df)}\n")
        
        # 写入模值（标量：磁场强度）
        f.write("SCALARS B_Magnitude double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for _, r in df.iterrows():
            # 计算磁化强度 M = sqrt(Mx^2 + My^2 + Mz^2)
            mag = np.sqrt(r['Mx']**2 + r['My']**2 + r['Mz']**2)
            f.write(f"{mag:.6f}\n")
            
        # 写入矢量（磁化方向）
        f.write("\nVECTORS M_Vector double\n")
        for _, r in df.iterrows():
            f.write(f"{r['Mx']} {r['My']} {r['Mz']}\n")

    print(f"✅ 转换完成！网格大小: {nx}x{ny}x{nz}")

# 调用示例
csv_to_structured_vtk('./inverted_M_results_test.csv', './M_structured_test.vtk')
