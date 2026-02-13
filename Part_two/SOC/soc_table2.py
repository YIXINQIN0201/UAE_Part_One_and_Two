# import rasterio
# import pandas as pd
# import numpy as np
# import os
#
# # ================= 配置部分 =================
# # 文件夹路径
# folder = r"D:\deep\fromFeishu\Land_Cover_Each_Year"
# file_2000 = os.path.join(folder, "C3S_2018_reclass_300m.tif")
# file_2015 = os.path.join(folder, "C3S_2024_reclass_300m.tif")
#
# # 类别对应表 (根据你提供的)
# class_mapping = {
#     1: "Forest",  # Forest Land
#     2: "Grassland",
#     3: "Cropland",
#     4: "Wetland",
#     5: "Artificial",
#     6: "Other",  # Other Land
#     7: "Water"  # Water bodies
# }
#
# # 【关键】你指定的关注清单 (只计算这些!)
# # 格式: (2000年代码, 2015年代码)
# # Cropland(3), Artificial(5), Other(6), Grassland(2), Water(7), Wetland(4), Forest(1)
# #  baseline
# # target_transitions = [
# #     (3, 5),  # Cropland -> Artificial
# #     (3, 6),  # Cropland -> Other
# #     (2, 6),  # Grassland -> Other
# #     (6, 5),  # Other -> Artificial
# #     (6, 3),  # Other -> Cropland
# #     (6, 2),  # Other -> Grassland
# #     (7, 5),  # Water -> Artificial
# #     (7, 1),  # Water -> Forest
# #     (7, 2),  # Water -> Grassland
# #     (7, 6),  # Water -> Other
# #     (7, 4)  # Water -> Wetland
# # ]
#
# # reporting period
# target_transitions = [
#     (5, 3),  # Artificial -> Cropland
#     (5, 2),  # Artificial -> Grassland
#     (5, 6),  # Artificial -> Other
#     (3, 5),  # Cropland -> Artificial
#     (3, 2),  # Cropland -> Grassland
#     (3, 6),  # Cropland -> Other
#     (2, 5),  # Grassland -> Artificial
#     (2, 3),  # Grassland -> Cropland
#     (2, 1),  # Grassland -> Forest
#     (2, 6),  # Grassland -> Other
#     (6, 5),  # Other -> Artificial
#     (6, 3),  # Other -> Cropland
# (6, 2),  # Other -> Grassland
# (6, 7),  # Other -> Water
# (7, 5),  # Water -> Artificial
# (7, 2),  # Water -> Grassland
# (7, 6),  # Water -> Other
# (7, 4),  # Water -> Wetland
# (4, 7),  # Wetland -> Water
#
#
# ]
#
# # ================= 计算过程 =================
# print("正在读取数据...")
# with rasterio.open(file_2000) as src0, rasterio.open(file_2015) as src1:
#     d2018 = src0.read(1).flatten()
#     d2024 = src1.read(1).flatten()
#
#     # 自动获取像素分辨率
#     res_x = src0.res[0]  # 通常是 300
#     res_y = src0.res[1]  # 通常是 300
#
#     # 判断是否为经纬度 (如果是经纬度，res_x会很小，如0.0027)
#     if src0.crs.is_geographic:
#         # C3S 默认通常是经纬度 WGS84
#         # 在阿联酋纬度，300m大约对应 0.09 km2 (简单估算)
#         # 或者使用更精确的平均纬度校正
#         pixel_area_km2 = (300 * 300) / 1e6  # 强制假设它是300m数据
#         print("提示: 数据是经纬度坐标，按标称300m分辨率计算面积 (0.09 km²/像素)")
#     else:
#         pixel_area_km2 = (res_x * res_y) / 1e6
#         print(f"提示: 投影坐标，像素面积 {pixel_area_km2:.4f} km²")
#
# print("正在筛选特定转换类型...")
#
# results = []
#
# for start_code, end_code in target_transitions:
#     # 1. 找到符合条件的像素索引
#     # 逻辑: 2000年等于start 且 2015年等于end
#     count = np.sum((d2018 == start_code) & (d2024 == end_code))
#
#     # 2. 计算面积
#     area_km2 = count * pixel_area_km2
#
#     # 3.以此构建记录
#     name_from = class_mapping.get(start_code, str(start_code))
#     name_to = class_mapping.get(end_code, str(end_code))
#
#     results.append({
#         "Land Conversion": f"{name_from} -> {name_to}",
#         "Net area change (km²)": round(area_km2, 4),  # 保留4位小数
#         "Pixel Count": count
#     })
#
# # ================= 输出结果 =================
# df = pd.DataFrame(results)
# print("-" * 30)
# print(df)
# print("-" * 30)
#
# # 保存
# df.to_csv("D:\deep\SoilGrids_Data\\appendix_table2\\reporting_Step1_Area_Change.csv", index=False)
# print("第一列计算完成！结果已保存为 reporting_Step1_Area_Change.csv")

######################################################################
# 算Initial SOC Stock (t/ha)和Final SOC stock (t/ha)
import rasterio
import numpy as np
import pandas as pd
from rasterio.warp import reproject, Resampling
import os

# ================= 1. 文件与路径配置 =================
folder_lc = r"D:\deep\fromFeishu\Land_Cover_Each_Year"
file_lc_2000 = os.path.join(folder_lc, "C3S_2018_reclass_300m.tif")
file_lc_2015 = os.path.join(folder_lc, "C3S_2024_reclass_300m.tif")
file_ocs = "D:\deep\SoilGrids_Data\\appendix_table2\\UAE_SoilGrids_OCS_0-30cm_mean.tif"  # 假设在当前目录

# 结果保存名
output_csv = "D:\deep\SoilGrids_Data\\appendix_table2\\reporting_Step2_SOC_Stocks_per_ha.csv"

# ================= 2. 参数设置 =================
# IPCC Factors (Flu) - Tropical Dry Climate
factors = {
    1: 1.00,  # Forest
    2: 1.00,  # Grassland
    3: 0.58,  # Cropland
    4: 1.00,  # Wetland
    5: 0.50,  # Artificial
    6: 0.80,  # Other Land
    7: 0.00  # Water
}

# 类别名称映射
class_names = {
    1: "Forest", 2: "Grassland", 3: "Cropland",
    4: "Wetland", 5: "Artificial", 6: "Other", 7: "Water"
}

# # 【目标转换清单】(与您要求的一致)
# target_transitions = [
#     (3, 5),  # Cropland -> Artificial
#     (3, 6),  # Cropland -> other
#     (2, 6),  # grassland -> other
#     (6, 5),  # other -> artificial
#     (6, 3),  # other -> cropland
#     (6, 2),  # other -> grassland
#     (7, 5),  # water -> artificial
#     (7, 1),  # water -> forest
#     (7, 2),  # water -> grassland
#     (7, 6),  # water -> other
#     (7, 4)  # water -> wetland
# ]

# reporting period
target_transitions = [
    (5, 3),  # Artificial -> Cropland
    (5, 2),  # Artificial -> Grassland
    (5, 6),  # Artificial -> Other
    (3, 5),  # Cropland -> Artificial
    (3, 2),  # Cropland -> Grassland
    (3, 6),  # Cropland -> Other
    (2, 5),  # Grassland -> Artificial
    (2, 3),  # Grassland -> Cropland
    (2, 1),  # Grassland -> Forest
    (2, 6),  # Grassland -> Other
    (6, 5),  # Other -> Artificial
    (6, 3),  # Other -> Cropland
(6, 2),  # Other -> Grassland
(6, 7),  # Other -> Water
(7, 5),  # Water -> Artificial
(7, 2),  # Water -> Grassland
(7, 6),  # Water -> Other
(7, 4),  # Water -> Wetland
(4, 7),  # Wetland -> Water
]

# ================= 3. 核心处理函数 =================

def calculate_step2():
    print("1. 读取 Land Cover 数据...")
    with rasterio.open(file_lc_2000) as src0:
        d2000 = src0.read(1)
        profile = src0.profile
        # 获取用于重采样的目标参数
        dst_transform = src0.transform
        dst_crs = src0.crs
        dst_shape = d2000.shape

    with rasterio.open(file_lc_2015) as src1:
        d2015 = src1.read(1)

    print("2. 读取并对齐 SoilGrids OCS 数据...")
    # 初始化一个与 LandCover 一模一样的空数组
    ocs_aligned = np.zeros(dst_shape, dtype=np.float32)

    with rasterio.open(file_ocs) as src_ocs:
        reproject(
            source=rasterio.band(src_ocs, 1),
            destination=ocs_aligned,
            src_transform=src_ocs.transform,
            src_crs=src_ocs.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear
        )
    print("   对齐完成。")

    print("3. 开始计算各类型的碳密度 (t/ha)...")
    results = []

    for start_code, end_code in target_transitions:
        # 1. 生成掩膜：找到所有属于该转换类型的像素
        mask = (d2000 == start_code) & (d2015 == end_code)

        # 2. 提取这些像素的 SOC_Ref 值
        soc_ref_values = ocs_aligned[mask]

        # 如果没有找到任何像素（该转换在地图上不存在），填 0 或 NaN
        if soc_ref_values.size == 0:
            print(f"   警告: 未在地图上发现转换 {start_code}->{end_code}")
            avg_ref = 0
        else:
            # 计算该区域 SOC_Ref 的平均值
            avg_ref = np.mean(soc_ref_values)

        # 3. 计算 Initial 和 Final Stock (t/ha)
        # 公式: Stock = Ref_Mean * Flu_Factor
        # 注意：这里我们假设 Fmg 和 Fi 都是 1

        # 初始因子
        f_start = factors.get(start_code, 1.0)
        # 结束因子
        f_end = factors.get(end_code, 1.0)

        initial_stock_tha = avg_ref * f_start
        final_stock_tha = avg_ref * f_end

        # 构建结果行
        name_from = class_names.get(start_code, str(start_code))
        name_to = class_names.get(end_code, str(end_code))

        results.append({
            "Land Conversion": f"{name_from} -> {name_to}",
            "Initial SOC Stock (t/ha)": round(initial_stock_tha, 2),
            "Final SOC stock (t/ha)": round(final_stock_tha, 2),
            "Avg SOC_Ref (Check)": round(avg_ref, 2)  # 仅供检查用
        })

    # ================= 4. 输出 =================
    df = pd.DataFrame(results)
    print("-" * 30)
    print(df[["Land Conversion", "Initial SOC Stock (t/ha)", "Final SOC stock (t/ha)"]])
    print("-" * 30)

    df.to_csv(output_csv, index=False)
    print(f"第二步计算完成！结果已保存至: {output_csv}")


if __name__ == "__main__":
    calculate_step2()