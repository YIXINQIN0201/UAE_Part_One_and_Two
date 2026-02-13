# import rasterio
# from rasterio.mask import mask
# import geopandas as gpd
# import numpy as np
# import pandas as pd
# import os
# import matplotlib.pyplot as plt
# from matplotlib.colors import ListedColormap
#
# # ================= 配置路径 =================
# # 1. 输入文件路径
# path_lc = r"D:\deep\fromFeishu\hmq\Land_Cover_Degradation\status_2018_2024_deg_stable_imp.tif"
# path_prod = r"D:\deep\fromFeishu\qyx\land_productivity\combined_degradation_2024_table45.tif"
# path_soc = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2018_2024.tif"
#
# # 2. 辅助文件
# path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
#
# # 3. 输出文件路径
# output_dir = r"D:\deep\SoilGrids_Data\appendix_table4"
# output_tif = os.path.join(output_dir, "SDG_15_3_1_Indicator_2024.tif")
# output_csv = os.path.join(output_dir, "National_Estimates_SDG_15_3_1.csv")
#
#
# # ================= 辅助函数 =================
# def calculate_area_weighted(data, transform, mask_array, crs):
#     """计算加权面积 (km²)"""
#     res_x = transform[0]
#     res_y = -transform[4]
#     height, width = data.shape
#
#     # 简单的 WGS84 判断
#     is_wgs84 = (crs.to_epsg() == 4326) or (crs.is_geographic)
#
#     if is_wgs84:
#         R = 6371.0
#         y_origin = transform[5]
#         pixel_height = transform[4]
#         y_indices = np.arange(height) + 0.5
#         latitudes = y_origin + y_indices * pixel_height
#         lat_rad = np.radians(latitudes)
#         pixel_width_km = np.abs(res_x) * (np.pi / 180) * R * np.cos(lat_rad)
#         pixel_height_km = np.abs(res_y) * (np.pi / 180) * R
#         row_pixel_area = pixel_width_km * pixel_height_km
#         count_per_row = np.sum(mask_array, axis=1)
#         return np.sum(count_per_row * row_pixel_area)
#     else:
#         pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
#         return np.sum(mask_array) * pixel_area_km2
#
#
# def ensure_alignment(raster_path, shp_gdf, target_shape=None):
#     """读取栅格，强制与 Shapefile 坐标系对齐，并裁剪"""
#     with rasterio.open(raster_path) as src:
#         # 1. 关键修复：检查坐标系是否一致
#         if str(src.crs) != str(shp_gdf.crs):
#             print(f"   [提示] 坐标系不一致，正在转换 Shapefile 匹配 {os.path.basename(raster_path)}...")
#             print(f"   Raster CRS: {src.crs}, Shapefile CRS: {shp_gdf.crs}")
#             shp_gdf = shp_gdf.to_crs(src.crs)
#
#         # 2. 执行裁剪
#         try:
#             out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
#             data = out_image[0]
#         except ValueError as e:
#             print(f"   [错误] 裁剪失败: {e}")
#             print(f"   Raster Bounds: {src.bounds}")
#             print(f"   Shapefile Bounds: {shp_gdf.total_bounds}")
#             raise e
#
#         # 3. 如果需要强制形状对齐 (处理1-2个像素的误差)
#         if target_shape is not None and data.shape != target_shape:
#             # 简单处理：取左上角对齐的最小公共区域
#             h = min(data.shape[0], target_shape[0])
#             w = min(data.shape[1], target_shape[1])
#             data = data[:h, :w]
#
#         return data, out_transform, src.crs, shp_gdf
#
#
# def main():
#     print("--- 开始 SDG 15.3.1 综合计算 (Fix版) ---")
#
#     # 1. 读取 Shapefile
#     print("1. 读取 Shapefile...")
#     shp = gpd.read_file(path_shp)
#
#     # 2. 读取基准数据 (Land Cover) 并更新 Shapefile 坐标系
#     print(f"2. 读取 Land Cover: {os.path.basename(path_lc)}")
#     arr_lc, out_transform, crs, shp_aligned = ensure_alignment(path_lc, shp)
#     base_shape = arr_lc.shape
#     print(f"   基准形状: {base_shape}")
#
#     # 3. 读取其他数据 (使用已对齐的 shapefile)
#     print(f"3. 读取 Productivity: {os.path.basename(path_prod)}")
#     # 注意：这里我们传入 base_shape 以确保尺寸严格一致
#     arr_prod, _, _, _ = ensure_alignment(path_prod, shp_aligned, target_shape=base_shape)
#
#     print(f"4. 读取 SOC: {os.path.basename(path_soc)}")
#     arr_soc, _, _, _ = ensure_alignment(path_soc, shp_aligned, target_shape=base_shape)
#
#     # 再次检查尺寸 (因为 ensure_alignment 可能会剪切 arr_prod/soc，也可能需要反过来剪切 arr_lc)
#     # 最稳妥的方法是取三者最小的共同尺寸
#     min_h = min(arr_lc.shape[0], arr_prod.shape[0], arr_soc.shape[0])
#     min_w = min(arr_lc.shape[1], arr_prod.shape[1], arr_soc.shape[1])
#
#     arr_lc = arr_lc[:min_h, :min_w]
#     arr_prod = arr_prod[:min_h, :min_w]
#     arr_soc = arr_soc[:min_h, :min_w]
#
#     print(f"   最终计算尺寸: {arr_lc.shape}")
#
#     # 4. 数据清洗 (处理 NoData)
#     valid_mask_lc = np.isin(arr_lc, [1, 2, 3])
#     valid_mask_prod = np.isin(arr_prod, [1, 2, 3])
#     valid_mask_soc = np.isin(arr_soc, [1, 2, 3])
#
#     # 只要任意图层有数据，就纳入计算
#     analyzed_mask = valid_mask_lc | valid_mask_prod | valid_mask_soc
#
#     print("5. 执行 'One Out, All Out' 逻辑...")
#     final_sdg = np.zeros(arr_lc.shape, dtype=np.uint8)
#
#     # 逻辑 1: Degraded (1) - 任意一个为1
#     mask_degraded = (
#             ((arr_lc == 1) & valid_mask_lc) |
#             ((arr_prod == 1) & valid_mask_prod) |
#             ((arr_soc == 1) & valid_mask_soc)
#     )
#
#     # 逻辑 2: Improved (3) - 没有1，且至少有一个3
#     mask_improved_raw = (
#             ((arr_lc == 3) & valid_mask_lc) |
#             ((arr_prod == 3) & valid_mask_prod) |
#             ((arr_soc == 3) & valid_mask_soc)
#     )
#     mask_improved = mask_improved_raw & (~mask_degraded) & analyzed_mask
#
#     # 逻辑 3: Stable (2) - 剩下的
#     mask_stable = analyzed_mask & (~mask_degraded) & (~mask_improved)
#
#     final_sdg[mask_stable] = 2
#     final_sdg[mask_improved] = 3
#     final_sdg[mask_degraded] = 1
#
#     # 5. 计算面积
#     print("6. 计算国家估算值...")
#     area_degraded = calculate_area_weighted(final_sdg, out_transform, final_sdg == 1, crs)
#     area_stable = calculate_area_weighted(final_sdg, out_transform, final_sdg == 2, crs)
#     area_improved = calculate_area_weighted(final_sdg, out_transform, final_sdg == 3, crs)
#
#     total_analyzed_area = area_degraded + area_stable + area_improved
#
#     # 计算比例
#     prop_degraded = (area_degraded / total_analyzed_area * 100) if total_analyzed_area > 0 else 0
#     prop_stable = (area_stable / total_analyzed_area * 100) if total_analyzed_area > 0 else 0
#     prop_improved = (area_improved / total_analyzed_area * 100) if total_analyzed_area > 0 else 0
#
#     # 输出结果
#     df_result = pd.DataFrame({
#         "Category": ["Degraded", "Stable", "Improved", "Total Analyzed"],
#         "Area (km2)": [area_degraded, area_stable, area_improved, total_analyzed_area],
#         "Proportion (%)": [prop_degraded, prop_stable, prop_improved, 100.0]
#     })
#
#     print("\n====== National Estimates (SDG 15.3.1) ======")
#     print(df_result.round(2))
#
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#     df_result.to_csv(output_csv, index=False)
#
#     # 6. 保存最终 TIF
#     meta = {
#         'driver': 'GTiff',
#         'height': final_sdg.shape[0],
#         'width': final_sdg.shape[1],
#         'transform': out_transform,
#         'crs': crs,
#         'count': 1,
#         'dtype': 'uint8',
#         'nodata': 0,
#         'compress': 'lzw'
#     }
#     with rasterio.open(output_tif, "w", **meta) as dst:
#         dst.write(final_sdg, 1)
#     print(f"最终结果已保存: {output_tif}")
#
#     # 7. 生成预览图
#     print("生成预览图...")
#     plt.figure(figsize=(10, 8))
#     cmap = ListedColormap(['white', 'red', '#ffcc00', 'green'])
#     plt.imshow(final_sdg, cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
#     plt.title(f"SDG 15.3.1 Result\nDegraded: {prop_degraded:.2f}%")
#     plt.axis('off')
#     plt.show()
#
#
# if __name__ == "__main__":
#     main()



#==========================================================
#Correct
#==========================================================
import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# ================= 1. 配置路径 (已更新为你提供的新路径) =================

# reporting period
# path_lc = r"D:\deep\fromFeishu\hmq\Land_Cover_Degradation\status_2018_2024_deg_stable_imp.tif"
# path_prod = r"D:\deep\fromFeishu\qyx\land_productivity_20260122\combined_degradation_2024_table45.tif"
# path_soc = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2018_2024.tif"

# baseline period
path_lc = r"D:\deep\fromFeishu\hmq\status_2000_2015_deg_stable_imp.tif"
path_prod = r"D:\deep\fromFeishu\qyx\land_productivity_20260122\combined_degradation_2015_table45.tif"
path_soc = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2000_2015.tif"

#mask
path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"

#  reporting period
# output_dir = r"D:\deep\SoilGrids_Data\appendix_table4\reporting"
# output_tif = os.path.join(output_dir, "SDG_15_3_1_Indicator_Final_Corrected.tif")
# output_csv = os.path.join(output_dir, "National_Estimates_SDG_15_3_1_Corrected.csv")

# baseline period
output_dir = r"D:\deep\SoilGrids_Data\appendix_table4\baseline"
output_tif = os.path.join(output_dir, "SDG_15_3_1_Indicator_Final.tif")
output_csv = os.path.join(output_dir, "National_Estimates_SDG_15_3_1.csv")


# ================= 2. 辅助函数 =================

def remap_array(data, mapping_dict, name="Unknown"):
    """
    将数组中的值根据字典进行映射
    例如: {0: 2, 1: 3, 2: 1} 表示把0变成2，1变成3...
    """
    print(f"   正在重映射 {name} 的像素值...")
    # 创建一个由 0 填充的同样大小的数组
    new_data = np.zeros_like(data)

    # 遍历映射字典
    for old_val, new_val in mapping_dict.items():
        # 找到旧值的位置
        mask_loc = (data == old_val)
        new_data[mask_loc] = new_val
        count = np.sum(mask_loc)
        print(f"     - 原值 {old_val} -> 新值 {new_val} (处理了 {count} 个像素)")

    return new_data


def ensure_alignment(raster_path, shp_gdf, target_shape=None):
    """读取栅格，强制对齐，并返回数据"""
    with rasterio.open(raster_path) as src:
        # 坐标系检查
        if str(src.crs) != str(shp_gdf.crs):
            shp_gdf = shp_gdf.to_crs(src.crs)

        # 裁剪
        out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
        data = out_image[0]

        # 形状对齐
        if target_shape is not None and data.shape != target_shape:
            h = min(data.shape[0], target_shape[0])
            w = min(data.shape[1], target_shape[1])
            data = data[:h, :w]

        return data, out_transform, src.crs, shp_gdf


def calculate_area_weighted(data, transform, mask_array, crs):
    """计算加权面积 (km²)"""
    res_x = transform[0]
    res_y = -transform[4]
    height, width = data.shape
    is_wgs84 = (crs.to_epsg() == 4326) or (crs.is_geographic)

    if is_wgs84:
        R = 6371.0
        y_origin = transform[5]
        pixel_height = transform[4]
        y_indices = np.arange(height) + 0.5
        latitudes = y_origin + y_indices * pixel_height
        lat_rad = np.radians(latitudes)
        pixel_width_km = np.abs(res_x) * (np.pi / 180) * R * np.cos(lat_rad)
        pixel_height_km = np.abs(res_y) * (np.pi / 180) * R
        row_pixel_area = pixel_width_km * pixel_height_km
        count_per_row = np.sum(mask_array, axis=1)
        return np.sum(count_per_row * row_pixel_area)
    else:
        pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
        return np.sum(mask_array) * pixel_area_km2


# ================= 3. 主流程 =================
def main():
    print("--- 开始 SDG 15.3.1 最终计算 (含数值重映射修正) ---")

    # 1. 读取 Shapefile
    try:
        shp = gpd.read_file(path_shp)
    except Exception as e:
        print(f"Shapefile 读取失败: {e}")
        return

    # 2. 读取并对齐数据
    print("1. 读取并对齐 Land Cover (基准)...")
    arr_lc_raw, out_transform, crs, shp_aligned = ensure_alignment(path_lc, shp)
    base_shape = arr_lc_raw.shape

    print("2. 读取 Productivity 和 SOC...")
    arr_prod_raw, _, _, _ = ensure_alignment(path_prod, shp_aligned, target_shape=base_shape)
    arr_soc_raw, _, _, _ = ensure_alignment(path_soc, shp_aligned, target_shape=base_shape)

    # 统一裁剪到最小尺寸
    min_h = min(arr_lc_raw.shape[0], arr_prod_raw.shape[0], arr_soc_raw.shape[0])
    min_w = min(arr_lc_raw.shape[1], arr_prod_raw.shape[1], arr_soc_raw.shape[1])

    arr_lc_raw = arr_lc_raw[:min_h, :min_w]
    arr_prod_raw = arr_prod_raw[:min_h, :min_w]
    arr_soc_raw = arr_soc_raw[:min_h, :min_w]

    # ================= 关键修正：重映射数值 =================
    print("\n3. 执行数值标准化 (统一为: 1=Degraded, 2=Stable, 3=Improved)...")

    # 映射规则：原值(键) -> 标准值(值)
    # 原始 Land Cover/Productivity: 0=Stable, 1=Improved, 2=Degraded
    # 目标 SDG 标准: 2=Stable, 3=Improved, 1=Degraded
    mapping_rule = {
        0: 2,  # Stable
        1: 3,  # Improved
        2: 1  # Degraded
    }

    arr_lc = remap_array(arr_lc_raw, mapping_rule, name="Land Cover")
    arr_prod = remap_array(arr_prod_raw, mapping_rule, name="Productivity")

    # SOC 已经是标准的 (1,2,3)，不需要重映射，但为了保险，过滤掉杂值
    print("   SOC 数据无需重映射 (假设已为标准格式)...")
    arr_soc = arr_soc_raw.copy()
    # 确保只保留 0,1,2,3
    arr_soc[~np.isin(arr_soc, [0, 1, 2, 3])] = 0

    # ================= SDG 逻辑计算 =================
    print("\n4. 执行 '一票否决' (One Out, All Out) 逻辑...")

    # 有效性掩膜 (只要有任意非0值)
    valid_mask_lc = np.isin(arr_lc, [1, 2, 3])
    valid_mask_prod = np.isin(arr_prod, [1, 2, 3])
    valid_mask_soc = np.isin(arr_soc, [1, 2, 3])

    analyzed_mask = valid_mask_lc | valid_mask_prod | valid_mask_soc

    final_sdg = np.zeros(base_shape[:2], dtype=np.uint8)  # 确保维度正确
    # 截取到min尺寸
    final_sdg = final_sdg[:min_h, :min_w]

    # 1. Degraded (1): 任意一个为 1
    mask_degraded = (
            ((arr_lc == 1) & valid_mask_lc) |
            ((arr_prod == 1) & valid_mask_prod) |
            ((arr_soc == 1) & valid_mask_soc)
    )

    # 2. Improved (3): 没有1，且至少有一个 3
    mask_improved_raw = (
            ((arr_lc == 3) & valid_mask_lc) |
            ((arr_prod == 3) & valid_mask_prod) |
            ((arr_soc == 3) & valid_mask_soc)
    )
    mask_improved = mask_improved_raw & (~mask_degraded) & analyzed_mask

    # 3. Stable (2): 剩下的有效区域
    mask_stable = analyzed_mask & (~mask_degraded) & (~mask_improved)

    final_sdg[mask_stable] = 2
    final_sdg[mask_improved] = 3
    final_sdg[mask_degraded] = 1

    # ================= 边界清洗 =================
    print("5. 应用边界清洗 (去除海洋)...")
    geo_mask = geometry_mask(
        geometries=shp_aligned.geometry,
        out_shape=(min_h, min_w),
        transform=out_transform,
        invert=True
    )
    final_sdg[~geo_mask] = 0

    # ================= 统计输出 =================
    print("6. 计算最终统计数据...")
    area_degraded = calculate_area_weighted(final_sdg, out_transform, final_sdg == 1, crs)
    area_stable = calculate_area_weighted(final_sdg, out_transform, final_sdg == 2, crs)
    area_improved = calculate_area_weighted(final_sdg, out_transform, final_sdg == 3, crs)

    total_area = area_degraded + area_stable + area_improved

    df_result = pd.DataFrame({
        "Category": ["Degraded", "Stable", "Improved", "Total Analyzed"],
        "Area (km2)": [area_degraded, area_stable, area_improved, total_area],
        "Proportion (%)": [
            area_degraded / total_area * 100 if total_area else 0,
            area_stable / total_area * 100 if total_area else 0,
            area_improved / total_area * 100 if total_area else 0,
            100.0
        ]
    })

    print("\n====== SDG 15.3.1 National Estimates (Corrected) ======")
    print(df_result.round(2))
    df_result.to_csv(output_csv, index=False)

    # 保存 TIF
    meta = {
        'driver': 'GTiff',
        'height': min_h,
        'width': min_w,
        'transform': out_transform,
        'crs': crs,
        'count': 1,
        'dtype': 'uint8',
        'nodata': 0,
        'compress': 'lzw'
    }
    with rasterio.open(output_tif, "w", **meta) as dst:
        dst.write(final_sdg, 1)

    print(f"\n结果已保存至: {output_tif}")

    # 绘图
    plt.figure(figsize=(10, 8))
    # 0:Trans, 1:Red, 2:Yellow, 3:Green
    cmap_vals = np.zeros((4, 4))
    cmap_vals[0] = [1, 1, 1, 0]
    cmap_vals[1] = [1, 0, 0, 1]
    cmap_vals[2] = [1, 0.8, 0, 1]
    cmap_vals[3] = [0, 0.5, 0, 1]

    plt.imshow(final_sdg, cmap=ListedColormap(cmap_vals), vmin=0, vmax=3, interpolation='nearest')
    plt.title(f"SDG 15.3.1 Final (Standardized)\nDegraded: {df_result.iloc[0, 2]:.2f}%")
    plt.axis('off')
    plt.show()


if __name__ == "__main__":
    main()