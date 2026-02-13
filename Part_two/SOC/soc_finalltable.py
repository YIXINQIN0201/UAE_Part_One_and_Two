"""
采用联合国可持续发展目标（SDG）规定的 “一票否决制” (One Out, All Out, OOAO) 逻辑，
将之前计算的土地覆盖（LC）、生产力（Productivity）和土壤有机碳（SOC）三个子指标合而为一，
生成最终的 15.3.1 土地退化综合统计结果


"""

import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os

# ================= 1. 配置路径 =================
# 输入文件夹
dir_lc = r"D:\deep\fromFeishu\hmq\Land_Cover_Degradation"
dir_prod = r"D:\deep\fromFeishu\qyx\land_productivity_20260122"
# 注意：这里使用刚才生成的“自定义标签版”SOC文件夹
dir_soc = r"D:\deep\SoilGrids_Data\soc_eachyeartif\yearly_soc_tifs"

# Shapefile (用于面积加权和边界限制)
path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"

# 输出文件
output_csv = r"D:\deep\SoilGrids_Data\chapter5_maintable\Final_SDG_15_3_1_Annual_Summary_2019_2024.csv"

# 年份
years = [2019, 2020, 2021, 2022, 2023, 2024]

# 你的自定义标签值
VAL_STABLE = 0
VAL_IMPROVED = 1
VAL_DEGRADED = 2
VAL_NODATA = 7


# ================= 2. 辅助函数 =================
def load_and_align(raster_path, shp_gdf, target_shape=None):
    """读取栅格，裁剪并强制对齐"""
    with rasterio.open(raster_path) as src:
        # 坐标系对齐
        if str(src.crs) != str(shp_gdf.crs):
            shp_gdf = shp_gdf.to_crs(src.crs)

        # 裁剪
        out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
        data = out_image[0]

        # 形状对齐 (处理1像素误差)
        if target_shape is not None and data.shape != target_shape:
            h = min(data.shape[0], target_shape[0])
            w = min(data.shape[1], target_shape[1])
            data = data[:h, :w]

        return data, out_transform, src.crs, shp_gdf


def calculate_area_weighted(data, transform, mask_array, crs):
    """计算加权面积 (km²)"""
    if np.sum(mask_array) == 0:
        return 0.0

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

        # 广播到 mask
        row_pixel_area = pixel_width_km * pixel_height_km
        count_per_row = np.sum(mask_array, axis=1)
        return np.sum(count_per_row * row_pixel_area)
    else:
        pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
        return np.sum(mask_array) * pixel_area_km2


# ================= 3. 主流程 =================
def main():
    print("--- 开始计算 SDG 15.3.1 综合指标 (2019-2024) ---")
    print(f"   标签标准: 0=Stable, 1=Improved, 2=Degraded")

    # 1. 读取 Shapefile
    try:
        shp = gpd.read_file(path_shp)
    except Exception as e:
        print(f"错误: 无法读取 Shapefile - {e}")
        return

    stats_list = []

    for year in years:
        print(f"\n=== 正在处理年份: {year} ===")

        # --- 构建文件名 ---
        # Land Cover: status_2018_{year}_deg_stable_imp.tif
        f_lc = f"status_2018_{year}_deg_stable_imp.tif"
        p_lc = os.path.join(dir_lc, f_lc)

        # Productivity: combined_degradation_{year}_table45.tif
        f_prod = f"combined_degradation_{year}_table45.tif"
        p_prod = os.path.join(dir_prod, f_prod)

        # SOC: soc_status_2018_{year}_custom.tif
        f_soc = f"soc_status_2018_{year}_custom.tif"
        p_soc = os.path.join(dir_soc, f_soc)

        # --- 检查文件存在性 ---
        if not (os.path.exists(p_lc) and os.path.exists(p_prod) and os.path.exists(p_soc)):
            print(f"   [警告] 缺失文件，跳过 {year} 年")
            print(
                f"   LC存在: {os.path.exists(p_lc)}, Prod存在: {os.path.exists(p_prod)}, SOC存在: {os.path.exists(p_soc)}")
            continue

        # --- 加载数据 (以 LC 为形状基准) ---
        arr_lc, trans, crs, shp_aligned = load_and_align(p_lc, shp)
        base_shape = arr_lc.shape

        arr_prod, _, _, _ = load_and_align(p_prod, shp_aligned, target_shape=base_shape)
        arr_soc, _, _, _ = load_and_align(p_soc, shp_aligned, target_shape=base_shape)

        # 统一裁剪到最小尺寸
        min_h = min(arr_lc.shape[0], arr_prod.shape[0], arr_soc.shape[0])
        min_w = min(arr_lc.shape[1], arr_prod.shape[1], arr_soc.shape[1])

        arr_lc = arr_lc[:min_h, :min_w]
        arr_prod = arr_prod[:min_h, :min_w]
        arr_soc = arr_soc[:min_h, :min_w]

        # --- 数据预处理 ---
        # 确保只处理有效值 (0, 1, 2)。其他值视为无效。
        # 你的 SOC 文件 NoData 是 7， LC/Prod 可能也有其他值
        valid_lc = np.isin(arr_lc, [VAL_STABLE, VAL_IMPROVED, VAL_DEGRADED])
        valid_prod = np.isin(arr_prod, [VAL_STABLE, VAL_IMPROVED, VAL_DEGRADED])
        valid_soc = np.isin(arr_soc, [VAL_STABLE, VAL_IMPROVED, VAL_DEGRADED])

        # 只要任意一个图层有有效数据，我们就在这个像素进行计算
        analyzed_mask = valid_lc | valid_prod | valid_soc

        # --- 应用 One Out, All Out 逻辑 ---

        # 1. 判断 Degraded (2)
        # 逻辑：任意一个有效图层显示为 Degraded
        is_degraded = (
                ((arr_lc == VAL_DEGRADED) & valid_lc) |
                ((arr_prod == VAL_DEGRADED) & valid_prod) |
                ((arr_soc == VAL_DEGRADED) & valid_soc)
        )

        # 2. 判断 Improved (1)
        # 逻辑：不是 Degraded，且至少有一个有效图层显示为 Improved
        is_improved_raw = (
                ((arr_lc == VAL_IMPROVED) & valid_lc) |
                ((arr_prod == VAL_IMPROVED) & valid_prod) |
                ((arr_soc == VAL_IMPROVED) & valid_soc)
        )
        # 必须排除掉已经标记为 degraded 的区域
        is_improved = is_improved_raw & (~is_degraded) & analyzed_mask

        # 3. 判断 Stable (0)
        # 逻辑：在分析区域内，既不是 Degraded 也不是 Improved
        is_stable = analyzed_mask & (~is_degraded) & (~is_improved)

        # --- 计算面积 ---
        area_deg = calculate_area_weighted(arr_lc, trans, is_degraded, crs)
        area_imp = calculate_area_weighted(arr_lc, trans, is_improved, crs)
        area_stab = calculate_area_weighted(arr_lc, trans, is_stable, crs)

        total = area_deg + area_imp + area_stab

        print(f"   [结果] Degraded: {area_deg:.2f} km², Improved: {area_imp:.2f} km², Stable: {area_stab:.2f} km²")

        stats_list.append({
            "Year": year,
            "Area of Degraded Land (km2)": area_deg,
            "Area of Improved Land (km2)": area_imp,
            "Area of Stable Land (km2)": area_stab,
            "Total Analyzed Area (km2)": total,
            "Percent Degraded (%)": (area_deg / total * 100) if total > 0 else 0
        })

    # ================= 保存结果 =================
    if stats_list:
        df = pd.DataFrame(stats_list)
        print("\n=== 最终汇总表 (SDG 15.3.1) ===")
        print(df.round(2))

        df.to_csv(output_csv, index=False)
        print(f"\n结果已保存至: {output_csv}")
    else:
        print("\n没有生成任何数据，请检查文件路径或文件名是否正确。")


if __name__ == "__main__":
    main()