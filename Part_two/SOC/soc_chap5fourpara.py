"""
基于 SDG 15.3.1（土地退化中性/LDN）标准，自动计算阿布扎比（Abu Dhabi）地区在 2019-2024 年间的土地退化和改善的面积统计。

Land Cover (土地覆盖): 读取预处理好的土地覆盖变化栅格，根据 SDG 映射规则（退化、稳定、改善）统计面积
Productivity (土地生产力): 分析生产力变化趋势，统计退化和改善的区域。
SOC (土壤有机碳): 这是一个动态计算模块。它对比 2018 年基准年与当前年份的土地覆盖类型，利用内置的 soc_factor_map（系数表）来估算碳储量变化。如果变化率超过 $\pm 10\%$，则判定为改善或退化。
"""

import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os

# ================= 配置路径 =================
# 1. 基础路径
dir_lc = r"D:\deep\fromFeishu\hmq\Land_Cover_Degradation"
dir_prod = r"D:\deep\fromFeishu\qyx\land_productivity_20260122"
dir_raw_lc = r"D:\deep\fromFeishu\Land_Cover_Each_Year"  # 用于动态计算 SOC

path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
output_csv = r"D:\deep\SoilGrids_Data\chapter5_maintable\Annual_Degradation_Stats_2019_2024.csv"

# 2. 年份列表
years = [2019, 2020, 2021, 2022, 2023, 2024]

# 3. 因子映射 (SOC计算用)
# 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
soc_factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])


# ================= 辅助函数 =================
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

        # 将行面积广播到整个 mask
        row_pixel_area = pixel_width_km * pixel_height_km
        count_per_row = np.sum(mask_array, axis=1)
        return np.sum(count_per_row * row_pixel_area)
    else:
        pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
        return np.sum(mask_array) * pixel_area_km2


def load_and_crop(raster_path, shp_gdf):
    with rasterio.open(raster_path) as src:
        if str(src.crs) != str(shp_gdf.crs):
            shp_gdf = shp_gdf.to_crs(src.crs)
        out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
        return out_image[0], out_transform, src.crs, src.nodata


def get_stats_from_layer(data, transform, crs, mapping_dict=None):
    """
    计算 Degraded 和 Improved 面积
    mapping_dict: {旧值: 新SDG值}。如果为 None，假设已经是标准值(1=Deg, 3=Imp)
    SDG标准: 1=Degraded, 2=Stable, 3=Improved
    """
    # 如果需要重映射
    if mapping_dict:
        # 创建临时数组以免修改原数据
        remapped = np.zeros_like(data)
        for old, new in mapping_dict.items():
            remapped[data == old] = new
        target_data = remapped
    else:
        target_data = data

    # 统计
    mask_deg = (target_data == 1)
    mask_imp = (target_data == 3)

    area_deg = calculate_area_weighted(target_data, transform, mask_deg, crs)
    area_imp = calculate_area_weighted(target_data, transform, mask_imp, crs)

    return area_deg, area_imp


# ================= 主流程 =================
def main():
    print("正在初始化...")
    shp = gpd.read_file(path_shp)
    results = []

    # 预加载 2018 年基准数据 (用于 SOC 计算)
    path_2018_raw = os.path.join(dir_raw_lc, "C3S_2018_reclass_300m.tif")
    if os.path.exists(path_2018_raw):
        print(f"加载 SOC 基准年份: {os.path.basename(path_2018_raw)}")
        base_soc_data, base_soc_transform, base_soc_crs, _ = load_and_crop(path_2018_raw, shp)
    else:
        print("警告: 找不到 2018 原始数据，将跳过 SOC 计算。")
        base_soc_data = None

    for year in years:
        print(f"\n=== 处理年份: {year} ===")

        # ---------------- 1. Land Cover ----------------
        # 文件名模式: status_2018_{year}_deg_stable_imp.tif
        f_lc = f"status_2018_{year}_deg_stable_imp.tif"
        p_lc = os.path.join(dir_lc, f_lc)

        if os.path.exists(p_lc):
            data, trans, crs, _ = load_and_crop(p_lc, shp)
            # 映射: 0(Stable)->2, 1(Improved)->3, 2(Degraded)->1
            mapping = {0: 2, 1: 3, 2: 1}
            deg, imp = get_stats_from_layer(data, trans, crs, mapping)
            results.append({"Year": year, "Indicator": "Land Cover", "Degraded (km2)": deg, "Improved (km2)": imp})
            print(f"  [Land Cover]   Degraded: {deg:.2f}, Improved: {imp:.2f}")
        else:
            print(f"  [Land Cover]   文件缺失: {f_lc}")

        # ---------------- 2. Productivity ----------------
        # 文件名模式: combined_degradation_{year}_table45.tif
        f_prod = f"combined_degradation_{year}_table45.tif"
        p_prod = os.path.join(dir_prod, f_prod)

        if os.path.exists(p_prod):
            data, trans, crs, _ = load_and_crop(p_prod, shp)
            # 映射: 0(Stable)->2, 1(Improved)->3, 2(Degraded)->1
            mapping = {0: 2, 1: 3, 2: 1}
            deg, imp = get_stats_from_layer(data, trans, crs, mapping)
            results.append({"Year": year, "Indicator": "Productivity", "Degraded (km2)": deg, "Improved (km2)": imp})
            print(f"  [Productivity] Degraded: {deg:.2f}, Improved: {imp:.2f}")
        else:
            print(f"  [Productivity] 文件缺失: {f_prod}")

        # ---------------- 3. SOC (动态计算) ----------------
        # 逻辑: (Factor_Year - Factor_2018) / Factor_2018
        if base_soc_data is not None:
            f_raw = f"C3S_{year}_reclass_300m.tif"
            p_raw = os.path.join(dir_raw_lc, f_raw)

            if os.path.exists(p_raw):
                curr_data, _, _, _ = load_and_crop(p_raw, shp)

                # 确保尺寸一致
                if curr_data.shape == base_soc_data.shape:
                    # 清洗
                    valid_mask = (base_soc_data >= 1) & (base_soc_data <= 7) & (curr_data >= 1) & (curr_data <= 7)
                    f_start = soc_factor_map[np.where(valid_mask, base_soc_data, 0)]
                    f_end = soc_factor_map[np.where(valid_mask, curr_data, 0)]

                    calc_mask = (f_start > 0) & (~np.isnan(f_start)) & (~np.isnan(f_end))

                    # 计算变化率
                    soc_change = np.zeros_like(f_start)
                    soc_change[calc_mask] = (f_end[calc_mask] - f_start[calc_mask]) / f_start[calc_mask] * 100

                    # 分类
                    mask_deg = (soc_change < -10) & calc_mask
                    mask_imp = (soc_change > 10) & calc_mask

                    # 面积
                    deg = calculate_area_weighted(curr_data, base_soc_transform, mask_deg, base_soc_crs)
                    imp = calculate_area_weighted(curr_data, base_soc_transform, mask_imp, base_soc_crs)

                    results.append({"Year": year, "Indicator": "SOC", "Degraded (km2)": deg, "Improved (km2)": imp})
                    print(f"  [SOC]          Degraded: {deg:.2f}, Improved: {imp:.2f}")
                else:
                    print(f"  [SOC]          尺寸不匹配，跳过。")
            else:
                print(f"  [SOC]          原始数据缺失: {f_raw}")

    # ================= 保存结果 =================
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(by=["Year", "Indicator"])

        # 整理格式
        df = df[["Year", "Indicator", "Degraded (km2)", "Improved (km2)"]]

        print("\n=== 最终统计结果 ===")
        print(df.round(2))

        df.to_csv(output_csv, index=False)
        print(f"\n结果已保存至: {output_csv}")
    else:
        print("未计算出任何结果。")


if __name__ == "__main__":
    main()