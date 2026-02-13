import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os

# --- 1. 核心配置 ---
# 必须使用 .shp 进行裁剪
path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
data_dir = r"D:\deep\fromFeishu\Land_Cover_Each_Year"

# 基线文件 (2015)
path_baseline = os.path.join(data_dir, "C3S_2018_reclass_300m.tif")

# 要计算的目标年份
target_years = [2019, 2020, 2021, 2022, 2023, 2024]

# 您的因子表 (Factors)
# 索引0设为NaN, 索引1对应Forest(1.0)...
# 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])


def load_and_crop(raster_path, shp_gdf):
    """
    读取 Raster 并用 Shapefile 进行裁剪 (切掉灰色部分)
    """
    with rasterio.open(raster_path) as src:
        # 1. 坐标系对齐
        if str(src.crs) != str(shp_gdf.crs):
            shp_gdf = shp_gdf.to_crs(src.crs)

        # 2. 执行裁剪 (crop=True)
        # 这一步会把 Shapefile 范围外的所有像素设为 NoData
        out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
        data = out_image[0]

        # 获取 NoData 值
        nodata = src.nodata if src.nodata is not None else 0

        return data, out_transform, nodata, src.crs


def calculate_area_weighted(mask_arr, transform, crs):
    """计算面积 (兼容经纬度加权)"""
    res_x = transform[0]
    res_y = -transform[4]

    if crs.to_epsg() == 4326:  # WGS84 经纬度
        height, width = mask_arr.shape
        y_indices = np.arange(height) + 0.5
        latitudes = transform[5] + y_indices * transform[4]
        lat_rad = np.radians(latitudes)

        # 每一行的像素面积 (km2)
        pixel_width_km = np.abs(res_x) * (np.pi / 180) * 6371.0 * np.cos(lat_rad)
        pixel_height_km = np.abs(res_y) * (np.pi / 180) * 6371.0
        row_area = pixel_width_km * pixel_height_km

        count_per_row = np.sum(mask_arr, axis=1)
        return np.sum(count_per_row * row_area)
    else:
        # 投影坐标 (单位米)
        pixel_area_km2 = (abs(res_x * res_y)) / 1_000_000
        return np.sum(mask_arr) * pixel_area_km2


def main():
    print("1. 读取 Shapefile (用于裁剪)...")
    shp_gdf = gpd.read_file(path_shp)

    print("2. 读取基线数据 (2015)...")
    if not os.path.exists(path_baseline):
        print(f"错误: 找不到基线文件 {path_baseline}")
        return

    # 加载并裁剪基线数据
    base_data, transform, base_nodata, crs = load_and_crop(path_baseline, shp_gdf)

    # 预处理基线: 将 Land Cover ID (1,2,3...) 转换为 Factor (1.0, 0.58...)
    base_data_clamped = np.clip(base_data, 0, 7)
    f_base = factor_map[base_data_clamped]

    # 确定基线有效区域 (排除水体 Factor=0 和 NoData)
    # 这是分母，必须大于0
    base_valid_mask = (~np.isnan(f_base)) & (f_base > 0)

    results = []

    print("\n3. 开始批量计算 (2017-2024)...")
    for year in target_years:
        fname = f"C3S_{year}_reclass_300m.tif"
        fpath = os.path.join(data_dir, fname)

        if not os.path.exists(fpath):
            print(f"  [跳过] 文件不存在: {fname}")
            continue

        print(f"  正在计算: 2015 (基线) vs {year} ...")

        # 加载并裁剪当前年份数据
        curr_data, _, curr_nodata, _ = load_and_crop(fpath, shp_gdf)

        # 检查形状
        if base_data.shape != curr_data.shape:
            print("  [错误] 栅格形状不一致，无法计算。")
            continue

        # 1. 转换当前年份 Land Cover -> Factor
        curr_data_clamped = np.clip(curr_data, 0, 7)
        f_curr = factor_map[curr_data_clamped]

        # 2. 确定计算区域 (基线和当前年份都必须有效)
        curr_valid_mask = (~np.isnan(f_curr))
        calc_mask = base_valid_mask & curr_valid_mask & (f_curr >= 0)

        # 3. 计算变化率: (F_year - F_base) / F_base * 100
        soc_change_pct = np.full(base_data.shape, np.nan)

        vals_base = f_base[calc_mask]
        vals_curr = f_curr[calc_mask]

        # 仅对分母大于0的像素计算 (虽已过滤，加保险)
        safe_div = vals_base > 0

        changes = (vals_curr[safe_div] - vals_base[safe_div]) / vals_base[safe_div] * 100

        # 填充结果
        temp_result = np.full(vals_base.shape, np.nan)
        temp_result[safe_div] = changes
        soc_change_pct[calc_mask] = temp_result

        # 4. 判定退化/改良 (阈值 10%)
        # Degraded: Change < -10%
        mask_deg = (soc_change_pct < -10)

        # Improved: Change > 10%
        mask_imp = (soc_change_pct > 10)

        # 5. 统计面积 (km2)
        area_deg = calculate_area_weighted(mask_deg, transform, crs)
        area_imp = calculate_area_weighted(mask_imp, transform, crs)

        print(f"    -> Degraded: {area_deg:.2f} km² | Improved: {area_imp:.2f} km²")

        results.append({
            "Year": year,
            "Baseline": "2018",
            "Area of Degraded Land (km²)": round(area_deg, 2),
            "Area of Improved Land (km²)": round(area_imp, 2)
        })

    # 保存结果
    if results:
        df_res = pd.DataFrame(results)
        out_file = "D:\deep\SoilGrids_Data\main_table\\SOC_Degradation_2019_2024.csv"
        df_res.to_csv(out_file, index=False)
        print(f"\n所有结果已保存至: {out_file}")


if __name__ == "__main__":
    main()