import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os

# --- 1. 配置参数 ---
# 包含所有分区的 Shapefile 路径
path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
# 分区名称所在的列名 (您确认是 'name')
name_col = "NAME"

# 栅格数据路径
data_dir = r"D:\deep\fromFeishu\Land_Cover_Each_Year"
path_baseline = os.path.join(data_dir, "C3S_2017_reclass_300m.tif")  # 基线 2017
path_target = os.path.join(data_dir, "C3S_2024_reclass_300m.tif")  # 目标 2019

# 因子映射 (IPCC Tropical Dry)
# 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])


def calculate_area_weighted(mask_arr, transform, crs):
    """计算面积 (兼容 WGS84 经纬度加权 和 投影坐标)"""
    res_x, res_y = transform[0], -transform[4]

    # 判断是否为经纬度坐标 (EPSG:4326)
    if crs.to_epsg() == 4326:
        height, width = mask_arr.shape
        # 计算每一行的中心纬度
        y_indices = np.arange(height) + 0.5
        latitudes = transform[5] + y_indices * transform[4]
        lat_rad = np.radians(latitudes)

        # 计算该行单个像素的面积 (km2)
        # 地球半径 R ≈ 6371 km
        pixel_width_km = np.abs(res_x) * (np.pi / 180) * 6371.0 * np.cos(lat_rad)
        pixel_height_km = np.abs(res_y) * (np.pi / 180) * 6371.0
        row_area = pixel_width_km * pixel_height_km

        # 统计每一行 True 的数量并乘以该行像素面积
        count_per_row = np.sum(mask_arr, axis=1)
        return np.sum(count_per_row * row_area)
    else:
        # 投影坐标 (假设单位是米)
        pixel_area_km2 = (abs(res_x * res_y)) / 1_000_000
        return np.sum(mask_arr) * pixel_area_km2


def process_region(region_name, geometry, src_base, src_target):
    """处理单个区域的 SOC 变化计算"""

    # 1. 裁剪基线 (2015)
    try:
        out_img_base, out_trans_base = mask(src_base, [geometry], crop=True)
        data_base = out_img_base[0]
    except ValueError:
        # 几何体不在栅格范围内
        return None

    # 2. 裁剪目标年 (2019)
    try:
        out_img_target, _ = mask(src_target, [geometry], crop=True)
        data_target = out_img_target[0]
    except ValueError:
        return None

    # 形状检查
    if data_base.shape != data_target.shape:
        print(f"  [警告] {region_name} 裁剪后形状不匹配，跳过。")
        return None

    # 3. 映射为因子
    # clip防止索引越界 (虽然理论上都是1-7)
    f_base = factor_map[np.clip(data_base, 0, 7)]
    f_target = factor_map[np.clip(data_target, 0, 7)]

    # 4. 确定有效计算区域
    # 基线必须有效且非水体 (>0)，目标年必须有效
    valid_mask = (~np.isnan(f_base)) & (f_base > 0) & (~np.isnan(f_target))

    if np.sum(valid_mask) == 0:
        return 0.0, 0.0

    # 5. 计算变化率 (Pixel-based)
    # 先初始化全图为 NaN
    change_pct = np.full(data_base.shape, np.nan)

    # 提取有效值计算
    v_base = f_base[valid_mask]
    v_target = f_target[valid_mask]

    # 计算公式: (Target - Base) / Base * 100
    calc_values = (v_target - v_base) / v_base * 100

    # 填回矩阵
    change_pct[valid_mask] = calc_values

    # 6. 判定退化/改良 (阈值 10%)
    mask_deg = (change_pct < -10)
    mask_imp = (change_pct > 10)

    # 7. 计算面积
    area_deg = calculate_area_weighted(mask_deg, out_trans_base, src_base.crs)
    area_imp = calculate_area_weighted(mask_imp, out_trans_base, src_base.crs)

    return area_deg, area_imp


def main():
    print("1. 读取 Shapefile...")
    try:
        gdf = gpd.read_file(path_shp)
        print(f"   成功读取。共包含 {len(gdf)} 个区域。")

        if name_col not in gdf.columns:
            print(f"   [错误] 列名 '{name_col}' 不存在。可用列名: {gdf.columns.tolist()}")
            return

    except Exception as e:
        print(f"   [错误] 读取 Shapefile 失败: {e}")
        return

    print("2. 准备栅格数据...")
    if not os.path.exists(path_baseline) or not os.path.exists(path_target):
        print("   [错误] 找不到 2015 或 2019 的 TIF 文件，请检查路径。")
        return

    results = []

    # 保持栅格文件打开，提高效率
    with rasterio.open(path_baseline) as src_base, rasterio.open(path_target) as src_target:

        # 坐标系对齐: 强制将 Shapefile 转为 Raster 的坐标系
        if str(gdf.crs) != str(src_base.crs):
            print(f"   正在将 Shapefile 坐标系转为 {src_base.crs} ...")
            gdf = gdf.to_crs(src_base.crs)

        print("\n3. 开始分区域计算 (2019 vs 2015)...")

        # 遍历 Shapefile 中的每一行
        for idx, row in gdf.iterrows():
            region_name = row[name_col]
            geometry = row.geometry

            print(f"   正在处理: {region_name} ...")

            res = process_region(region_name, geometry, src_base, src_target)

            if res:
                deg, imp = res
                print(f"     -> 退化: {deg:.2f} km² | 改良: {imp:.2f} km²")
                results.append({
                    "Emirate": region_name,
                    "Year": 2019,
                    "Area of Degraded Land (km²)": round(deg, 2),
                    "Area of Improved Land (km²)": round(imp, 2)
                })
            else:
                print("     -> 无有效数据 (可能区域在栅格范围外)")
                results.append({
                    "Emirate": region_name,
                    "Year": 2019,
                    "Area of Degraded Land (km²)": 0,
                    "Area of Improved Land (km²)": 0
                })

    # --- 输出结果 ---
    if results:
        df = pd.DataFrame(results)
        # 按名称排序
        df = df.sort_values(by="Emirate")

        print("\n=== Areal distribution of degraded and improved land (2019) ===")
        print(df)

        out_csv = "D:\deep\SoilGrids_Data\main_table\\SOC_Distribution_By_Region_2024.csv"
        #df.to_csv(out_csv, index=False)
        print(f"\n结果已保存至: {out_csv}")
    else:
        print("未生成任何结果。")


if __name__ == "__main__":
    main()