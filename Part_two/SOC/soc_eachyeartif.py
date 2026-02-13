import rasterio
from rasterio.mask import mask
from rasterio.features import geometry_mask
import geopandas as gpd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# ================= 1. 配置路径 =================
# 原始 Land Cover 所在的文件夹
dir_raw_lc = r"D:\deep\fromFeishu\Land_Cover_Each_Year"

# Shapefile 路径
path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"

# 输出文件夹
output_dir = r"D:\deep\SoilGrids_Data\soc_eachyeartif\yearly_soc_tifs"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# # reporting目标年份
# years = [2019, 2020, 2021, 2022, 2023, 2024]

# baseline目标年份
years = [2015]

# SOC 因子映射 (保持不变)
# 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])


# ================= 2. 辅助函数 =================
def load_and_crop(raster_path, shp_gdf):
    """读取并裁剪栅格"""
    with rasterio.open(raster_path) as src:
        if str(src.crs) != str(shp_gdf.crs):
            shp_gdf = shp_gdf.to_crs(src.crs)
        out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
        return out_image[0], out_transform, src.crs, shp_gdf


# ================= 3. 主程序 =================
def main():
    print("--- 开始生成年度 SOC 变化 TIF (自定义标签版) ---")
    print("   标签定义: 0=Stable, 1=Improved, 2=Degraded, 7=NoData")

    # 1. 读取 Shapefile
    try:
        shp = gpd.read_file(path_shp)
    except Exception as e:
        print(f"错误: 无法读取 Shapefile - {e}")
        return

    # 2. 读取并处理基准年份 (2018)
    path_base = os.path.join(dir_raw_lc, "C3S_2000_reclass_300m.tif")
    if not os.path.exists(path_base):
        print(f"错误: 找不到基准年份文件 {path_base}")
        return

    print(f"1. 加载基准年份 (2018)...")
    base_data, base_transform, crs, shp_aligned = load_and_crop(path_base, shp)

    # 生成 2018 的 SOC 因子图
    base_data_safe = base_data.copy()
    base_data_safe[(base_data < 0) | (base_data > 7)] = 0
    soc_base = factor_map[base_data_safe]

    # 3. 循环处理每一年
    for year in years:
        filename = f"C3S_{year}_reclass_300m.tif"
        file_path = os.path.join(dir_raw_lc, filename)

        print(f"\n--- 处理年份: {year} ---")
        if not os.path.exists(file_path):
            print(f"   [警告] 文件不存在，跳过")
            continue

        # 读取当年数据
        curr_data, curr_transform, _, _ = load_and_crop(file_path, shp_aligned)

        # 尺寸检查
        if curr_data.shape != base_data.shape:
            print(f"   [错误] 尺寸不匹配，跳过。")
            continue

        # 转换为 SOC 因子
        curr_data_safe = curr_data.copy()
        curr_data_safe[(curr_data < 0) | (curr_data > 7)] = 0
        soc_curr = factor_map[curr_data_safe]

        # --- 计算变化率 ---
        calc_mask = (soc_base > 0) & (~np.isnan(soc_base)) & (soc_curr > 0) & (~np.isnan(soc_curr))
        soc_change_pct = np.zeros_like(soc_base)
        soc_change_pct[calc_mask] = (soc_curr[calc_mask] - soc_base[calc_mask]) / soc_base[calc_mask] * 100

        # --- 关键修改：自定义标签分类 ---
        # 初始化全为 7 (NoData)
        final_raster = np.full(base_data.shape, 7, dtype=np.uint8)

        # 1. Degraded (2): < -10%
        mask_deg = (soc_change_pct < -10) & calc_mask

        # 2. Improved (1): > 10%
        mask_imp = (soc_change_pct > 10) & calc_mask

        # 3. Stable (0): -10% <= x <= 10%
        mask_stable = (soc_change_pct >= -10) & (soc_change_pct <= 10) & calc_mask

        # 赋值
        final_raster[mask_stable] = 0
        final_raster[mask_imp] = 1
        final_raster[mask_deg] = 2

        # --- 边界清洗 ---
        # Shapefile 外部强制设为 7
        geo_mask = geometry_mask(
            geometries=shp_aligned.geometry,
            out_shape=base_data.shape,
            transform=base_transform,
            invert=True
        )
        final_raster[~geo_mask] = 7

        # --- 保存 ---
        # out_name = f"soc_status_2018_{year}_custom.tif"
        out_name = f"soc_status_2000_{year}_custom.tif"
        out_path = os.path.join(output_dir, out_name)

        meta = {
            'driver': 'GTiff',
            'height': final_raster.shape[0],
            'width': final_raster.shape[1],
            'transform': base_transform,
            'crs': crs,
            'count': 1,
            'dtype': 'uint8',
            'nodata': 7,
            'compress': 'lzw'
        }

        with rasterio.open(out_path, 'w', **meta) as dst:
            dst.write(final_raster, 1)

        print(f"   已生成: {out_name}")

    print("\n所有文件处理完成！")

    # --- 可视化检查 (验证标签颜色) ---
    print("正在生成 2015 年预览图 (验证标签)...")
    plt.figure(figsize=(10, 8))

    # 构建颜色映射以匹配你的标签
    # 0: Stable -> 黄色
    # 1: Improved -> 绿色
    # 2: Degraded -> 红色
    # 3-6: 未使用 (透明)
    # 7: NoData -> 白色

    # 初始化一个 8 色的 colormap (0-7)
    cmap_colors = np.zeros((8, 4))

    cmap_colors[0] = [1, 0.8, 0, 1]  # 0: Stable (Yellow)
    cmap_colors[1] = [0, 0.5, 0, 1]  # 1: Improved (Green)
    cmap_colors[2] = [1, 0, 0, 1]  # 2: Degraded (Red)
    # 3-6 保持透明
    cmap_colors[7] = [1, 1, 1, 0]  # 7: NoData (Transparent/White)

    custom_cmap = ListedColormap(cmap_colors)

    plt.imshow(final_raster, cmap=custom_cmap, vmin=0, vmax=7, interpolation='nearest')
    plt.title("SOC Status Check (Custom Labels: 0=Stable, 1=Imp, 2=Deg)")

    # 自定义图例
    cbar = plt.colorbar(ticks=[0.4, 1.3, 2.2, 6.6], shrink=0.6)
    cbar.ax.set_yticklabels(['0: Stable', '1: Improved', '2: Degraded', '7: NoData'])

    plt.axis('off')
    plt.show()


if __name__ == "__main__":
    main()