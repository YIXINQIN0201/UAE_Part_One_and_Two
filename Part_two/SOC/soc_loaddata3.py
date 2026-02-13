# # # baseline period
# import rasterio
# from rasterio.mask import mask
# import geopandas as gpd
# import numpy as np
# import pandas as pd
# import os
#
# # --- 1. 配置路径 ---
# # 你的 Shapefile 路径 (用于裁剪)
# path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
#
# # 你的 Raster 路径
# path_2000 = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2000_reclass_300m.tif"
# path_2015 = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2015_reclass_300m.tif"
#
# # 因子映射 (保持不变)
# # 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
# factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])
#
#
# def load_and_crop_raster(raster_path, shp_gdf):
#     """
#     读取并根据 Shapefile 裁剪 Raster
#     返回: 裁剪后的数据(numpy), 新的Transform, 新的NoData值
#     """
#     with rasterio.open(raster_path) as src:
#         # 1. 检查坐标系并重投影 Shapefile (如果需要)
#         if str(src.crs) != str(shp_gdf.crs):
#             print(f"检测到坐标系不一致，正在将 Shapefile 从 {shp_gdf.crs} 转换为 {src.crs} ...")
#             shp_gdf = shp_gdf.to_crs(src.crs)
#
#         # 2. 执行裁剪 (crop=True 表示把图片尺寸也切小到 Shapefile 的边界)
#         # 这一步会自动把 Shapefile 外面的像素设为 src.nodata
#         out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
#
#         # out_image 形状是 (Bands, Height, Width)，我们只需要第一波段
#         data = out_image[0]
#
#         # 获取 NoData 值 (如果源文件没定义，通常设为 0 或 255，这里要小心)
#         nodata_val = src.nodata
#         if nodata_val is None:
#             nodata_val = 0  # 假设 0 是背景
#
#         return data, out_transform, nodata_val, src.crs
#
#
# def calculate_accurate_area_masked(data, transform, mask_array, crs):
#     """
#     计算裁剪后数据的真实面积 (纬度加权)
#     """
#     # 获取分辨率
#     res_x = transform[0]  # 像素宽度 (度)
#     res_y = -transform[4]  # 像素高度 (度，通常Transform里是负数代表向下)
#
#     height, width = data.shape
#
#     # 判断是否为 WGS84 (经纬度)
#     is_wgs84 = (crs.to_epsg() == 4326)
#
#     if is_wgs84:
#         # 地球半径 (km)
#         R = 6371.0
#         # 计算每一行的中心纬度
#         # Transform: x_origin, pixel_width, 0, y_origin, 0, pixel_height
#         y_origin = transform[5]
#         pixel_height = transform[4]  # 负值
#
#         # 生成每一行的行号
#         y_indices = np.arange(height) + 0.5
#         latitudes = y_origin + y_indices * pixel_height
#
#         # 转弧度
#         lat_rad = np.radians(latitudes)
#
#         # 计算每一行的单像素面积 (km²)
#         pixel_width_km = np.abs(res_x) * (np.pi / 180) * R * np.cos(lat_rad)
#         pixel_height_km = np.abs(res_y) * (np.pi / 180) * R
#         row_pixel_area = pixel_width_km * pixel_height_km
#
#         # 统计 mask 中每一行的有效像素数，乘以该行的单像素面积
#         count_per_row = np.sum(mask_array, axis=1)
#         total_area = np.sum(count_per_row * row_pixel_area)
#         return total_area
#
#     else:
#         # 投影坐标，直接乘
#         pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
#         return np.sum(mask_array) * pixel_area_km2
#
#
# def main_workflow():
#     print("1. 正在读取 Shapefile...")
#     try:
#         shp_gdf = gpd.read_file(path_shp)
#     except Exception as e:
#         print(f"读取 Shapefile 失败，请检查路径或安装 geopandas: {e}")
#         return
#
#     print("2. 正在裁剪 2000年 数据...")
#     data_2000, transform_2000, nodata_2000, crs = load_and_crop_raster(path_2000, shp_gdf)
#
#     print("3. 正在裁剪 2015年 数据...")
#     data_2015, transform_2015, nodata_2015, _ = load_and_crop_raster(path_2015, shp_gdf)
#
#     # 确保两个数组形状一致 (裁剪后应该是一样的)
#     if data_2000.shape != data_2015.shape:
#         print("错误：两年份数据裁剪后形状不一致，无法计算。")
#         return
#
#     print(f"裁剪后矩阵形状: {data_2000.shape}")
#
#     # --- 数据清理 ---
#     # 将 NoData 区域、以及 0 值区域 (通常是背景) 统一标记为无效
#     # 注意：你的类别是 1-7。如果裁剪后出现的 NoData 是 0 或 255，需要处理防止索引越界
#
#     # 设定有效范围 mask: 值必须在 1-7 之间
#     valid_range_2000 = (data_2000 >= 1) & (data_2000 <= 7)
#     valid_range_2015 = (data_2015 >= 1) & (data_2015 <= 7)
#
#     # 组合有效性: 必须都在 Shapefile 内部且是有效类别
#     # 注意：crop=True 后，shapefile 外面被填成了 nodata，这里会自动过滤掉
#     final_valid_pixels = valid_range_2000 & valid_range_2015
#
#     # 准备索引数组 (防止脏数据导致 index error)
#     # 创建副本进行计算
#     idx_2000 = data_2000.copy()
#     idx_2015 = data_2015.copy()
#
#     # 把无效值设为 0 (对应 factor_map[0] = NaN)
#     idx_2000[~final_valid_pixels] = 0
#     idx_2015[~final_valid_pixels] = 0
#
#     # --- 映射因子 ---
#     f_start = factor_map[idx_2000]
#     f_end = factor_map[idx_2015]
#
#     # --- 计算变化率 ---
#     # 核心计算逻辑：只算因子 > 0 的区域 (排除水体和背景)
#     calc_mask = (f_start > 0) & (~np.isnan(f_start)) & (~np.isnan(f_end))
#
#     soc_change = np.full(data_2000.shape, np.nan)
#     # 避免除以 0
#     soc_change[calc_mask] = (f_end[calc_mask] - f_start[calc_mask]) / f_start[calc_mask] * 100
#
#     # --- 分类统计 ---
#     # Degraded: < -10%
#     mask_degraded = (soc_change < -10) & calc_mask
#
#     # Non-degraded: >= -10%
#     mask_non_degraded = (soc_change >= -10) & calc_mask
#
#     # No Data:
#     # 这里的定义是：Shapefile 范围内，但无法计算 SOC 变化的区域 (如水体、或者原始数据缺失)
#     # 我们用 shapefile 的几何掩膜来确定 "Total Land Area" 应该是多少
#     # rasterio.mask 裁剪出的矩阵是矩形的，包含 shapefile 边界外的 nodata。
#     # 我们只关心 data_2000 != nodata_2000 的部分 (即 Shapefile 内部)
#
#     inside_shp_mask = (data_2000 != nodata_2000)
#
#     # 修正 NoData 定义：在 Shapefile 内部，但没有被归类为 Degraded 或 Non-degraded 的
#     mask_nodata_area = inside_shp_mask & (~mask_degraded) & (~mask_non_degraded)
#
#     print("4. 正在计算真实面积...")
#     area_degraded = calculate_accurate_area_masked(data_2000, transform_2000, mask_degraded, crs)
#     area_non_degraded = calculate_accurate_area_masked(data_2000, transform_2000, mask_non_degraded, crs)
#     area_nodata = calculate_accurate_area_masked(data_2000, transform_2000, mask_nodata_area, crs)
#
#     total_area = area_degraded + area_non_degraded + area_nodata
#
#     # --- 输出表 ---
#     df = pd.DataFrame({
#         "Category": ["Degraded SOC", "Non-degraded SOC", "No Data / Water"],
#         "Area (km²)": [area_degraded, area_non_degraded, area_nodata],
#         "Percent (%)": [
#             area_degraded / total_area * 100 if total_area else 0,
#             area_non_degraded / total_area * 100 if total_area else 0,
#             area_nodata / total_area * 100 if total_area else 0
#         ]
#     })
#
#     print("\n=== National estimates of SOC stock degradation (Cropped by Shapefile) ===")
#     print(df.round(2))
#     print(f"Total Area analyzed: {total_area:.2f} km²")
#     df.to_csv("D:\deep\SoilGrids_Data\\appendix_table3\\SOC_baseline_2000_2015.csv", index=False)
#
#
# if __name__ == "__main__":
#     main_workflow()

#=========================================================================做baseline期的csv和tif
# import rasterio
# from rasterio.mask import mask
# import geopandas as gpd
# import numpy as np
# import pandas as pd
# import os
#
# # --- 1. 配置路径 ---
# # 你的 Shapefile 路径 (用于裁剪)
# path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
#
# # 你的 Raster 路径 (Baseline: 2000 vs 2015)
# path_2000 = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2000_reclass_300m.tif"
# path_2015 = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2015_reclass_300m.tif"
#
# # 输出 TIF 文件的路径
# output_tif_path = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2000_2015.tif"
# # 输出 CSV 文件的路径
# output_csv_path = r"D:\deep\SoilGrids_Data\appendix_table3\SOC_baseline_2000_2015.csv"
#
# # 因子映射 (保持不变)
# # 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
# factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])
#
#
# def load_and_crop_raster(raster_path, shp_gdf):
#     """
#     读取并根据 Shapefile 裁剪 Raster
#     返回: 裁剪后的数据(numpy), 新的Transform, 新的NoData值, CRS
#     """
#     with rasterio.open(raster_path) as src:
#         # 1. 检查坐标系并重投影 Shapefile (如果需要)
#         if str(src.crs) != str(shp_gdf.crs):
#             print(f"检测到坐标系不一致，正在将 Shapefile 从 {shp_gdf.crs} 转换为 {src.crs} ...")
#             shp_gdf = shp_gdf.to_crs(src.crs)
#
#         # 2. 执行裁剪 (crop=True 表示把图片尺寸也切小到 Shapefile 的边界)
#         out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
#
#         # out_image 形状是 (Bands, Height, Width)，我们只需要第一波段
#         data = out_image[0]
#
#         # 获取 NoData 值
#         nodata_val = src.nodata
#         if nodata_val is None:
#             nodata_val = 0
#
#         return data, out_transform, nodata_val, src.crs
#
#
# def calculate_accurate_area_masked(data, transform, mask_array, crs):
#     """
#     计算裁剪后数据的真实面积 (纬度加权)
#     """
#     res_x = transform[0]
#     res_y = -transform[4]
#     height, width = data.shape
#     is_wgs84 = (crs.to_epsg() == 4326)
#
#     if is_wgs84:
#         R = 6371.0
#         y_origin = transform[5]
#         pixel_height = transform[4]
#         y_indices = np.arange(height) + 0.5
#         latitudes = y_origin + y_indices * pixel_height
#         lat_rad = np.radians(latitudes)
#
#         pixel_width_km = np.abs(res_x) * (np.pi / 180) * R * np.cos(lat_rad)
#         pixel_height_km = np.abs(res_y) * (np.pi / 180) * R
#         row_pixel_area = pixel_width_km * pixel_height_km
#
#         count_per_row = np.sum(mask_array, axis=1)
#         total_area = np.sum(count_per_row * row_pixel_area)
#         return total_area
#     else:
#         pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
#         return np.sum(mask_array) * pixel_area_km2
#
#
# def main_workflow():
#     print("1. 正在读取 Shapefile...")
#     try:
#         shp_gdf = gpd.read_file(path_shp)
#     except Exception as e:
#         print(f"读取 Shapefile 失败: {e}")
#         return
#
#     print("2. 正在裁剪 2000年 数据...")
#     data_2000, transform_2000, nodata_2000, crs = load_and_crop_raster(path_2000, shp_gdf)
#
#     print("3. 正在裁剪 2015年 数据...")
#     data_2015, transform_2015, nodata_2015, _ = load_and_crop_raster(path_2015, shp_gdf)
#
#     if data_2000.shape != data_2015.shape:
#         print("错误：两年份数据裁剪后形状不一致，无法计算。")
#         return
#
#     print(f"裁剪后矩阵形状: {data_2000.shape}")
#
#     # --- 数据清理 ---
#     valid_range_2000 = (data_2000 >= 1) & (data_2000 <= 7)
#     valid_range_2015 = (data_2015 >= 1) & (data_2015 <= 7)
#     final_valid_pixels = valid_range_2000 & valid_range_2015
#
#     idx_2000 = data_2000.copy()
#     idx_2015 = data_2015.copy()
#     idx_2000[~final_valid_pixels] = 0
#     idx_2015[~final_valid_pixels] = 0
#
#     # --- 映射因子 ---
#     f_start = factor_map[idx_2000]
#     f_end = factor_map[idx_2015]
#
#     # --- 计算变化率 ---
#     # 核心计算逻辑：只算因子 > 0 的区域
#     calc_mask = (f_start > 0) & (~np.isnan(f_start)) & (~np.isnan(f_end))
#
#     soc_change = np.full(data_2000.shape, np.nan)
#     soc_change[calc_mask] = (f_end[calc_mask] - f_start[calc_mask]) / f_start[calc_mask] * 100
#
#     # --- 分类逻辑 ---
#     # Degraded: < -10%
#     mask_degraded = (soc_change < -10) & calc_mask
#     # Non-degraded: >= -10%
#     mask_non_degraded = (soc_change >= -10) & calc_mask
#
#     inside_shp_mask = (data_2000 != nodata_2000)
#     mask_nodata_area = inside_shp_mask & (~mask_degraded) & (~mask_non_degraded)
#
#     # --- [新增功能] 生成并保存 TIF 文件 ---
#     print(f"4. 正在生成状态 TIF 文件: {output_tif_path}")
#
#     # 创建一个全0的矩阵，类型为 uint8 (节省空间)
#     # 0: No Data
#     status_raster = np.zeros(data_2000.shape, dtype=np.uint8)
#
#     # 赋值 (Coding Scheme):
#     # 1: Degraded
#     # 2: Non-degraded (Stable)
#     status_raster[mask_degraded] = 1
#     status_raster[mask_non_degraded] = 2
#
#     # 确保输出目录存在
#     out_dir = os.path.dirname(output_tif_path)
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
#
#     # 定义 Metadata
#     out_meta = {
#         'driver': 'GTiff',
#         'height': status_raster.shape[0],
#         'width': status_raster.shape[1],
#         'transform': transform_2000,
#         'crs': crs,
#         'count': 1,  # 1个波段
#         'dtype': 'uint8',  # 数据类型
#         'nodata': 0,  # NoData 值
#         'compress': 'lzw'  # 压缩格式
#     }
#
#     # 写入文件
#     with rasterio.open(output_tif_path, "w", **out_meta) as dest:
#         dest.write(status_raster, 1)
#
#     print("   -> TIF 文件保存成功。分类代码: 1=Degraded, 2=Non-degraded, 0=No Data")
#
#     # --- 计算面积与输出表格 ---
#     print("5. 正在计算真实面积...")
#     area_degraded = calculate_accurate_area_masked(data_2000, transform_2000, mask_degraded, crs)
#     area_non_degraded = calculate_accurate_area_masked(data_2000, transform_2000, mask_non_degraded, crs)
#     area_nodata = calculate_accurate_area_masked(data_2000, transform_2000, mask_nodata_area, crs)
#
#     total_area = area_degraded + area_non_degraded + area_nodata
#
#     # --- 输出表 ---
#     df = pd.DataFrame({
#         "Category": ["Degraded SOC", "Non-degraded SOC", "No Data / Water"],
#         "Area (km²)": [area_degraded, area_non_degraded, area_nodata],
#         "Percent (%)": [
#             area_degraded / total_area * 100 if total_area else 0,
#             area_non_degraded / total_area * 100 if total_area else 0,
#             area_nodata / total_area * 100 if total_area else 0
#         ]
#     })
#
#     print("\n=== National estimates of SOC stock degradation (Cropped by Shapefile) ===")
#     print(df.round(2))
#     print(f"Total Area analyzed: {total_area:.2f} km²")
#
#     #df.to_csv(output_csv_path, index=False)
#     #print(f"CSV 表格已保存至: {output_csv_path}")
#
#
# if __name__ == "__main__":
#     main_workflow()


#=================================================================================做reporting期的csv
# # reporting period
# import rasterio
# from rasterio.mask import mask
# import geopandas as gpd
# import numpy as np
# import pandas as pd
# import os
#
# # --- 1. 配置路径 (请确认 2018 和 2024 的文件名是否正确) ---
# # Shapefile 路径
# path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
#
# # Reporting Period 的起止年份文件
# path_start_year = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2018_reclass_300m.tif"
# path_end_year = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2024_reclass_300m.tif"
#
# # 输出路径 (新增)
# output_tif_path = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2018_2024_v2.tif"
#
# # 因子映射 (保持不变)
# # 0:NoData, 1:Forest, 2:Grassland, 3:Cropland, 4:Wetland, 5:Artificial, 6:Other, 7:Water
# factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])
#
#
# def load_and_crop_raster(raster_path, shp_gdf):
#     """
#     读取并根据 Shapefile 裁剪 Raster
#     """
#     with rasterio.open(raster_path) as src:
#         # 1. 检查坐标系并重投影 Shapefile
#         if str(src.crs) != str(shp_gdf.crs):
#             print(f"检测到坐标系不一致，正在转换 Shapefile 坐标系...")
#             shp_gdf = shp_gdf.to_crs(src.crs)
#
#         # 2. 执行裁剪 (crop=True)
#         out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
#         data = out_image[0]
#
#         # 获取 NoData 值
#         nodata_val = src.nodata
#         if nodata_val is None:
#             nodata_val = 0
#
#         return data, out_transform, nodata_val, src.crs
#
#
# def calculate_accurate_area_masked(data, transform, mask_array, crs):
#     """
#     计算裁剪后数据的真实面积 (纬度加权校正)
#     """
#     res_x = transform[0]
#     res_y = -transform[4]
#     height, width = data.shape
#
#     # WGS84 经纬度判断
#     is_wgs84 = (crs.to_epsg() == 4326)
#
#     if is_wgs84:
#         R = 6371.0
#         y_origin = transform[5]
#         pixel_height = transform[4]
#
#         y_indices = np.arange(height) + 0.5
#         latitudes = y_origin + y_indices * pixel_height
#         lat_rad = np.radians(latitudes)
#
#         pixel_width_km = np.abs(res_x) * (np.pi / 180) * R * np.cos(lat_rad)
#         pixel_height_km = np.abs(res_y) * (np.pi / 180) * R
#         row_pixel_area = pixel_width_km * pixel_height_km
#
#         count_per_row = np.sum(mask_array, axis=1)
#         total_area = np.sum(count_per_row * row_pixel_area)
#         return total_area
#     else:
#         # 投影坐标 (米)
#         pixel_area_km2 = (abs(res_x) * abs(res_y)) / 1_000_000
#         return np.sum(mask_array) * pixel_area_km2
#
#
# def main_workflow():
#     print("1. 正在读取 Shapefile...")
#     try:
#         shp_gdf = gpd.read_file(path_shp)
#     except Exception as e:
#         print(f"错误: {e}")
#         return
#
#     print(
#         f"2. 正在裁剪 Reporting Period 数据 ({os.path.basename(path_start_year)} -> {os.path.basename(path_end_year)})...")
#
#     # 加载数据
#     data_start, transform, nodata_start, crs = load_and_crop_raster(path_start_year, shp_gdf)
#     data_end, _, nodata_end, _ = load_and_crop_raster(path_end_year, shp_gdf)
#
#     if data_start.shape != data_end.shape:
#         print("错误：两年份数据形状不一致。")
#         return
#
#     # --- 数据清理 ---
#     # 有效像素必须在 1-7 之间
#     valid_range_start = (data_start >= 1) & (data_start <= 7)
#     valid_range_end = (data_end >= 1) & (data_end <= 7)
#     final_valid_pixels = valid_range_start & valid_range_end
#
#     idx_start = data_start.copy()
#     idx_end = data_end.copy()
#     idx_start[~final_valid_pixels] = 0
#     idx_end[~final_valid_pixels] = 0
#
#     # 映射因子
#     f_start = factor_map[idx_start]
#     f_end = factor_map[idx_end]
#
#     # --- 计算变化率 ---
#     # 仅当起始因子 > 0 时计算 (排除水体和背景)
#     calc_mask = (f_start > 0) & (~np.isnan(f_start)) & (~np.isnan(f_end))
#
#     soc_change = np.full(data_start.shape, np.nan)
#     soc_change[calc_mask] = (f_end[calc_mask] - f_start[calc_mask]) / f_start[calc_mask] * 100
#
#     # --- 分类统计 (四类) ---
#     # 1. Degraded: < -10%
#     mask_degraded = (soc_change < -10) & calc_mask
#
#     # 2. Stable: -10% <= change <= 10%
#     mask_stable = (soc_change >= -10) & (soc_change <= 10) & calc_mask
#
#     # 3. Improved: > 10%
#     mask_improved = (soc_change > 10) & calc_mask
#
#     # 4. No Data:
#     # 定义：在 Shapefile 内部 (data_start != nodata)，但不在上述三个mask里的区域
#     # 这包括了水体 (Factor=0) 和原始数据本身就是 NoData 的部分
#     inside_shp_mask = (data_start != nodata_start)
#     mask_nodata_area = inside_shp_mask & (~mask_degraded) & (~mask_stable) & (~mask_improved)
#
#     print("3. 正在计算各类面积...")
#     area_degraded = calculate_accurate_area_masked(data_start, transform, mask_degraded, crs)
#     area_stable = calculate_accurate_area_masked(data_start, transform, mask_stable, crs)
#     area_improved = calculate_accurate_area_masked(data_start, transform, mask_improved, crs)
#     area_nodata = calculate_accurate_area_masked(data_start, transform, mask_nodata_area, crs)
#
#     total_area = area_degraded + area_stable + area_improved + area_nodata
#
#     # --- 输出表 ---
#     df = pd.DataFrame({
#         "Category": [
#             "Land area with degraded SOC",
#             "Land area with stable SOC",
#             "Land area with improved SOC",
#             "Land area with no SOC data"
#         ],
#         "Area (km²)": [area_degraded, area_stable, area_improved, area_nodata],
#         "Percent (%)": [
#             area_degraded / total_area * 100 if total_area else 0,
#             area_stable / total_area * 100 if total_area else 0,
#             area_improved / total_area * 100 if total_area else 0,
#             area_nodata / total_area * 100 if total_area else 0
#         ]
#     })
#
#     print("\n=== Reporting Period (2018-2024) SOC Estimates ===")
#     print(df.round(2))
#     print(f"Total Area analyzed: {total_area:.2f} km²")
#     #df.to_csv("D:\deep\SoilGrids_Data\\appendix_table3\\SOC_Reporting_Period_2018_2024.csv", index=False)
#
#
# if __name__ == "__main__":
#     main_workflow()

#============================================================================绘制2018-2024的.tif图
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np
import pandas as pd
import os

# --- 1. 配置路径 ---
path_shp = r"D:\deep\fromFeishu\abu_dhabi_all\abu_dhabi_all.shp"
path_start_year = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2018_reclass_300m.tif"
path_end_year = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2024_reclass_300m.tif"

# 输出路径 (新增)
output_tif_path = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2018_2024_v2.tif"

# 因子映射
factor_map = np.array([np.nan, 1.00, 1.00, 0.58, 1.00, 0.50, 0.80, 0.00])


def load_and_crop_raster(raster_path, shp_gdf):
    """读取并根据 Shapefile 裁剪 Raster"""
    with rasterio.open(raster_path) as src:
        if str(src.crs) != str(shp_gdf.crs):
            shp_gdf = shp_gdf.to_crs(src.crs)
        out_image, out_transform = mask(src, shp_gdf.geometry, crop=True)
        data = out_image[0]
        nodata_val = src.nodata if src.nodata is not None else 0
        return data, out_transform, nodata_val, src.crs


def calculate_accurate_area_masked(data, transform, mask_array, crs):
    """计算面积 (用于统计表)"""
    res_x = transform[0]
    res_y = -transform[4]
    height, width = data.shape
    is_wgs84 = (crs.to_epsg() == 4326)

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


def main_workflow():
    print("1. 正在读取 Shapefile...")
    try:
        shp_gdf = gpd.read_file(path_shp)
    except Exception as e:
        print(f"错误: {e}")
        return

    print("2. 正在加载和处理栅格数据...")
    data_start, transform, nodata_start, crs = load_and_crop_raster(path_start_year, shp_gdf)
    data_end, _, nodata_end, _ = load_and_crop_raster(path_end_year, shp_gdf)

    if data_start.shape != data_end.shape:
        print("错误：两年份数据形状不一致。")
        return

    # --- 数据清理与计算 ---
    valid_range_start = (data_start >= 1) & (data_start <= 7)
    valid_range_end = (data_end >= 1) & (data_end <= 7)
    final_valid_pixels = valid_range_start & valid_range_end

    idx_start = data_start.copy()
    idx_end = data_end.copy()
    idx_start[~final_valid_pixels] = 0
    idx_end[~final_valid_pixels] = 0

    f_start = factor_map[idx_start]
    f_end = factor_map[idx_end]

    # 计算变化率
    calc_mask = (f_start > 0) & (~np.isnan(f_start)) & (~np.isnan(f_end))
    soc_change = np.full(data_start.shape, np.nan)
    soc_change[calc_mask] = (f_end[calc_mask] - f_start[calc_mask]) / f_start[calc_mask] * 100

    # --- 分类逻辑 ---
    # Degraded: < -10%
    mask_degraded = (soc_change < -10) & calc_mask
    # Stable: -10% <= change <= 10%
    mask_stable = (soc_change >= -10) & (soc_change <= 10) & calc_mask
    # Improved: > 10%
    mask_improved = (soc_change > 10) & calc_mask

    # 剩余区域为 NoData

    # =================================================================================
    # 新增部分：生成并保存 SOC Status TIFF
    # =================================================================================
    print(f"3. 正在生成 SOC Status GeoTIFF: {os.path.basename(output_tif_path)}")

    # 创建输出矩阵，默认填充 0 (NoData)
    # 使用 uint8 节省空间 (0-255)
    soc_status_raster = np.zeros(data_start.shape, dtype=np.uint8)

    # 赋值
    # 1: Degraded
    soc_status_raster[mask_degraded] = 1
    # 2: Stable
    soc_status_raster[mask_stable] = 2
    # 3: Improved
    soc_status_raster[mask_improved] = 3

    # 确保输出目录存在
    output_dir = os.path.dirname(output_tif_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 保存文件
    meta = {
        'driver': 'GTiff',
        'height': soc_status_raster.shape[0],
        'width': soc_status_raster.shape[1],
        'count': 1,
        'dtype': 'uint8',
        'crs': crs,
        'transform': transform,
        'nodata': 0,
        'compress': 'lzw'  # 推荐使用压缩
    }

    with rasterio.open(output_tif_path, 'w', **meta) as dst:
        dst.write(soc_status_raster, 1)

    print(f"   已保存: {output_tif_path}")
    print(f"   像素值说明: 1=Degraded, 2=Stable, 3=Improved, 0=NoData")

    # =================================================================================
    # 继续原本的统计计算
    # =================================================================================
    print("4. 正在计算统计表...")
    inside_shp_mask = (data_start != nodata_start)
    mask_nodata_area = inside_shp_mask & (~mask_degraded) & (~mask_stable) & (~mask_improved)

    area_degraded = calculate_accurate_area_masked(data_start, transform, mask_degraded, crs)
    area_stable = calculate_accurate_area_masked(data_start, transform, mask_stable, crs)
    area_improved = calculate_accurate_area_masked(data_start, transform, mask_improved, crs)
    area_nodata = calculate_accurate_area_masked(data_start, transform, mask_nodata_area, crs)

    total_area = area_degraded + area_stable + area_improved + area_nodata

    df = pd.DataFrame({
        "Category": [
            "Land area with degraded SOC",
            "Land area with stable SOC",
            "Land area with improved SOC",
            "Land area with no SOC data"
        ],
        "Area (km²)": [area_degraded, area_stable, area_improved, area_nodata],
        "Percent (%)": [
            area_degraded / total_area * 100 if total_area else 0,
            area_stable / total_area * 100 if total_area else 0,
            area_improved / total_area * 100 if total_area else 0,
            area_nodata / total_area * 100 if total_area else 0
        ]
    })

    print("\n=== Reporting Period (2018-2024) SOC Estimates ===")
    print(df.round(2))
    #csv_path = r"D:\deep\SoilGrids_Data\appendix_table3\SOC_Reporting_Period_2018_2024.csv"
    #df.to_csv(csv_path, index=False)
    #print(f"统计表已保存至: {csv_path}")


if __name__ == "__main__":
     main_workflow()

# #========================================================================分别check三个图
# import rasterio
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.colors import ListedColormap
#
# # 设置你的文件路径
# # tif_path = r"D:\deep\SoilGrids_Data\appendix_table3\soc_status_2018_2024.tif" # SOC
# #tif_path = r"D:\deep\fromFeishu\hmq\Land_Cover_Degradation\status_2018_2024_deg_stable_imp.tif" # land_cover
# tif_path = r"D:\deep\fromFeishu\qyx\land_productivity_20260122\combined_degradation_2024_table45.tif" # productivity
#
# def check_soc_tif(file_path):
#     print(f"正在检查文件: {file_path}")
#
#     try:
#         with rasterio.open(file_path) as src:
#             # 读取第一个波段
#             data = src.read(1)
#             profile = src.profile
#
#             print("\n--- 1. 文件基础信息 ---")
#             print(f"尺寸: {src.width} x {src.height}")
#             print(f"坐标系: {src.crs}")
#             print(f"数据类型: {profile['dtype']}")
#             print(f"NoData值: {src.nodata}")
#
#             # 统计唯一值
#             unique_vals, counts = np.unique(data, return_counts=True)
#
#             print("\n--- 2. 像素值统计 (验证关键) ---")
#             # 创建一个字典方便查看
#             stats = dict(zip(unique_vals, counts))
#
#             # # 定义标准对照  soc
#             # labels = {
#             #     0: "No Data (背景/水体)",
#             #     1: "Degraded (退化)",
#             #     2: "Stable (稳定)",
#             #     3: "Improved (改善)"
#             # }
#
#             # 定义标准对照  land_cover
#             labels = {
#                 0: "Stable",
#                 1: "Improved ",
#                 2: "Degraded",
#                 3: "no label"
#             }
#
#             total_valid_pixels = 0
#             for val in unique_vals:
#                 label = labels.get(val, "未知值")
#                 count = stats[val]
#                 print(f"值 {val} [{label}]: {count} 个像素")
#
#                 if val in [1, 2, 3]:
#                     total_valid_pixels += count
#
#             print("-" * 30)
#             if total_valid_pixels > 0:
#                 print("✅ 检测结果: 成功！文件中包含有效的分类数据 (1, 2, 3)。")
#             else:
#                 print("❌ 检测结果: 失败。文件中只有 NoData (0)，可能是计算过程被全部过滤了。")
#
#             # --- 3. 可视化预览 ---
#             # 为了让你看清楚，我们手动上色
#             # 0: 白色(透明), 1: 红色, 2: 浅黄色, 3: 绿色
#             print("\n--- 3. 正在生成预览图... ---")
#
#             # 创建自定义色表
#             # 颜色顺序对应值: 0, 1, 2, 3
#             colors = ['white', 'red', 'gold', 'green']
#             cmap = ListedColormap(colors)
#
#             # 限制数据范围在0-3之间以防万一
#             plot_data = data.copy()
#             plot_data[plot_data > 3] = 0
#
#             plt.figure(figsize=(10, 8))
#             im = plt.imshow(plot_data, cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
#
#             # 创建图例
#             cbar = plt.colorbar(im, ticks=[0.375, 1.125, 1.875, 2.625])
#             cbar.ax.set_yticklabels(['No Data', 'Degraded (Red)', 'Stable (Yellow)', 'Improved (Green)'])
#
#             plt.title("SOC Status Verification Check")
#             plt.axis('off')
#             plt.show()
#
#     except Exception as e:
#         print(f"读取文件出错: {e}")
#
#
# if __name__ == "__main__":
#     check_soc_tif(tif_path)


#===============================================================================
"""
输出为：
[!] 警告: 数据是 WGS84 (经纬度) 坐标系。
    不能直接相乘计算面积。需要使用纬度加权算法。
"""
# import rasterio
#
# # 你的文件路径
# path_2000 = r"D:\deep\fromFeishu\Land_Cover_Each_Year\C3S_2000_reclass_300m.tif"
#
# try:
#     with rasterio.open(path_2000) as src:
#         print("=== TIF 文件元数据检查 ===")
#         print(f"CRS (坐标系): {src.crs}")
#         print(f"Units (单位): {src.crs.linear_units if src.crs else 'Unknown'}")
#         print(f"Transform (仿射变换): \n{src.transform}")
#         print(f"Resolution (分辨率): {src.res}")
#
#         # 判断逻辑
#         if src.crs and src.crs.to_epsg() == 4326:
#             print("\n[!] 警告: 数据是 WGS84 (经纬度) 坐标系。")
#             print("    不能直接相乘计算面积。需要使用纬度加权算法。")
#         elif src.crs and src.crs.linear_units == 'metres':
#             print("\n[OK] 数据是投影坐标系 (米)。可以直接计算面积。")
#         else:
#             print(f"\n[?] 未知或特殊的单位: {src.crs.linear_units}")
#
# except Exception as e:
#     print(f"无法读取文件: {e}")