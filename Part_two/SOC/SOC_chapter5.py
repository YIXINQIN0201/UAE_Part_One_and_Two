"""
基于 SoilGrids 全球土壤数据库，精确计算并生成 0-30cm 深度范围内“土壤有机碳储量（SOC Stock）”的基准地图，并按年度统计不同土地利用类型下的碳储量变化。
$$SOC\ Stock\ (t/ha) = SOC \times BDOD \times Depth \times (1 - CFVO) \times 0.1$$

SOC (土壤有机碳含量): 每克土壤中有多少碳 ($dg/kg$)。
BDOD (容重/体密度): 土壤的紧实程度 ($cg/cm^3$)。
CFVO (粗碎屑/砾石含量): 扣除掉不含碳的石头部分。
Depth (深度): 代码分三层计算（0-5cm, 5-15cm, 15-30cm）并求和，得到 0-30cm 的总储量。
"""


import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from tqdm import tqdm

# ================= 配置区域 =================

# 1. 您的本地数据路径
# SoilGrids 9个 TIFF 文件所在的文件夹
SOIL_DIR = r"D:\deep\SoilGrids_Data"

# 土地利用数据所在的文件夹 (2000-2024年的 tif)
LC_DIR = r"D:\deep\fromFeishu\Land_Cover_Each_Year"

# 结果输出位置
OUTPUT_CSV = r"D:\deep\soc_stock_results_2000_2024.csv"

# 【新增】输出 SOC 储量栅格图的路径
OUTPUT_IMAGE_PATH = r"D:\deep\SOC_Stock_0-30cm_UAE_Baseline.tif"

# 2. 土地利用分类字典 (已更新为您确认的 1-7 对应关系)
CLASS_MAPPING = {
    1: "Forest Land",
    2: "Grassland",
    3: "Cropland",
    4: "Wetland",
    5: "Artificial",
    6: "Other Land",
    7: "Water bodies"
}


# ===========================================

def get_resampled_layer(path, reference_profile, reference_transform, reference_shape):
    """
    读取并重采样 SoilGrids 图层以匹配 Land Cover 的分辨率和范围
    """
    with rasterio.open(path) as src:
        destination = np.zeros(reference_shape, dtype=np.float32)

        reproject(
            source=rasterio.band(src, 1),
            destination=destination,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=reference_transform,
            dst_crs=reference_profile['crs'],
            resampling=Resampling.bilinear  # SOC 是连续变量，必须用双线性插值
        )
        return destination


def calculate_soc_stock_map(reference_tif_path):
    """
    计算 0-30cm SOC Stock (t/ha) 并对齐到 Land Cover 网格
    """
    print(">>> [Step 1/2] 正在计算基准 SOC 储量分布图 (0-30cm)...")

    # 读取参考影像 (Land Cover) 的元数据
    with rasterio.open(reference_tif_path) as ref:
        ref_profile = ref.profile
        ref_transform = ref.transform
        ref_shape = (ref.height, ref.width)

    # 定义分层信息
    layers_info = [
        {'label': '0-5cm', 'depth': 5},
        {'label': '5-15cm', 'depth': 10},
        {'label': '15-30cm', 'depth': 15}
    ]

    # 初始化总储量数组
    total_soc_stock = np.zeros(ref_shape, dtype=np.float32)

    for layer in layers_info:
        label = layer['label']
        depth = layer['depth']

        # 拼接文件路径
        f_soc = os.path.join(SOIL_DIR, f"soc_{label}_mean.tif")
        f_bd = os.path.join(SOIL_DIR, f"bdod_{label}_mean.tif")
        f_cf = os.path.join(SOIL_DIR, f"cfvo_{label}_mean.tif")

        # 检查文件是否存在
        for f in [f_soc, f_bd, f_cf]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"缺失文件: {f}")

        # 读取并重采样
        soc_raw = get_resampled_layer(f_soc, ref_profile, ref_transform, ref_shape)
        bd_raw = get_resampled_layer(f_bd, ref_profile, ref_transform, ref_shape)
        cf_raw = get_resampled_layer(f_cf, ref_profile, ref_transform, ref_shape)

        # --- 核心计算公式 ---
        # Stock (t/ha) = SOC(g/kg) * BD(g/cm3) * Depth(cm) * (1 - Gravel%) * 0.1

        # 1. 单位转换
        soc_g_kg = soc_raw * 0.1  # dg/kg -> g/kg
        bd_g_cm3 = bd_raw * 0.01  # cg/cm3 -> g/cm3
        cf_frac = cf_raw * 0.001  # permille -> fraction (0-1)
        cf_frac = np.clip(cf_frac, 0, 1)  # 修正异常值(>1的情况)

        # 2. 计算该层储量
        layer_stock = soc_g_kg * bd_g_cm3 * depth * (1 - cf_frac) * 0.1

        # 3. 累加到总储量
        total_soc_stock += layer_stock
        print(f"    -> 已合并图层: {label}")

    return total_soc_stock


def save_soc_map_as_tif(soc_map_array, template_path, output_path):
    """【新增函数】将计算好的 numpy 数组保存为 GeoTIFF 文件"""
    print(f"\n>>> [Step 2/3] 正在保存 SOC 储量栅格图至: {output_path} ...")

    # 读取模板文件的元数据以获取空间参考信息
    with rasterio.open(template_path) as src:
        profile = src.profile

    # 更新元数据以匹配我们的 SOC 数据类型 (float32)
    profile.update(
        dtype=rasterio.float32,
        count=1,
        nodata=np.nan,  # 设置无效值为 NaN
        compress='lzw'
    )

    # 写入文件
    with rasterio.open(output_path, 'w', **profile) as dst:
        # 将可能存在的非常小的负数（计算误差）修正为0，并处理NaN
        soc_map_array[soc_map_array < 0] = 0
        # 将 numpy 中的 nan 填充为 profile 中定义的 nodata 值
        soc_map_array = np.nan_to_num(soc_map_array, nan=profile['nodata'])
        dst.write(soc_map_array, 1)
    print(">>> 保存完成！")


def main():
    years = range(2000, 2025)  # 2000-2024

    # 1. 获取参考模板 (使用2000年的文件) 来确定网格大小
    template_path = os.path.join(LC_DIR, "C3S_2000_reclass_300m.tif")
    if not os.path.exists(template_path):
        print(f"错误: 找不到模板文件 {template_path}")
        return

    # 2. 生成 SOC 分布图 (这一步最耗时)
    soc_stock_map = calculate_soc_stock_map(template_path)

    # 【新增】3. 将计算结果保存为 TIFF 文件
    #save_soc_map_as_tif(soc_stock_map, template_path, OUTPUT_IMAGE_PATH)

    # 3. 逐年统计
    print("\n>>> [Step 2/2] 开始逐年统计 (2000-2024)...")
    all_results = []

    for year in tqdm(years, desc="Progress"):
        tif_name = f"C3S_{year}_reclass_300m.tif"
        tif_path = os.path.join(LC_DIR, tif_name)

        if not os.path.exists(tif_path):
            continue

        with rasterio.open(tif_path) as src:
            lc_data = src.read(1)

            # 尺寸校验
            if lc_data.shape != soc_stock_map.shape:
                print(f"警告: {year}年数据尺寸不符，跳过。")
                continue

            # 遍历 1-7 类
            for pixel_val, class_name in CLASS_MAPPING.items():

                # 创建掩膜 (只选当前类别的像素)
                mask = (lc_data == pixel_val)

                if np.any(mask):
                    # 提取对应位置的 SOC 值
                    soc_values = soc_stock_map[mask]

                    # 去除 SOC 无效值 (NaN)
                    valid_soc = soc_values[~np.isnan(soc_values)]

                    if valid_soc.size > 0:
                        mean_val = np.mean(valid_soc)
                        pixel_count = valid_soc.size

                        all_results.append({
                            "Year": year,
                            "Land Cover Class": class_name,
                            "Mean SOC Stock (t/ha)": round(mean_val, 4),
                            "Pixel Count": pixel_count
                        })
                    else:
                        # 有地类像素，但 SOC 数据缺失
                        all_results.append({
                            "Year": year,
                            "Land Cover Class": class_name,
                            "Mean SOC Stock (t/ha)": 0.0,
                            "Pixel Count": 0
                        })

    # 4. 导出 CSV
    if all_results:
        df = pd.DataFrame(all_results)

        # 为了让 CSV 按照 Forest, Grassland... 的顺序排列，我们创建一个辅助排序列
        # 将 Class Name 映射回 ID (1-7)
        name_to_id = {v: k for k, v in CLASS_MAPPING.items()}
        df['Class_ID'] = df['Land Cover Class'].map(name_to_id)

        # 排序: 先按年份，再按 Class ID
        df = df.sort_values(by=['Year', 'Class_ID'])

        # 删除辅助列
        df = df.drop(columns=['Class_ID'])

        # 保存
        #df.to_csv(OUTPUT_CSV, index=False, encoding='utf-8-sig')
        print(f"\n成功！结果已保存至: {OUTPUT_CSV}")
    else:
        print("未生成任何结果，请检查数据路径。")


if __name__ == "__main__":
    main()