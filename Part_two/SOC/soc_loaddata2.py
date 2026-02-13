import ee
import os

# ================= 1. 初始化 =================
try:
    ee.Initialize(project='soc-uae', http_transport=None)
    print("GEE 初始化成功！")
except Exception as e:
    print(f"初始化失败: {e}")
    # ee.Authenticate(force=True)

# ================= 2. 定义研究区域 =================
print("获取阿联酋边界...")
uae_roi = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017") \
    .filter(ee.Filter.eq("country_na", "United Arab Emirates"))

# ================= 3. 获取物理属性 (Sand & Clay) =================
print("正在计算沙和粘土的 0-30cm 平均值...")


def get_soil_prop_safe(asset_path, var_name):
    """
    获取土壤属性，计算 0-30cm 加权平均，并强制转换为 Float 以避免类型冲突。
    """
    try:
        source = ee.Image(asset_path)

        # 尝试获取分层数据
        # 权重: 0-5cm(5), 5-15cm(10), 15-30cm(15) -> 总和 30
        l1 = source.select(f'{var_name}_0-5cm_mean')
        l2 = source.select(f'{var_name}_5-15cm_mean')
        l3 = source.select(f'{var_name}_15-30cm_mean')

        # 计算加权平均
        weighted_mean = l1.multiply(5).add(l2.multiply(10)).add(l3.multiply(15)).divide(30)

        # 【关键修正】强制转换为 Float，并裁剪
        return weighted_mean.toFloat().rename(f'{var_name}_0-30cm_mean').clip(uae_roi)

    except Exception as e:
        print(f"警告: 无法获取 {var_name} 的分层数据 ({e})")
        return None


# 获取数据
sand_img = get_soil_prop_safe("projects/soilgrids-isric/sand_mean", "sand")
clay_img = get_soil_prop_safe("projects/soilgrids-isric/clay_mean", "clay")

# ================= 4. 合并与导出 =================
if sand_img and clay_img:
    print("合并波段...")
    # 将两个浮点数图像合并
    final_image = sand_img.addBands(clay_img)

    print("创建导出任务...")
    task = ee.batch.Export.image.toDrive(
        image=final_image,
        description='UAE_SoilGrids_Sand_Clay_0-30cm',
        folder='GEE_Soil_Downloads',
        region=uae_roi.geometry(),
        scale=250,
        crs='EPSG:4326',
        maxPixels=1e13,
        fileFormat='GeoTIFF'
    )

    task.start()
    print(f"任务已启动! ID: {task.id}")
    print("包含波段: sand_0-30cm_mean, clay_0-30cm_mean")
    print("请前往 Google Drive 查看结果。")
else:
    print("错误: 未能成功加载沙或粘土数据，无法导出。")