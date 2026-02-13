# appendix tabel1 0-30 soc
"""
- 使用 Google Earth Engine (GEE) Python API 进行认证与初始化，并强制绑定到指定 GEE 项目 'soc-uae'；
- 以 “United Arab Emirates（阿联酋）” 为研究区 ROI；
- 从 SoilGrids（ISRIC）资产中批量导出 3 个土壤变量在 3 个深度层的均值栅格（共 9 张 GeoTIFF）到 Google Drive。
"""
import ee
import os
import shutil

# ================= 1. 自动清理旧凭证 (防弹操作) =================
# 获取凭证文件路径
credentials_path = os.path.expanduser("C:\\Users\86180/.config/earthengine/credentials")
if os.path.exists(credentials_path):
    try:
        os.remove(credentials_path)
        print(f"已自动删除旧凭证文件: {credentials_path}")
        print("请在弹出的浏览器中重新登录。")
    except Exception as e:
        print(f"无法删除凭证文件，请手动删除: {e}")

# ================= 2. 初始化与强制绑定 =================
try:
    # 强制重新验证
    ee.Authenticate(force=True)

    # 初始化项目
    # http_transport=None 有时能解决SSL库冲突
    ee.Initialize(project='soc-uae', http_transport=None)

    # 【核心补丁】强制检查并修正内部项目设置
    # 这步是为了防止 GEE 库在导出时“忘记”项目 ID
    try:
        # 尝试使用新版 API 设置
        if hasattr(ee.data, 'setProject'):
            ee.data.setProject('soc-uae')
        # 暴力检查内部变量 (针对旧版库的 Bug)
        if hasattr(ee.data, '_cloud_api_resource'):
            current_resource = ee.data._cloud_api_resource
            if '517222506229' in current_resource or current_resource is None:
                print("检测到项目 ID 未正确绑定，正在强制修正...")
                ee.data._cloud_api_resource = 'projects/soc-uae'
    except Exception as patch_error:
        print(f"应用补丁时出现非致命警告: {patch_error}")

    print("初始化成功！当前绑定的项目: soc-uae")

except Exception as e:
    print(f"初始化严重失败: {e}")
    print("建议尝试：pip install --upgrade earthengine-api")
    exit(1)


def export_soilgrids_layers():
    # ================= 配置区域 =================

    # 1. 定义研究区 (ROI)
    countries = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
    roi = countries.filter(ee.Filter.eq('country_na', 'United Arab Emirates')).geometry()

    # 2. 导出参数
    # scale=250 是 SoilGrids 的原生分辨率
    SCALE = 250
    FOLDER_NAME = "SoilGrids_For_SOC"

    # ===========================================

    # 定义我们要下载的 9 个组合
    variables = ['soc', 'bdod', 'cfvo']
    depths = ['0-5cm', '5-15cm', '15-30cm']

    print(f"开始提交 9 个导出任务到 Google Drive (文件夹: {FOLDER_NAME})...")

    for var in variables:
        image_id = f"projects/soilgrids-isric/{var}_mean"
        for depth in depths:

            # 2. 选择对应的波段
            # 波段名称规律: soc_0-5cm_mean
            band_name = f"{var}_{depth}_mean"

            # 3. 获取图像并裁剪
            img = ee.Image(image_id).select(band_name).clip(roi)

            # 4. 构建输出文件名 (必须与之前的处理代码匹配)
            # 例如: soc_0-5cm_mean.tif
            file_name_prefix = f"{var}_{depth}_mean"

            # 5. 创建导出任务
            task = ee.batch.Export.image.toDrive(
                image=img,
                description=f"Export_{file_name_prefix}",
                folder=FOLDER_NAME,
                fileNamePrefix=file_name_prefix,
                region=roi,
                scale=SCALE,
                crs='EPSG:4326',  # WGS84 坐标系，与 C3S 兼容
                maxPixels=1e13,  # 允许导出大文件
                fileFormat='GeoTIFF'
            )

            # 6. 启动任务
            task.start()
            print(f" -> 已提交任务: {file_name_prefix}")



if __name__ == "__main__":
    export_soilgrids_layers()