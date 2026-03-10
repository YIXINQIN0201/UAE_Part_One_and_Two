#计算低生产力像元的掩膜
from __future__ import annotations
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from pathlib import Path

from dataclasses import dataclass
from rasterio.crs import CRS
from pathlib import Path   
import numpy as np
from typing import Any,Mapping,Union,Optional
from dataclasses import field
import geopandas as gpd
from rasterio.features import geometry_mask
import pandas as pd
from tqdm import tqdm
import math
import rasterio
from rasterio.windows import Window

from sympy import HadamardPower

# 可选：卫星底图（需要联网）
USE_BASEMAP = True
try:
    import contextily as cx  # pip install contextily
except Exception:
    USE_BASEMAP = False

@dataclass(frozen=False)
class MaskArgs:
    mask: bool = False
    slope_threshold: Union[float, list[float]] = 0.025


@dataclass(frozen=False)
class NppArgs:
    npp_path: Path
    height: int
    width: int
    nodata: int
    npp_profile: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=False)
class BasicArgs:
    # ---- required (无默认值的必须全放前面) ----
    metrics: str
    name_field: str
    dst_crs: CRS
    dst_res: int
    reporting_year: int
    out_dir: Path
    shp_path: Path
    lecu_path:Path = Path("../data/LCEU/lecu.tif")

    # ---- optional ----
    z_edges: list[float] = field(default_factory=list)
    mask_low_productivity: bool = False

    # ---- nested optional args ----
    # NppArgs 不能 default_factory=NppArgs（因为它需要必填参数）
    npp_args: Optional[NppArgs] = None

    trend_args: MaskArgs = field(default_factory=MaskArgs)
    state_args: MaskArgs = field(default_factory=MaskArgs)

    performance_threshold: Union[float, list[float]] = field(
        default_factory=lambda: [0.5, 1.0]
    )
def initalize_args(metrics:str,reporting_year:int) -> BasicArgs:
    """初始化评估参数配置。

    根据给定的指标类型和报告年份，创建并配置 BasicArgs 对象，
    包含评估所需的所有默认参数设置。

    Args:
        metrics: 评估指标类型，可以是 "trend"、"state" 或 "performance"。
        reporting_year: 报告年份，用于确定输入数据和输出目录路径。

    Returns:
        BasicArgs: 包含完整评估参数配置的数据类实例。

    Note:
        此函数设置了默认的参数值，包括：
        - 名称字段：NAME_1
        - 目标坐标系：EPSG:6933
        - 目标分辨率：250米
        - Z值分箱边界：[-inf, -1.96, -0.5, 0.5, 1.96, inf]
        - 趋势和状态参数的默认掩膜设置
        - NPP数据路径（基于报告年份计算）
        - 输出目录路径

    Example:
        >>> args = initalize_args("trend", 2023)
        >>> print(args.metrics)
        'trend'
    """
    # basic arguments
    name_field = "NAME_1"
    dst_crs = "EPSG:6933" 
    dst_res = 250
    z_edges = [-np.inf, -1.96, -0.5, 0.5, 1.96, np.inf]
    trend_args = MaskArgs(mask = False,slope_threshold = 0.025)
    state_args = MaskArgs(mask = False,slope_threshold = 0.025)
    npp_path = Path(f"../data/NPP_PROXY/annual_auc_toBands_{reporting_year-15}_{reporting_year}.tif")  # 改成你的tif
    npp_args = prepare_profile(npp_path)

    out_dir = Path(f"../output/{metrics}/{reporting_year}/")  # 想输出到哪里
    out_dir.mkdir(parents=True, exist_ok=True)

    #input path
    shp_path  = Path("../abu_dhabi_all/gadm41_ARE_1.shp")   # 7个酋长国边界
    return BasicArgs(
    metrics=metrics,
    name_field=name_field,
    dst_crs=dst_crs,
    dst_res=dst_res,
    z_edges=z_edges,
    reporting_year=reporting_year,
    out_dir=out_dir,
    shp_path=shp_path,
    mask_low_productivity=False,
    npp_args=npp_args,
    trend_args=trend_args,
    state_args=state_args
)
#-------------------------------------------------------------------------------
def prepare_profile(npp_path: Path) -> NppArgs:
    """从NPP栅格文件中提取元数据并创建NppArgs对象。

    读取NPP栅格文件的元数据信息，包括高度、宽度、nodata值和profile配置，
    用于后续处理时的参数传递。

    Args:
        npp_path: NPP栅格文件的路径。

    Returns:
        NppArgs: 包含NPP文件元数据的数据类实例，包括：
            - npp_path: 文件路径
            - height: 栅格高度（行数）
            - width: 栅格宽度（列数）
            - nodata: 无数据值
            - npp_profile: 完整的栅格profile配置字典

    Example:
        >>> npp_path = Path("../data/NPP_PROXY/annual_auc_toBands_2008_2023.tif")
        >>> npp_args = prepare_profile(npp_path)
        >>> print(npp_args.height, npp_args.width)
    """
    with rasterio.open(npp_path) as src:
        height = src.height
        width = src.width
        nodata = src.nodata
        npp_args = NppArgs(npp_path=npp_path,height = height, width = width, nodata = nodata,npp_profile=src.profile.copy(),
        )  # 读到原始 profile
    return npp_args
#-------------------------------------------------------------------------------
class PlotArgs:
    USE_BASEMAP: bool = True
    OUTDIR: Path = Path("../output/plots/emirate_local_figs")
    RADIUS_M: int = 30000
# 最基本的参数设定--------------------------------


def mask_low_productivity(y_prod: np.ndarray) -> np.ndarray:
    """生成低生产力像元的掩膜。

    基于生产力数据的时间序列，识别并掩膜低生产力像元。
    使用均值、90分位数和有效年份数作为筛选条件。

    Args:
        y_prod: 生产力数据数组，形状为 (T, H, W)，其中T为时间维度（年份），
                H为高度，W为宽度。数据类型为float32，NaN表示无效值。

    Returns:
        np.ndarray: 布尔掩膜数组，形状为 (H, W)。True表示有效的高生产力像元，
                    False表示低生产力像元。掩膜标准包括：
                    - 有效年份数 >= 12
                    - 平均生产力 >= 30分位数阈值
                    - 90分位数 >= 30分位数阈值

    Note:
        - 使用30%分位数作为阈值，去除底部几乎没生产力的像元
        - 至少需要12个有效年份才认为像元可靠
        - 打印掩膜百分比和阈值信息用于调试

    Example:
        >>> y_prod = np.random.rand(16, 100, 100).astype(np.float32)
        >>> mask = mask_low_productivity(y_prod)
        >>> print(f"有效像元比例: {mask.mean() * 100}%")
    """
    # y_prod: (T, H, W) float32, NaN=invalid
    valid_prod = np.isfinite(y_prod)
    N_prod = valid_prod.sum(axis=0)

    mean_prod = np.nanmean(y_prod, axis=0)
    p90_prod  = np.nanpercentile(y_prod, 90, axis=0)

    # ---- 选择阈值（推荐：分位数阈值）----
    valid_mean = mean_prod[np.isfinite(mean_prod)]
    valid_p90  = p90_prod[np.isfinite(p90_prod)]

    # 去掉“几乎没生产力”的底部像元：你可以先用 30% 试试
    mean_thr = float(np.nanpercentile(valid_mean, 30))
    p90_thr  = float(np.nanpercentile(valid_p90, 30))

    # 至少多少个有效年份才算可靠
    valid_years_min = 12  # 2000-2015 共16年，建议 10~14 之间

    prod_mask = (N_prod >= valid_years_min) & (mean_prod >= mean_thr) & (p90_prod >= p90_thr)

    print("productive mask %:", float(np.mean(prod_mask) * 100))
    print("mean_thr:", mean_thr, "p90_thr:", p90_thr, "valid_years_min:", valid_years_min)
    return prod_mask

def read_block_to_float(src, win, nodata):
    """从栅格文件中读取指定窗口的数据块并转换为浮点数。

    读取栅格数据的指定窗口区域，将nodata值转换为NaN，
    以便进行数值计算。

    Args:
        src: rasterio DatasetReader对象，打开的栅格数据源。
        win: rasterio.windows.Window对象，指定要读取的窗口区域。
        nodata: 无数据值。如果为None或NaN，则不进行转换。

    Returns:
        np.ndarray: 读取的数据数组，形状为 (T, h, w)，其中T为波段数，
                    h和w为窗口的高度和宽度。数据类型为float32。
                    nodata值已转换为np.nan。

    Example:
        >>> with rasterio.open("data.tif") as src:
        ...     win = Window(0, 0, 512, 512)
        ...     data = read_block_to_float(src, win, nodata=-9999)
        ...     print(data.shape, data.dtype)
    """
    data = src.read(window=win).astype(np.float32)  # (T,h,w)
    if nodata is not None and not np.isnan(nodata):
        # 把nodata 变成nan
        data = np.where(data == nodata, np.nan, data)
    return data

#projection
def reproject_to_equal_area(
    src_path: Path,
    dst_path: Path,
    dst_crs: str,
    dst_res: float,
    resampling: Resampling = Resampling.bilinear,
    *,
    force_float32: bool = True,
    dst_nodata: float = np.nan,
    copy_band_descriptions: bool = True,
    block_size: int = 512,
    compress: str = "lzw",
) -> None:
    """将栅格（单波段或多波段）重投影到等面积坐标系并写入文件。

    统一版本：兼容你原来的 `reproject_to_equal_area` 与 `reproject_to_equal_area_for_npp`。
    - 自动保留所有波段（count=src.count）
    - 输出 float32，nodata 统一为 NaN（默认）
    - 将源 nodata 值统一转换为 NaN，再进行重投影
    - 可选复制波段描述（band descriptions）

    Args:
        src_path: 源栅格路径（单/多波段均可）。
        dst_path: 输出栅格路径。
        dst_crs: 目标 CRS，如 "EPSG:6933"。
        dst_res: 目标分辨率（单位与 dst_crs 一致，通常为米）。
        resampling: 重采样方法（默认 bilinear）。
        force_float32: 是否强制输出 float32（默认 True）。
        dst_nodata: 输出 nodata 值（默认 NaN）。
        copy_band_descriptions: 是否复制 src.descriptions 到输出。
        block_size: 分块大小（默认 512）。
        compress: 压缩方式（默认 LZW）。

    Returns:
        None（直接写文件）。

    Notes:
        - 输出：tiled + blockxsize/blockysize + LZW（可调）
        - 如果源数据有 nodata（如 -9999、0 等），会在重投影前被替换为 NaN。
        - 若源 nodata 本身就是 NaN 或 None，也能正常工作。
    """
    src_path = Path(src_path)
    dst_path = Path(dst_path)

    with rasterio.open(src_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, resolution=dst_res
        )

        out_dtype = "float32" if force_float32 else src.dtypes[0]

        profile = src.profile.copy()
        profile.update(
            crs=dst_crs,
            transform=transform,
            width=width,
            height=height,
            compress=compress,
            tiled=True,
            blockxsize=block_size,
            blockysize=block_size,
            dtype=out_dtype,
            nodata=dst_nodata,
            count=src.count,
        )

        # 处理源 nodata：尽量在重投影前把 nodata 统一成 NaN
        src_nodata = src.nodata
        src_has_finite_nodata = src_nodata is not None and np.isfinite(src_nodata)

        with rasterio.open(dst_path, "w", **profile) as dst:
            for b in range(1, src.count + 1):
                src_data = src.read(b)

                # 强制 float32 以支持 NaN nodata
                if force_float32:
                    src_data = src_data.astype(np.float32, copy=False)

                # src nodata -> NaN（只在 nodata 是有限数值时做替换）
                if src_has_finite_nodata:
                    src_data = np.where(src_data == src_nodata, np.nan, src_data).astype(
                        np.float32 if force_float32 else src_data.dtype,
                        copy=False,
                    )

                dst_data = np.full((height, width), dst_nodata, dtype=np.float32 if force_float32 else src_data.dtype)

                reproject(
                    source=src_data,
                    destination=dst_data,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=resampling,
                    # 关键：我们已经把 src_data 的 nodata 转成 NaN，因此这里用 NaN
                    src_nodata=np.nan if force_float32 else src_nodata,
                    dst_nodata=dst_nodata,
                )

                dst.write(dst_data, b)

            if copy_band_descriptions and src.descriptions:
                dst.descriptions = src.descriptions
#*********************************************************************************


def pixel_area_km2(src: rasterio.io.DatasetReader) -> float:
    """计算栅格像元的面积（平方公里）。

    基于栅格的仿射变换参数计算单个像元的面积。
    仅适用于等面积投影的栅格数据。

    Args:
        src: rasterio DatasetReader对象，已打开的栅格数据源。

    Returns:
        float: 单个像元的面积（单位：平方公里）。

    Raises:
        ValueError: 如果栅格使用地理坐标系（度为单位）而不是投影坐标系，
                   会抛出错误，提示需要先重投影到等面积坐标系。

    Note:
        - transform.a 表示像元宽度，transform.e 表示像元高度（通常为负值）
        - 仅当像元尺寸大于0.1时认为使用投影坐标系（米）
        - 计算结果 = |width * height| / 1e6 （转换为km²）

    Example:
        >>> with rasterio.open("data.tif") as src:
        ...     area = pixel_area_km2(src)
        ...     print(f"像元面积: {area} km²")
    """
    # transform.a = pixel width, transform.e = pixel height (通常为负)
    px_w = src.transform.a
    px_h = src.transform.e
    # 粗略判断是否米单位：绝对值像素尺寸一般在几十到几千之间（比如 250m）
    if abs(px_w) < 0.1 or abs(px_h) < 0.1:
        raise ValueError(
            "Raster appears to be in degrees (geographic CRS). "
            "For accurate area, reproject to an equal-area projection first."
        )
    return abs(px_w * px_h) / 1e6  # km^2

def mask_low_productivity(y_prod: np.ndarray) -> np.ndarray:
# y_prod: (T, H, W) float32, NaN=invalid
    valid_prod = np.isfinite(y_prod)
    N_prod = valid_prod.sum(axis=0)

    mean_prod = np.nanmean(y_prod, axis=0)
    p90_prod  = np.nanpercentile(y_prod, 90, axis=0)

    # ---- 选择阈值（推荐：分位数阈值）----
    valid_mean = mean_prod[np.isfinite(mean_prod)]
    valid_p90  = p90_prod[np.isfinite(p90_prod)]

    # 去掉“几乎没生产力”的底部像元：你可以先用 30% 试试
    mean_thr = float(np.nanpercentile(valid_mean, 30))
    p90_thr  = float(np.nanpercentile(valid_p90, 30))

    # 至少多少个有效年份才算可靠
    valid_years_min = 12  # 2000-2015 共16年，建议 10~14 之间

    prod_mask = (N_prod >= valid_years_min) & (mean_prod >= mean_thr) & (p90_prod >= p90_thr)

    print("productive mask %:", float(np.mean(prod_mask) * 100))
    print("mean_thr:", mean_thr, "p90_thr:", p90_thr, "valid_years_min:", valid_years_min)
    return prod_mask

def mask_small_diff(s: np.ndarray,threshold: float = 0.0025) -> np.ndarray:
    """掩膜差异值较小的像元。

    基于阈值过滤掉变化幅度较小的像元，通常用于去除
    统计上不显著的差异。

    Args:
        s: 差异值数组，可以是斜率、差值等任何数值数组。
        threshold: 差异阈值，默认0.0025。只有大于此值的像元才会被保留。

    Returns:
        np.ndarray: 布尔掩膜数组，与输入数组形状相同。
                   True表示差异显著（大于阈值）的像元。

    Example:
        >>> slope = np.array([0.001, 0.005, 0.01])
        >>> mask = mask_small_diff(slope, threshold=0.0025)
        >>> print(mask)  # [False, True, True]
    """
    return s>threshold

def mask_small_slope(s: np.ndarray, threshold: float) -> np.ndarray:
    """掩膜斜率绝对值较小的像元。

    基于阈值过滤掉斜率绝对值较小的像元，用于去除
    变化幅度较小的趋势。

    Args:
        s: 斜率数组，可以是任意形状的numpy数组。
        threshold: 斜率阈值。绝对值小于此值的像元将被掩膜（返回0），
                  绝对值大于等于此值的像元将被保留（返回1）。

    Returns:
        np.ndarray: 掩膜数组，与输入数组形状相同，数据类型为float32。
                   0表示斜率绝对值小于阈值（被掩膜），
                   1表示斜率绝对值大于等于阈值（保留）。

    Example:
        >>> slope = np.array([-0.01, -0.001, 0.001, 0.01])
        >>> mask = mask_small_slope(slope, threshold=0.0025)
        >>> print(mask)  # [1.0, 0.0, 0.0, 1.0]
    """
    mask = np.abs(s) >= threshold
    return mask

def classify_bins(z: np.ndarray,Z_EDGES) -> np.ndarray:
    """将Z值数组按照指定的分箱边界进行分类。

    将连续的Z值（如z-score）转换为离散的分箱类别，
    用于后续的退化/改善土地分类。

    Args:
        z: Z值数组（如z-score），可以是任意形状的numpy数组。
           可能已经过生产力和小差异筛选，包含NaN值。
        Z_EDGES: 分箱边界列表，如 [-inf, -1.96, -0.5, 0.5, 1.96, inf]。
                 分箱数量 = len(Z_EDGES) - 1。

    Returns:
        np.ndarray: 分箱ID数组，数据类型为int16，与输入数组形状相同。
                   编码规则：
                   - 0: 掩膜/无效值（未使用的编码）
                   - 1..N: 分箱编号（N为分箱数量，由Z_EDGES决定）
                   - 99: 无效值（NaN）

    Note:
        - 使用np.digitize进行分箱，返回1到len(Z_EDGES)-1的值
        - 无效值（NaN）标记为99
        - 输入z可能已经过生产力和小差异筛选

    Example:
        >>> z = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        >>> edges = [-np.inf, -1.96, -0.5, 0.5, 1.96, np.inf]
        >>> bins = classify_bins(z, edges)
        >>> print(bins)  # [1, 2, 3, 4, 5]
    """
    bin_id = np.zeros(z.shape, dtype=np.int16)

    gate = np.isfinite(z) 

    # gate 未通过 -> 99
    bin_id[~gate] = 99

    # gate 通过 -> 按 z 分箱 1..5
    # np.digitize 返回 1..len(edges)-1
    idx = np.digitize(z[gate], Z_EDGES)  # 1..5
    bin_id[gate] = idx.astype(np.int16)

    return bin_id
#-------------------------------------------------------------------------------

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path


def setup_logger(
    name: str = "app",
    log_file: str | Path = "app.log",
    level: int = logging.INFO,
    fmt: str = "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt: str = "%Y-%m-%d %H:%M:%S",
    max_bytes: int = 10 * 1024 * 1024,  # 10MB
    backup_count: int = 5,
    console: bool = True,
) -> logging.Logger:
    """创建同时输出到控制台和滚动日志文件的日志记录器。

    配置一个支持文件滚动和双重输出（文件+控制台）的日志记录器。
    防止重复配置（如在notebook中多次调用）。

    Args:
        name: 日志记录器名称，默认为 "app"。
        log_file: 日志文件路径（字符串或Path对象），默认为 "app.log"。
        level: 日志级别，默认为logging.INFO。
        fmt: 日志格式字符串，默认为包含时间、级别、名称和消息的格式。
        datefmt: 日期时间格式字符串，默认为 "%Y-%m-%d %H:%M:%S"。
        max_bytes: 单个日志文件的最大字节数，超过此值会触发滚动，
                  默认为10MB。
        backup_count: 保留的备份日志文件数量，默认为5。
        console: 是否同时输出到控制台，默认为True。

    Returns:
        logging.Logger: 配置好的日志记录器对象。

    Note:
        - 使用RotatingFileHandler实现文件滚动
        - 防止重复添加handler（通过_check_configured标记）
        - 设置propagate=False避免重复打印

    Example:
        >>> logger = setup_logger("my_app", "app.log", level=logging.DEBUG)
        >>> logger.info("程序启动")
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # 防止重复添加 handler（比如 notebook / 多次调用）
    if getattr(logger, "_configured", False):
        return logger

    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)

    # 文件 handler（滚动）
    fh = RotatingFileHandler(
        filename=str(log_path),
        maxBytes=max_bytes,
        backupCount=backup_count,
        encoding="utf-8",
    )
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # 控制台 handler
    if console:
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    logger.propagate = False  # 不把日志再传给 root，避免重复打印
    logger._configured = True
    return logger

#-------------------------------------------------------------------------------
#trend functions
#-------------------------------------------------------------------------------
def theil_sen_slope_block(y: np.ndarray) -> np.ndarray:
    """计算Theil-Sen斜率估计值。

    使用Theil-Sen方法计算时间序列的斜率，该方法对异常值
    具有鲁棒性。计算所有时间点对之间的斜率，取中位数。

    Args:
        y: 时间序列数据数组，形状为 (T, h, w)，其中T为时间维度，
           h和w为空间维度。数据类型为float32，NaN表示无效值。

    Returns:
        np.ndarray: 斜率数组，形状为 (h, w)，数据类型为float32。
                   单位为单位/年（与输入数据的单位一致）。
                   只有所有时间点都有有效值的像元才会计算斜率。

    Note:
        - Theil-Sen方法：计算所有时间点对 (i, j) 的斜率 (yj-yi)/(j-i)，
          取所有斜率的中位数
        - 只对每年都有有效值的像元进行计算
        - 对异常值具有鲁棒性，是线性趋势的稳健估计

    Example:
        >>> y = np.random.rand(16, 100, 100).astype(np.float32)
        >>> slope = theil_sen_slope_block(y)
        >>> print(f"斜率范围: [{slope.min():.4f}, {slope.max():.4f}]")
    """
    T, h, w = y.shape
    y = y.astype(np.float32, copy=False)
    # 只有这些像素才算：每年都有值
    valid_all = np.all(np.isfinite(y), axis=0)  # (h,w)

    out = np.full((h, w), np.nan, dtype=np.float32)
    if not np.any(valid_all):
        return out

    slopes = []
    for i in range(T - 1):
        yi = y[i]
        for j in range(i + 1, T):
            yj = y[j]
            s = (yj - yi) / (j - i)
            slopes.append(s)

    S = np.stack(slopes, axis=0)  # (P,h,w)

    med = np.median(S[:, valid_all], axis=0).astype(np.float32)
    out[valid_all] = med
    return out

def mann_kendall_S(y: np.ndarray) -> np.ndarray:
    """计算Mann-Kendall检验的S统计量。

    计算Mann-Kendall趋势检验的S统计量，用于检测时间序列的
    单调趋势。S是所有时间点对之间符号函数的总和。

    Args:
        y: 时间序列数据数组，形状为 (T, h, w)，其中T为时间维度，
           h和w为空间维度。数据类型为float32，NaN表示无效值。

    Returns:
        np.ndarray: S统计量数组，形状为 (h, w)，数据类型为float32。
                   S = sum(sign(yj - yi)) for all i < j
                   - 正值表示上升趋势
                   - 负值表示下降趋势
                   - 只有所有时间点都有有效值的像元才会计算

    Note:
        - sign函数：+1（上升）、0（相等）、-1（下降）
        - 只对每年都有有效值的像元进行计算
        - 通常需要进一步转换为Z值进行显著性检验

    Example:
        >>> y = np.random.rand(16, 100, 100).astype(np.float32)
        >>> s = mann_kendall_S(y)
        >>> print(f"S统计量范围: [{s.min():.2f}, {s.max():.2f}]")
    """
    T, h, w = y.shape
    y = y.astype(np.float32, copy=False)
    # 只有这些像素才算：每年都有值
    valid_all = np.all(np.isfinite(y), axis=0)  # (h,w)

    out = np.full((h, w), np.nan, dtype=np.float32)
    if not np.any(valid_all):
        return out

    sign = []
    for i in range(T - 1):
        yi = y[i]
        for j in range(i + 1, T):
            yj = y[j]
            s = np.sign(yj - yi)
            sign.append(s)

    S = np.stack(sign, axis=0)  # (P,h,w)

    med = np.sum(S[:, valid_all], axis=0).astype(np.float32)
    out[valid_all] = med
    return out
    
def mann_kendall_z(y: np.ndarray) -> np.ndarray:
    """计算Mann-Kendall检验的Z值（标准化统计量）。

    将Mann-Kendall的S统计量转换为标准正态分布的Z值，
    用于检验趋势的统计显著性。

    Args:
        y: 时间序列数据数组，形状为 (T, h, w)，其中T为时间维度，
           h和w为空间维度。数据类型为float32，NaN表示无效值。

    Returns:
        np.ndarray: Z值数组，形状为 (h, w)，数据类型为float32。
                    Z = S / sqrt(Var(S))
                    - |Z| > 1.96 表示在95%置信水平下显著
                    - 正值表示上升趋势，负值表示下降趋势

    Raises:
        ValueError: 如果输入数组不是3维（T, H, W）形状。

    Note:
        - 考虑了数据中的并列值（ties）对方差的影响
        - 使用每像元独特的并列值修正方差计算
        - 只对每年都有有效值的像元进行计算
        - Var(S) = [T(T-1)(2T+5) - sum(t(t-1)(2t+5))] / 18
          （t为并列值数量）

    Example:
        >>> y = np.random.rand(16, 100, 100).astype(np.float32)
        >>> z = mann_kendall_z(y)
        >>> significant = np.abs(z) > 1.96
        >>> print(f"显著趋势像元比例: {significant.sum() / significant.size * 100}%")
    """
    if y.ndim != 3:
        raise ValueError(f"y must be (T,H,W), got {y.shape}")
    T, h, w = y.shape
    y = y.astype(np.float32, copy=False)
    # 只有这些像素才算：每年都有值
    valid_all = np.all(np.isfinite(y), axis=0)  # (h,w)
    z = np.full((h, w), np.nan, dtype=np.float32)
    s = mann_kendall_S(y)
    base_num = T * (T - 1) * (2 * T + 5)  # 注意：这里先不除 18
    denom = 18.0

    # per-pixel Var(S) with ties
    varS = np.full((h, w), np.nan, dtype=np.float32)
    for r, c in np.argwhere(valid_all):
        vals = y[:, r, c]
        _, counts = np.unique(vals, return_counts=True)
        tie_sum = 0.0
        for t in counts:
            if t > 1:
                tie_sum += t * (t - 1) * (2 * t + 5)
        varS[r, c] = (base_num - tie_sum) / denom

    # 避免 varS<=0 导致 sqrt 问题（例如所有值都相等，varS=0）
    ok = valid_all & (varS > 0)
    z[ok] = s[ok] / np.sqrt(varS[ok])
    return z
def kendall_tau_b_z(y: np.ndarray) -> np.ndarray:
    """计算Kendall tau-b相关系数并转换为z-score。

    计算时间序列与时间的Kendall tau-b相关系数，然后使用Equation 6将其转换为z-score。
    时间被视为严格递增序列（无并列），y值可能包含并列值。

    Args:
        y: 时间序列数据数组，形状为 (T, h, w)，其中T为时间维度，
           h和w为空间维度。数据类型为float32，NaN表示无效值。

    Returns:
        np.ndarray: Z-score数组，形状为 (h, w)，数据类型为float32。
                   使用Equation 6将Kendall tau-b转换为z-score：
                   z = (3 * τ * sqrt(N * (N - 1))) / (2 * (2N + 5))
                   - |z| > 1.96 表示在95%置信水平下显著
                   - 正值表示上升趋势，负值表示下降趋势

    Raises:
        ValueError: 如果输入数组不是3维（T, H, W）形状。

    Note:
        - tau_b考虑了y值中的并列值
        - tau_b = S / sqrt(P * (P - Ty))，其中：
          - S: 一致的符号对减去不一致的符号对
          - P: 总对数 = T(T-1)/2
          - Ty: y值中的并列值对数
        - 使用Equation 6将tau-b转换为z-score：
          z = (3 * τ * sqrt(N * (N - 1))) / (2 * (2N + 5))
        - 只对每年都有有效值的像元进行计算
        - 当时间序列在时间域上单调（无并列值）时，此方法等价于Mann-Kendall Z值

    Example:
        >>> y = np.random.rand(16, 100, 100).astype(np.float32)
        >>> z = kendall_tau_b_z(y)
        >>> significant = np.abs(z) > 1.96
        >>> print(f"显著趋势像元比例: {significant.sum() / significant.size * 100}%")
    """
    if y.ndim != 3:
        raise ValueError(f"y must be (T,H,W), got {y.shape}")
    T, h, w = y.shape
    y = y.astype(np.float32, copy=False)

    valid_all = np.all(np.isfinite(y), axis=0)
    out = np.full((h, w), np.nan, dtype=np.float32)
    if not np.any(valid_all):
        return out

    # S per pixel (same as MK S)
    S = np.zeros((h, w), dtype=np.float32)
    for i in range(T - 1):
        yi = y[i]
        for j in range(i + 1, T):
            S += np.sign(y[j] - yi).astype(np.float32)

    P = T * (T - 1) / 2.0  # total pairs

    # per-pixel Ty
    Ty = np.zeros((h, w), dtype=np.float32)
    for r, c in np.argwhere(valid_all):
        vals = y[:, r, c]
        _, counts = np.unique(vals, return_counts=True)
        Ty_rc = 0.0
        for t in counts:
            if t > 1:
                Ty_rc += t * (t - 1) / 2.0
        Ty[r, c] = Ty_rc

    # 计算Kendall tau-b
    denom = np.sqrt(P * (P - Ty))  # (h,w)
    ok = valid_all & (denom > 0)
    tau_b = np.zeros((h, w), dtype=np.float32)
    tau_b[ok] = S[ok] / denom[ok]
    
    # 使用Equation 6将tau-b转换为z-score
    # z = (3 * τ * sqrt(N * (N - 1))) / (2 * (2N + 5))
    # 其中 N = T (样本数)
    N = float(T)
    numerator = 3.0 * tau_b * np.sqrt(N * (N - 1.0))
    denominator = 2.0 * (2.0 * N + 5.0)
    
    out[ok] = numerator[ok] / denominator
    
    # denom==0 的情况：所有值都一样 => tau 不定义，保持为 NaN
    return out
#-------------------------------------------------------------------------------
#state functions
#-------------------------------------------------------------------------------
def distribution_consistent(y: np.ndarray, min_n0: int = 10, min_n1: int = 3, eps: float = 1e-6):
    """计算两时期分布的统计一致性Z值。

    将时间序列分为两个时期（前13年和后3年），
    检验后一时期的均值是否与前一时期分布一致。

    Args:
        y: 时间序列数据数组，必须是16年的数据，形状为 (16, H, W)，
           其中H和W为空间维度。NaN表示缺失值。
        min_n0: 第一时期（前13年）所需的最小有效值数量，默认为10。
        min_n1: 第二时期（后3年）所需的最小有效值数量，默认为3。
        eps: 标准差的最小阈值，避免除以接近0的值，默认为1e-6。

    Returns:
        tuple: 包含两个数组的元组：
            - z (np.ndarray): Z值数组，形状为 (H, W)，数据类型为float32。
                             z = (mean1 - mean0) / (std0 / sqrt(n1))
                             用于检验后一时期均值是否显著偏离前一时期分布。
            - prod_diff (np.ndarray): 生产力差异数组，形状为 (H, W)，
                                      prod_diff = mean1 - mean0。

    Raises:
        ValueError: 如果输入数组不是3维或时间维度不是16年。

    Note:
        - H0（第一时期）: 前13年 (y[:13])
        - H1（第二时期）: 后3年 (y[13:])
        - Z值基于单样本t检验的标准化形式
        - 只对满足最小样本要求的像元进行计算

    Example:
        >>> y = np.random.rand(16, 100, 100).astype(np.float32)
        >>> z, diff = distribution_consistent(y)
        >>> significant = np.abs(z) > 1.96
        >>> print(f"显著变化像元: {significant.sum()}")
    """
    if y.ndim != 3:
        raise ValueError(f"y must be (T,H,W), got {y.shape}")
    T, H, W = y.shape
    if T != 16:
        raise ValueError(f"y must be 16 years, got {T}")

    y = y.astype(np.float32, copy=False)

    # 创建掩码：只对T维度上所有值都有效的像素进行计算
    # 检查每个像素在时间维度上是否所有值都有效（不是NaN且是有限数）
    valid_time_mask = np.all(np.isfinite(y), axis=0)  # (H, W)

    y0 = y[:13]
    y1 = y[13:]

    # 只对T全满的像素计算统计量，其他像素设为NaN
    mean0 = np.full((H, W), np.nan, dtype=np.float32)
    std0 = np.full((H, W), np.nan, dtype=np.float32)
    mean1 = np.full((H, W), np.nan, dtype=np.float32)
    
    if valid_time_mask.any():
        # 只对T全满的像素计算
        mean0[valid_time_mask] = np.nanmean(y0[:, valid_time_mask], axis=0)
        std0[valid_time_mask] = np.nanstd(y0[:, valid_time_mask], axis=0, ddof=0)
        mean1[valid_time_mask] = np.nanmean(y1[:, valid_time_mask], axis=0)

    prod_diff = mean1 - mean0

    # gate：保证分母和样本数合理，并且只对T全满的像素进行计算
    # T全满意味着n0=13, n1=3（对于16年的数据）
    gate = valid_time_mask & np.isfinite(std0) & (std0 > eps)

    z = np.full((H, W), np.nan, dtype=np.float32)
    # n1 = 3（后3年的数据），因为只对T全满的像素计算
    z[gate] = prod_diff[gate] / (std0[gate] / np.sqrt(3.0))

    return z, prod_diff
#-------------------------------------------------------------------------------
#performance functions
#-------------------------------------------------------------------------------
#对于每个像素点的16个序列，取平均 作为observed npp
#计算每一个分区的90 percentile的值，返回一个维度一样的矩阵
#里面存储每个像素点应该对应的90 percentile的值
def calculate_90_quantile(mean_npp: np.ndarray, lceu_id: np.ndarray, nodata_lceu: int = 0) -> np.ndarray:
    """计算每个LCEU分区的90分位数NPP值。

    为每个LCEU（土地覆盖生态单元）分区计算90分位数NPP值，
    并将该值赋给该分区内的所有像元。

    Args:
        mean_npp: 每像素的评估期多年平均生产力数组，
                 形状为 (H, W)，数据类型为float32。
        lceu_id: 每像素所属的LCEU分区编号数组，
                形状为 (H, W)，数据类型为整数。
        nodata_lceu: LCEU数据中的无数据值，默认为0。

    Returns:
        np.ndarray: 90分位数NPP数组，形状为 (H, W)，数据类型为float32。
                   每个像元的值等于其所属LCEU分区的90分位数。

    Raises:
        ValueError: 如果mean_npp不是2维数组，或lceu_id形状与mean_npp不匹配。

    Note:
        - LCEU分区用于定义生态相似的区域
        - 90分位数用于定义"高生产力"的阈值
        - 同一LCEU内的所有像元共享相同的90分位数值

    Example:
        >>> mean_npp = np.random.rand(100, 100).astype(np.float32)
        >>> lceu_id = np.random.randint(1, 5, (100, 100))
        >>> p90 = calculate_90_quantile(mean_npp, lceu_id)
        >>> print(f"P90值范围: [{p90.min():.3f}, {p90.max():.3f}]")
    """
    if mean_npp.ndim != 2:
        raise ValueError(f"mean_npp must be (H,W), got {mean_npp.shape}")
    if lceu_id.shape != mean_npp.shape:
        raise ValueError("lceu_id shape must match mean_npp")

    H, W = mean_npp.shape
    mean_flat = mean_npp.reshape(-1).astype(np.float32, copy=False)
    lceu_flat = lceu_id.reshape(-1).astype(np.int32, copy=False)

    # 有效像素：有 LCEU + 有 mean 值
    valid = (lceu_flat != nodata_lceu) & np.isfinite(mean_flat)
    mean_v = mean_flat[valid]
    lceu_v = lceu_flat[valid]

    # 逐 LCEU 计算 P90
    p90_by_gid = {}
    for gid in np.unique(lceu_v):
        vals = mean_v[lceu_v == gid]
        # 只对有足够有效值的LCEU计算分位数（至少2个值）
        if len(vals) >= 2:
            # 过滤掉NaN值，只使用有效值
            vals_valid = vals[np.isfinite(vals)]
            if len(vals_valid) >= 2:
                p90_by_gid[int(gid)] = float(np.nanpercentile(vals_valid, 90))

    # 回填到每个像素
    p90_flat = np.full(mean_flat.shape, np.nan, dtype=np.float32)
    for gid, p90 in p90_by_gid.items():
        p90_flat[lceu_flat == gid] = p90

    return p90_flat.reshape(H, W)
def performance_evaluation(y: np.ndarray,lceu_id: np.ndarray,    # (H,W) LCEU 分区编号
    nodata_lceu: int = 0,):
    """评估土地生产力性能指标。

    计算每个像元的平均NPP相对于其所在LCEU分区90分位数的比值，
    用于评估土地的生产力性能水平。

    Args:
        y: 时间序列NPP数据数组，形状为 (T, H, W)，其中T为时间维度，
           H和W为空间维度。数据类型为float32。
        lceu_id: 每像素所属的LCEU分区编号数组，形状为 (H, W)，
                数据类型为整数。
        nodata_lceu: LCEU数据中的无数据值，默认为0。

    Returns:
        np.ndarray: 性能比率数组，形状为 (H, W)，数据类型为float32。
                   Performance = mean_npp / P90_lceu
                   - > 1.0: 高于分区90分位数（高性能）
                   - < 1.0: 低于分区90分位数（低性能）
                   - LCEU为nodata的像元值为NaN

    Raises:
        ValueError: 如果输入数组不是3维（T, H, W）形状。

    Note:
        - mean_npp是时间维度的平均值
        - P90_lceu是每个LCEU分区的90分位数
        - 用于评估相对于同类生态区域的性能水平

    Example:
        >>> y = np.random.rand(16, 100, 100).astype(np.float32)
        >>> lceu_id = np.random.randint(1, 5, (100, 100))
        >>> perf = performance_evaluation(y, lceu_id)
        >>> high_perf = perf > 1.0
        >>> print(f"高性能像元比例: {high_perf.sum() / high_perf.size * 100}%")
    """
    if y.ndim != 3:
        raise ValueError(f"y must be (T,H,W), got {y.shape}")

    T, H, W = y.shape
    y_float = y.astype(np.float32, copy=False)
    
    # 创建掩码：只对T维度上所有值都有效的像素进行计算
    # 检查每个像素在时间维度上是否所有值都有效（不是NaN且是有限数）
    valid_time_mask = np.all(np.isfinite(y_float), axis=0)  # (H, W)
    
    # 只对T为满的像素计算mean_npp，其他像素设为NaN
    mean_npp = np.full((H, W), np.nan, dtype=np.float32)
    if valid_time_mask.any():
        mean_npp[valid_time_mask] = np.nanmean(y_float[:, valid_time_mask], axis=0)
    
    quantile_npp = calculate_90_quantile(mean_npp, lceu_id, nodata_lceu=nodata_lceu)  # (H,W)

    ratio = mean_npp / quantile_npp
    ratio[~np.isfinite(ratio)] = np.nan
    ratio[lceu_id == nodata_lceu] = np.nan
    return ratio.astype(np.float32)
#-------------------------------------------------------------------------------
#final report functions
#-------------------------------------------------------------------------------
def generate_z_or_ratio_results(args: BasicArgs):
    """生成Z值或性能比率结果文件。

    根据评估指标类型（trend/state/performance），生成相应的
    Z值或性能比率栅格文件。支持分块处理以节省内存。

    Args:
        args: BasicArgs对象，包含所有评估参数配置，包括：
            - metrics: 指标类型（"trend"/"state"/"performance"）
            - npp_args: NPP数据参数
            - trend_args/state_args: 趋势/状态参数配置
            - out_dir: 输出目录

    Returns:
        Union[Path, Tuple[Path, Path]]: 
            - 对于performance: 返回performance_ratio.tif路径
            - 对于trend/state: 返回(zscore.tif路径, slope.tif或prod_diff.tif路径)的元组

    Note:
        - trend指标: 生成zscore.tif（Mann-Kendall Z或Kendall tau-b）
                    和slope.tif（Theil-Sen斜率）
        - state指标: 生成zscore.tif（分布一致性Z）和prod_diff.tif（生产力差异）
        - performance指标: 生成performance_ratio.tif（性能比率），
                          需要读取lecu.tif作为输入
        - 使用512x512分块处理，支持大文件
        - 输出格式为float32，nodata为NaN

    Example:
        >>> args = initalize_args("trend", 2023)
        >>> z_path, slope_path = generate_z_or_ratio_results(args)
        >>> print(f"Z值文件: {z_path}, 斜率文件: {slope_path}")
    """
    metrics = args.metrics
    out_dir = args.out_dir

    if metrics == "trend":
        out_aux = out_dir / "slope.tif"
        out_z_path = out_dir / "zscore.tif"          # 输出：要写
        z_is_input = False
    elif metrics == "state":
        out_aux = out_dir / "prod_diff.tif"
        out_z_path = out_dir / "zscore.tif"          # 输出：要写
        z_is_input = False
    elif metrics == "performance":
        out_aux = out_dir / "performance_ratio.tif"  # 输出：ratio
        out_z_path = args.lecu_path                  # 输入：要读
        z_is_input = True
    else:
        raise ValueError(f"Unknown metrics: {metrics}")

    tif_path = args.npp_args.npp_path
    npp_args = args.npp_args
    profile = npp_args.npp_profile
    nodata = npp_args.nodata
    height = npp_args.height
    width = npp_args.width

    profile_float = profile.copy()
    profile_float.update(count=1, dtype="float32", nodata=np.nan)

    block_h, block_w = 512, 512
    n_rows = math.ceil(height / block_h)
    n_cols = math.ceil(width / block_w)

    # --- 1) 先打开 NPP 输入 ---
    with rasterio.open(tif_path) as src:
        # --- 2) 打开 ratio/slope/prod_diff 输出 ---
        with rasterio.open(out_aux, "w", **profile_float) as dst_output:
            # --- 3) 根据 metrics 决定 z 是输入还是输出 ---
            if z_is_input:
                # performance：z 是 lecu 输入，只读
                with rasterio.open(out_z_path, "r") as zsrc:
                    for br in tqdm(range(n_rows), desc="generate performance results"):
                        for bc in range(n_cols):
                            r0 = br * block_h
                            c0 = bc * block_w
                            h = min(block_h, height - r0)
                            w = min(block_w, width - c0)
                            win = Window(c0, r0, w, h)

                            y = read_block_to_float(src, win, nodata)
                            x = zsrc.read(1, window=win).astype(np.float32, copy=False)
                            ratio = performance_evaluation(y, x)
                            dst_output.write(ratio, 1, window=win)

            else:
                # trend/state：z 是输出，要写
                with rasterio.open(out_z_path, "w", **profile_float) as dst_z:
                    for br in tqdm(range(n_rows), desc=f"generate {metrics} results"):
                        for bc in range(n_cols):
                            r0 = br * block_h
                            c0 = bc * block_w
                            h = min(block_h, height - r0)
                            w = min(block_w, width - c0)
                            win = Window(c0, r0, w, h)

                            y = read_block_to_float(src, win, nodata)

                            if metrics == "trend":
                                slope = theil_sen_slope_block(y)
                                if not args.trend_args.mask:
                                    z = mann_kendall_z(y)
                                else:
                                    z = kendall_tau_b_z(y)

                                dst_z.write(z, 1, window=win)
                                dst_output.write(slope, 1, window=win)

                            elif metrics == "state":
                                z, prod_diff = distribution_consistent(y)
                                dst_z.write(z, 1, window=win)
                                dst_output.write(prod_diff, 1, window=win)


    if metrics == "performance":
        return out_aux   # out_aux 是 performance_ratio.tif
    else:
        return out_z_path, out_aux

#-------------------------------------------------------------------------------
# prod_diff_eq
def generate_report_tables(z_eq,npp_eq,mask_eq,args:BasicArgs):
    """生成评估报告表格和分类结果。

    基于重投影后的Z值、NPP和掩膜数据，计算各酋长国的
    退化和改善土地面积，并生成统计表格。

    Args:
        z_eq: 重投影后的Z值栅格文件路径（等面积坐标系）。
        npp_eq: 重投影后的NPP栅格文件路径（等面积坐标系）。
        mask_eq: 掩膜栅格文件路径（如斜率或差异掩膜），
                对于performance指标可能为None。
        args: BasicArgs对象，包含报告参数，包括：
            - metrics: 指标类型
            - reporting_year: 报告年份
            - shp_path: 酋长国边界shapefile路径
            - dst_crs: 目标坐标系
            - name_field: shapefile中的名称字段
            - z_edges: Z值分箱边界
            - mask_low_productivity: 是否应用生产力掩膜
            - trend_args/state_args: 掩膜参数配置

    Returns:
        tuple: 包含三个元素的元组：
            - degrading (np.ndarray): 退化土地掩膜数组，布尔类型，形状为(H,W)。
            - improving (np.ndarray): 改善土地掩膜数组，布尔类型，形状为(H,W)。
            - result (dict): 汇总结果字典，包含：
                - "Target year": 报告年份
                - "Area of Degrading Land (km²)": 总退化面积（平方公里）
                - "Area of Improving Land (km²)": 总改善面积（平方公里）

    Note:
        - 根据z_edges将Z值分箱，bin_id=1表示退化，bin_id=N-1表示改善
        - 应用多种掩膜（生产力掩膜、斜率掩膜、差异掩膜）筛选有效像元
        - 按酋长国统计退化和改善面积，打印详细表格
        - 使用等面积坐标系确保面积计算准确

    Example:
        >>> args = initalize_args("trend", 2023)
        >>> degrading, improving, result = generate_report_tables(
        ...     "z_eq.tif", "npp_eq.tif", "slope_eq.tif", args
        ... )
        >>> print(result)
    """
    metrics = args.metrics
    threshold = None
    if metrics == "trend":
        mask_slope_or_small_diff = args.trend_args.mask
    elif metrics == "state":
        mask_slope_or_small_diff = args.state_args.mask
        threshold = args.state_args.slope_threshold
     # --------- 读取 args（不重新初始化 args）---------
    reporting_year = args.reporting_year
    out_dir = args.out_dir
    shp_path = args.shp_path
    dst_crs = args.dst_crs
    NAME_FIELD = args.name_field
    Z_EDGES = args.z_edges  # 注意你 args 里是 z_edges（小写）
    mask_low_productivity = args.mask_low_productivity

    # --------- 读取矢量边界---------
    emirates = gpd.read_file(shp_path).to_crs(dst_crs)
    emirates["geom_area_km2"] = emirates.geometry.area / 1e6

    # 定义输出文件路径
    out_degrad_path = out_dir / "degrading.tif"
    out_improve_path = out_dir / "improving.tif"
    out_bins_path = out_dir / "bins.tif"
    out_overall_path = out_dir / "overall.tif"
    out_emirates_table_path = out_dir / f"{reporting_year}_emirates_table.csv"
    out_land_cover_table_path = out_dir / f"{reporting_year}_land_cover_table.csv"

    with rasterio.open(z_eq) as zsrc, rasterio.open(npp_eq) as npp_src:
        z = zsrc.read(1).astype(np.float32)
        px_area_km2 = pixel_area_km2(zsrc)
        # 获取profile用于保存输出文件
        profile_output = zsrc.profile.copy()
        profile_output.update(
            dtype="uint8",
            nodata=0,
            compress="lzw",
            tiled=True,
            blockxsize=512,
            blockysize=512,
        )
        profile_bins = zsrc.profile.copy()
        profile_bins.update(
            dtype="int16",
            nodata=99,
            compress="lzw",
            tiled=True,
            blockxsize=512,
            blockysize=512,
        )
        # nodata->nan（我们写的时候 dst_nodata=nan）
        z = np.where(np.isfinite(z), z, np.nan)
        valid = np.isfinite(z)
        # 加入 productive mask（关键）
        if mask_low_productivity:
            #返回0/1
            prod_mask = mask_low_productivity(npp_src.read())
            valid = valid & prod_mask
        if metrics == "trend":
            mask_slope = args.trend_args.mask
            slope_eq = mask_eq
            if mask_slope:
         
                slope_threshold = args.trend_args.slope_threshold
                if slope_eq is None:
                    raise ValueError(
                        "trend_args.mask=True 需要提供 slope_eq（用于 mask_small_slope）。"
                    )
                with rasterio.open(slope_eq) as ssrc:
                    s = ssrc.read(1).astype(np.float32)
                    prod_slop_mask = mask_small_slope(s,slope_threshold) 
                    valid = valid & prod_slop_mask.astype(bool)
        elif metrics == "state":
            mask_small_diff = args.state_args.mask
            threshold = args.state_args.slope_threshold
            if mask_small_diff:
                if mask_eq is None:
                    raise ValueError("state_args.mask=True 需要提供 mask_eq。")
                with rasterio.open(mask_eq) as mask_src:
                    small_diff_mask = mask_src.read(1) > threshold
                valid = valid & small_diff_mask.astype(bool)
        
        gate = valid
        #boolen 的 ndarray，也是在等面积坐标系底下生成的
        # gate 是 0/1 矩阵（或 bool），z 是同形状浮点数组
        new_z = np.where(gate == 1, z, np.nan).astype(np.float32, copy=False)  
        bins = classify_bins(new_z,Z_EDGES)
        degrading = (bins == 1)
        improving = (bins == len(Z_EDGES)-1)
        # 计算稳定区域（既不是退化也不是改善的有效像素）
        stable = (bins >= 2) & (bins < len(Z_EDGES)-1) & (bins != 99)
        # 计算无数据区域（无效像素）
        no_data = (bins == 99) | ~np.isfinite(z)
        nodata_code = np.uint8(99)
        overall = np.full(bins.shape, nodata_code, dtype=np.uint8)

        # 只在有效像素写入分类
        overall[degrading & ~no_data] = np.uint8(2)
        overall[stable   & ~no_data] = np.uint8(0)
        overall[improving & ~no_data] = np.uint8(1)
        
        rows = []
        total_degradation = 0
        total_improvement = 0
       
        for _, row in tqdm(emirates.iterrows(), total=len(emirates), desc="Emirates"):
            reg_mask = geometry_mask(
                [row.geometry],
                out_shape=(zsrc.height, zsrc.width),
                transform=zsrc.transform,
                invert=True
            )
            deg_km2 = float(np.count_nonzero(degrading & reg_mask) * px_area_km2)
            imp_km2 = float(np.count_nonzero(improving & reg_mask) * px_area_km2)
            total_degradation += deg_km2
            total_improvement += imp_km2
            
            rows.append({
            "Emirate": row[NAME_FIELD],
            "Degrading land (km²)": deg_km2,
            "Improving land (km²)": imp_km2,
        })
        
        # 计算全区域的土地覆盖变化统计
        # 创建全区域掩膜（所有酋长国的并集）
        total_mask = geometry_mask(
            emirates.geometry.tolist(),
            out_shape=(zsrc.height, zsrc.width),
            transform=zsrc.transform,
            invert=True
        )
        
        # 计算各类别的面积
        improved_area_km2 = float(np.count_nonzero(improving & total_mask) * px_area_km2)
        stable_area_km2 = float(np.count_nonzero(stable & total_mask) * px_area_km2)
        degraded_area_km2 = float(np.count_nonzero(degrading & total_mask) * px_area_km2)
        no_data_area_km2 = float(np.count_nonzero(no_data & total_mask) * px_area_km2)
        
        # 计算总土地面积（包括有效像素和无数据）
        total_area_km2 = improved_area_km2 + stable_area_km2 + degraded_area_km2 + no_data_area_km2
        
        # 计算百分比（基于总土地面积，包括无数据）
        if total_area_km2 > 0:
            improved_percent = (improved_area_km2 / total_area_km2) * 100
            stable_percent = (stable_area_km2 / total_area_km2) * 100
            degraded_percent = (degraded_area_km2 / total_area_km2) * 100
            no_data_percent = (no_data_area_km2 / total_area_km2) * 100
        else:
            improved_percent = 0.0
            stable_percent = 0.0
            degraded_percent = 0.0
            no_data_percent = 0.0
        
        # 生成土地覆盖变化统计表格
        land_cover_table = pd.DataFrame({
            "Land cover change category": [
                "Land area with improved land cover",
                "Land area with stable land cover",
                "Land area with degraded land cover",
                "Land area with no land cover data"
            ],
            "Area (km²)": [
                round(improved_area_km2, 1),
                round(stable_area_km2, 1),
                round(degraded_area_km2, 1),
                round(no_data_area_km2, 1)
            ],
            "Percent of total land area (%)": [
                round(improved_percent, 1),
                round(stable_percent, 1),
                round(degraded_percent, 1),
                round(no_data_percent, 1)
            ]
        })
        
        result = {
                "Target year": reporting_year,
                "Area of Degrading Land (km²)": total_degradation,
                "Area of Improving Land (km²)": total_improvement,
            }
        table36 = pd.DataFrame(rows).sort_values("Emirate").reset_index(drop=True)
        print(table36)
        table36.to_csv(out_emirates_table_path, index=False)
        land_cover_table.to_csv(out_land_cover_table_path, index=False)
        
        # 打印土地覆盖变化统计表格
        print("\n" + "=" * 80)
        print("Land Cover Change Statistics")
        print("=" * 80)
        print(land_cover_table.to_string(index=False))
        print("=" * 80)
        
        # 保存 degrading 和 improving 文件
        with rasterio.open(out_degrad_path, "w", **profile_output) as dst_output:
            dst_output.write(degrading.astype(np.uint8), 1)
        with rasterio.open(out_improve_path, "w", **profile_output) as dst_output:
            dst_output.write(improving.astype(np.uint8), 1)
        
        # 保存 bins 文件（用于后续组合分析）
        with rasterio.open(out_bins_path, "w", **profile_bins) as dst_output:
            dst_output.write(bins.astype(np.int16), 1)
        with rasterio.open(out_overall_path, "w", **profile_bins) as dst_output:
            dst_output.write(overall.astype(np.uint8), 1)
        return degrading, improving, result
#-------------------------------------------------------------------------------

def generate_report_period(metrics:str,args:BasicArgs):
    """生成完整的评估报告周期数据。

    执行完整的评估流程：重投影NPP数据、生成Z值/比率结果、
    重投影结果数据、计算退化和改善面积统计。

    Args:
        metrics: 评估指标类型，必须是 "trend"、"state" 或 "performance"。
        args: BasicArgs对象，包含所有评估参数配置。

    Returns:
        dict: 汇总结果字典，包含：
            - "Target year": 报告年份
            - "Area of Degrading Land (km²)": 总退化面积（平方公里）
            - "Area of Improving Land (km²)": 总改善面积（平方公里）

    Note:
        执行流程：
        1. 将NPP数据重投影到等面积坐标系
        2. 生成Z值/斜率/比率结果（基于原始坐标系）
        3. 将结果重投影到等面积坐标系
        4. 生成报告表格和统计结果

        输出文件：
        - NPP重投影文件：{npp_path.stem}NPP_PROXY_eqarea6933.tif
        - Z值重投影文件：{zscore.stem}_eqarea6933.tif
        - 斜率/差异重投影文件：{slope.stem}_slope_eqarea6933.tif

    Raises:
        AssertionError: 如果metrics不是 "trend"、"state" 或 "performance"。

    Example:
        >>> args = initalize_args("trend", 2023)
        >>> result = generate_report_period("trend", args)
        >>> print(f"退化面积: {result['Area of Degrading Land (km²)']} km²")
    """
    # basic arguments
    dst_crs = args.dst_crs
    dst_res = args.dst_res

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    
    npp_path = args.npp_args.npp_path
    # step one: reproject npp to equal area and generate productivity mask
    npp_eq = out_dir / (npp_path.stem + "NPP_PROXY_eqarea6933.tif")
    reproject_to_equal_area(npp_path, npp_eq, dst_crs, dst_res, Resampling.bilinear)
    #step two: generate slop,z-score,tau-a,N
    out_z:Path
    z_eq:Path
    assert metrics in ("trend","state","performance")
    output  = generate_z_or_ratio_results(args)
    s_or_prod_diff_eq = None
    if isinstance(output, tuple):
        out_z,out_slope_or_prod_diff_or_ratio = output
        z_eq = out_dir /(out_z.stem + "_eqarea6933.tif")
        reproject_to_equal_area(out_z, z_eq, dst_crs, dst_res, Resampling.bilinear)
        s_or_prod_diff_eq = out_dir /(out_slope_or_prod_diff_or_ratio.stem + "_slope_eqarea6933.tif")
        reproject_to_equal_area(out_slope_or_prod_diff_or_ratio, s_or_prod_diff_eq, dst_crs, dst_res, Resampling.bilinear)
    elif isinstance(output, Path):
        out_ratio = output
        z_eq = out_dir /(out_ratio.stem + "_eqarea6933.tif")
        reproject_to_equal_area(out_ratio, z_eq, dst_crs, dst_res, Resampling.bilinear)

    #step four: calculate the area of the degrading pixels
    degrading, improving,result = generate_report_tables(z_eq,npp_eq,s_or_prod_diff_eq,args)
   
    return result

#-------------------------------------------------------------------------------
def combine_metrics_degradation(
    trend_degrading_path: Path,
    state_degrading_path: Path,
    performance_degrading_path: Path,
    output_path: Path,
    table_version: str = "4-4",
) -> Path:
    """组合三个指标的退化判断结果，根据查找表（Table 4-4或4-5）判断每个像素是否最终退化。

    根据Table 4-4或Table 4-5的查找表规则，结合trend、state、performance三个指标
    的退化判断结果（degrading.tif），判断每个像素点是否最终退化。只对三个指标都有值的像素进行比较。

    Args:
        trend_degrading_path: Trend指标的退化判断结果文件路径（uint8，1=退化，0=非退化）。
        state_degrading_path: State指标的退化判断结果文件路径（uint8，1=退化，0=非退化）。
        performance_degrading_path: Performance指标的退化判断结果文件路径（uint8，1=退化，0=非退化）。
        output_path: 输出最终退化判断结果文件路径（uint8，1=退化，0=非退化）。
        table_version: 查找表版本，可选"4-4"或"4-5"，默认为"4-4"。

    Returns:
        Path: 输出文件路径。

    Note:
        - 查找表规则（Table 4-4）：
          - Class 1: Trend(Y), State(Y), Performance(Y) -> Degraded(Y)
          - Class 2: Trend(Y), State(Y), Performance(N) -> Degraded(Y)
          - Class 3: Trend(Y), State(N), Performance(Y) -> Degraded(Y)
          - Class 4: Trend(Y), State(N), Performance(N) -> Degraded(Y)
          - Class 5: Trend(N), State(Y), Performance(Y) -> Degraded(Y)
          - Class 6: Trend(N), State(Y), Performance(N) -> Degraded(N)
          - Class 7: Trend(N), State(N), Performance(Y) -> Degraded(N)
          - Class 8: Trend(N), State(N), Performance(N) -> Degraded(N)

        - 查找表规则（Table 4-5，相对于4-4的变化）：
          - Class 4: Trend(Y), State(N), Performance(N) -> Degraded(N) [变化]
          - Class 7: Trend(N), State(N), Performance(Y) -> Degraded(Y) [变化]

        - 只对三个指标都有有效值的像素进行比较（排除nodata值）。
        - 输入文件应为uint8类型，1表示退化，0表示非退化或无效值（nodata=0）。
        - 输出为uint8类型，1表示最终判定为退化，0表示非退化或无效值。

    Example:
        >>> trend_deg = Path("../output/trend/2023/degrading.tif")
        >>> state_deg = Path("../output/state/2023/degrading.tif")
        >>> perf_deg = Path("../output/performance/2023/degrading.tif")
        >>> output = Path("../output/combined_degradation_2023_table44.tif")
        >>> result_path = combine_metrics_degradation(
        ...     trend_deg, state_deg, perf_deg, output, table_version="4-5"
        ... )
    """
    # 读取三个degrading文件
    with rasterio.open(trend_degrading_path) as src_trend:
        trend_degrading = src_trend.read(1).astype(np.uint8)
        profile = src_trend.profile.copy()
        trend_nodata = src_trend.nodata if src_trend.nodata is not None else 0
        
        with rasterio.open(state_degrading_path) as src_state:
            state_degrading = src_state.read(1).astype(np.uint8)
            state_nodata = src_state.nodata if src_state.nodata is not None else 0
            
            with rasterio.open(performance_degrading_path) as src_perf:
                performance_degrading = src_perf.read(1).astype(np.uint8)
                perf_nodata = src_perf.nodata if src_perf.nodata is not None else 0
                
                # 检查尺寸是否一致
                if (trend_degrading.shape != state_degrading.shape or 
                    trend_degrading.shape != performance_degrading.shape):
                    raise ValueError(
                        f"Degrading arrays must have the same shape. "
                        f"Got: trend {trend_degrading.shape}, "
                        f"state {state_degrading.shape}, "
                        f"performance {performance_degrading.shape}"
                    )
                
                # 标记有效像素
                # 由于 degrading.tif 文件：1=退化，0=非退化或无效值，nodata=0
                # 我们无法直接区分 0 是表示非退化还是无效值
                # 但是，由于 degrading 和 improving 文件是在 generate_report_tables 中生成的，
                # 它们只对有效像素（gate=True）计算了 bins，无效像素的 bins=99
                # 对于 bins=99 的像素，degrading 和 improving 都是 False (0)
                # 
                # 为了判断有效像素，我们可以：
                # 1. 如果三个文件的 nodata 都是 0，则很难区分
                # 2. 更实际的方法：假设所有像素都有效（除非有明确的无效标记）
                # 3. 或者：检查是否有 improving 文件辅助判断（但用户只要求 degrading）
                #
                # 采用策略：假设所有像素都有效，因为 degrading.tif 本身就是从有效像素计算的
                # 如果某个方法对某个像素没有值，那么在生成 degrading.tif 时该像素已经是 0
                # 我们只对三个文件都有值的像素进行比较
                #
                # 实际上，我们可以更保守：只有当值在有效范围内（0 或 1）时认为有效
                # 但由于 uint8 的值只能是 0-255，我们需要检查是否有其他标记
                valid_pixels = np.ones(trend_degrading.shape, dtype=bool)
                
                # 如果 nodata 不是 None，则排除 nodata 值
                if trend_nodata is not None:
                    trend_valid = (trend_degrading != trend_nodata)
                else:
                    trend_valid = np.ones_like(trend_degrading, dtype=bool)
                    
                if state_nodata is not None:
                    state_valid = (state_degrading != state_nodata)
                else:
                    state_valid = np.ones_like(state_degrading, dtype=bool)
                    
                if perf_nodata is not None:
                    perf_valid = (performance_degrading != perf_nodata)
                else:
                    perf_valid = np.ones_like(performance_degrading, dtype=bool)
                
                # 只对三个指标都有有效值的像素进行比较（都不为nodata）
                valid_pixels = trend_valid & state_valid & perf_valid
                
                # 转换为退化状态（Y=True表示退化，N=False表示非退化）
                #  degrading == 0 表示退化(Y)，degrading == 1 表示stable(N)，degrading == 2 表示improving(Y)
                trend_y = (trend_degrading == 0) & valid_pixels
                state_y = (state_degrading == 0) & valid_pixels
                perf_y = (performance_degrading == 0) & valid_pixels

                trend_stable = (trend_degrading == 1) & valid_pixels
                state_stable = (state_degrading == 1) & valid_pixels
                perf_stable = (performance_degrading == 1) & valid_pixels

                trend_improving = (trend_degrading == 2) & valid_pixels
                state_improving = (state_degrading == 2) & valid_pixels
                perf_improving = (performance_degrading == 2) & valid_pixels
                
                # 初始化输出数组（默认0表示非退化或无效）
                final_degraded = np.zeros(trend_degrading.shape, dtype=np.uint8)
                final_class = np.full(trend_degrading.shape, 99, dtype=np.uint8)
                # 只对有效像素进行查找表判断
                if valid_pixels.any():
                    # 获取有效像素的退化状态
                    t_y = trend_y[valid_pixels]
                    s_y = state_y[valid_pixels]
                    p_y = perf_y[valid_pixels]

                    t_i = trend_improving[valid_pixels]
                    s_i = state_improving[valid_pixels]
                    p_i = perf_improving[valid_pixels]

                    # 根据查找表判断是否退化
                    if table_version == "4-4":
                        # Table 4-4规则
                        # Class 1-5: 退化
                        # Class 6-8: 非退化
                        degraded_mask = (
                            (t_y & s_y & p_y) |      # Class 1
                            (t_y & s_y & ~p_y) |     # Class 2
                            (t_y & ~s_y & p_y) |     # Class 3
                            (t_y & ~s_y & ~p_y) |    # Class 4
                            (~t_y & s_y & p_y)       # Class 5
                            # Class 6, 7, 8 都是非退化，不需要额外判断
                        )
                        
                    elif table_version == "4-5":
                        # Table 4-5规则（相对于4-4的变化）
                        # Class 1-3, 5, 7: 退化
                        # Class 4, 6, 8: 非退化
                        degraded_mask = (
                            (t_y & s_y & p_y) |      # Class 1: 退化
                            (t_y & s_y & ~p_y) |     # Class 2: 退化
                            (t_y & ~s_y & p_y) |     # Class 3: 退化
                            # Class 4: Trend(Y), State(N), Performance(N) -> 非退化（变化）
                            (~t_y & s_y & p_y) |     # Class 5: 退化
                            # Class 6: 非退化
                            (~t_y & ~s_y & p_y)      # Class 7: 退化（变化）
                            # Class 8: 非退化
                        )
                    else:
                        raise ValueError(f"Unknown table_version: {table_version}. Must be '4-4' or '4-5'")
                    # === 2) 对“非退化”再分：improving vs stable ===
                    non_degraded = ~degraded_mask
                    improving_mask = non_degraded & (t_i | s_i | p_i)   # 只要任一指标 improving 就算 improving
                    stable_mask    = non_degraded & (~improving_mask)   # 剩下就是 stable

                    # 写入最终分类：0/1/2
                    out = np.empty_like(degraded_mask, dtype=np.uint8)
                    out[degraded_mask]  = 0
                    out[stable_mask]    = 1
                    out[improving_mask] = 2

                    final_class[valid_pixels] = out
                    
                    # 将结果写入输出数组（1=退化，0=非退化）
                    
                
                # 准备输出profile
                output_profile = profile.copy()
                output_profile.update(
                    dtype="uint8",
                    nodata=99,  # 99表示非退化或无效值
                    compress="lzw",
                    tiled=True,
                    blockxsize=512,
                    blockysize=512,
                )
                
                # 写入输出文件
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with rasterio.open(output_path, "w", **output_profile) as dst:
                    dst.write(final_class, 1)
    
    return output_path

#-------------------------------------------------------------------------------
def productivity_land_cover_classification(
    land_cover_path: Path,
    productivity_bins: Path,
    land_cover_class_names: Optional[dict[int, str]] = None,
    dst_crs: str = "EPSG:6933",
    dst_res: float = 250.0,
) -> pd.DataFrame:
    """按土地覆盖类别统计生产力退化/非退化/无数据面积。

    读取土地覆盖分类文件和生产力bins文件，对每个土地覆盖类别统计
    退化、非退化和无数据的面积（km²）。

    Args:
        land_cover_path: 土地覆盖分类栅格文件路径（如 C3S 分类文件）。
        productivity_bins: 生产力bins分类栅格文件路径（如 trend/state/performance的bins.tif）。
        land_cover_class_names: 可选的字典，将土地覆盖类别数值映射到名称。
                                如果为None，则使用数值作为类别名称。
                                例如：{1: "Tree-covered areas", 2: "Grasslands", ...}
        dst_crs: 目标等面积坐标系，默认为"EPSG:6933"。
        dst_res: 目标分辨率（米），默认为250.0。

    Returns:
        pd.DataFrame: 统计表格，包含以下列：
            - Land Cover Class: 土地覆盖类别名称或数值
            - Degraded (km²): 退化土地面积
            - Non-Degraded (km²): 非退化土地面积（包括稳定和改善）
            - No Data (km²): 无数据面积

    Note:
        - bins分类规则：
          - bin_id == 1: 退化（Degraded）
          - bin_id 在 2 到 N-1 之间: 非退化（Non-Degraded，包括稳定和改善）
          - bin_id == 99: 无效值（No Data）
        - land_cover的nodata值也被视为No Data
        - 自动将两个文件都重投影到指定的等面积坐标系和分辨率
        - 使用两个文件的并集边界来确定输出范围，确保覆盖所有数据

    Example:
        >>> land_cover = Path("../data/Land_Cover_Each_Year/C3S_2015_reclass_300m.tif")
        >>> bins = Path("output/trend/2015/bins.tif")
        >>> class_names = {10: "Tree-covered areas", 20: "Grasslands"}
        >>> df = productivity_land_cover_classification(land_cover, bins, class_names)
        >>> print(df)
    """
    from rasterio.warp import reproject, Resampling, calculate_default_transform, transform_bounds
    
    # 步骤1: 读取两个文件的元数据，确定并集边界
    with rasterio.open(land_cover_path) as lc_src, rasterio.open(productivity_bins) as bins_src:
        # 获取源文件边界并转换为目标坐标系
        lc_bounds = transform_bounds(lc_src.crs, dst_crs, *lc_src.bounds)
        bins_bounds = transform_bounds(bins_src.crs, dst_crs, *bins_src.bounds)
        
        # 计算并集边界
        union_bounds = (
            min(lc_bounds[0], bins_bounds[0]),  # min left
            min(lc_bounds[1], bins_bounds[1]),  # min bottom
            max(lc_bounds[2], bins_bounds[2]),  # max right
            max(lc_bounds[3], bins_bounds[3]),  # max top
        )
        
        # 计算目标transform和尺寸（基于并集边界和指定分辨率）
        # 注意：如果两个文件已经在目标CRS中，我们需要手动创建transform
        try:
            dst_transform, dst_width, dst_height = calculate_default_transform(
                dst_crs, dst_crs, 
                width=None, height=None,
                left=union_bounds[0],
                bottom=union_bounds[1],
                right=union_bounds[2],
                top=union_bounds[3],
                resolution=dst_res
            )
        except Exception as e:
            # 如果calculate_default_transform失败，手动创建transform
            from rasterio.transform import from_bounds
            dst_transform = from_bounds(
                union_bounds[0], union_bounds[1], union_bounds[2], union_bounds[3],
                int((union_bounds[2] - union_bounds[0]) / dst_res),
                int((union_bounds[3] - union_bounds[1]) / dst_res)
            )
            dst_width = int((union_bounds[2] - union_bounds[0]) / dst_res)
            dst_height = int((union_bounds[3] - union_bounds[1]) / dst_res)
        
        # 验证transform是否正确创建
        if dst_transform is None:
            raise ValueError(f"Failed to create transform for CRS {dst_crs} with resolution {dst_res}")
        
        # 获取land_cover的nodata值
        lc_nodata = lc_src.nodata if lc_src.nodata is not None else -32768
        bins_nodata = bins_src.nodata if bins_src.nodata is not None else 99
    
    # 步骤2: 重投影land_cover到目标坐标系和分辨率
    with rasterio.open(land_cover_path) as lc_src:
        lc_data = lc_src.read(1)
        lc_data_reprojected = np.full((dst_height, dst_width), lc_nodata, dtype=lc_src.dtypes[0])
        
        reproject(
            source=lc_data,
            destination=lc_data_reprojected,
            src_transform=lc_src.transform,
            src_crs=lc_src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,  # 分类数据使用nearest
            src_nodata=lc_src.nodata,
            dst_nodata=lc_nodata
        )
    
    # 步骤3: 重投影bins到目标坐标系和分辨率
    with rasterio.open(productivity_bins) as bins_src:
        bins_data_src = bins_src.read(1).astype(np.int16)
        bins_data_reprojected = np.full((dst_height, dst_width), bins_nodata, dtype=np.int16)
        
        reproject(
            source=bins_data_src,
            destination=bins_data_reprojected,
            src_transform=bins_src.transform,
            src_crs=bins_src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,  # 分类数据使用nearest
            src_nodata=bins_src.nodata,
            dst_nodata=bins_nodata
        )
        
        # 确定bins的最大值（排除nodata=99）
        valid_bins = bins_data_reprojected[bins_data_reprojected != bins_nodata]
        if len(valid_bins) == 0:
            max_bin_id = 1
        else:
            max_bin_id = int(np.max(valid_bins))
    
    # 步骤4: 计算像元面积（基于等面积坐标系）
    # 直接从transform计算像元面积，避免创建临时文件
    if dst_transform is None:
        raise ValueError(f"Transform is None. Cannot calculate pixel area.")
    
    px_w = dst_transform.a
    px_h = dst_transform.e  # e通常是负数
    
    # 检查 transform 属性是否为 None
    if px_w is None or px_h is None:
        raise ValueError(
            f"Invalid transform: {dst_transform}. "
            f"px_w={px_w}, px_h={px_h}. "
            f"CRS: {dst_crs}, Resolution: {dst_res}"
        )
    
    # 取绝对值（e通常是负数）
    px_w_abs = abs(float(px_w))
    px_h_abs = abs(float(px_h))
    
    # 检查像元尺寸是否合理（避免度为单位的地理坐标）
    if px_w_abs < 0.1 or px_h_abs < 0.1:
        raise ValueError(
            f"Raster appears to be in degrees (geographic CRS). "
            f"Transform pixel size: {px_w_abs} x {px_h_abs}. "
            f"For accurate area, reproject to an equal-area projection first."
        )
    
    # 计算像元面积（km²）
    px_area_km2 = (px_w_abs * px_h_abs) / 1e6
    
    # 验证计算结果
    if px_area_km2 is None or not np.isfinite(px_area_km2):
        raise ValueError(
            f"Failed to calculate pixel area. "
            f"px_w_abs={px_w_abs}, px_h_abs={px_h_abs}, "
            f"px_area_km2={px_area_km2}"
        )
    
    # 获取所有土地覆盖类别（排除nodata）
    unique_lc_classes = np.unique(lc_data_reprojected)
    unique_lc_classes = unique_lc_classes[unique_lc_classes != lc_nodata]
    print(unique_lc_classes)
    
    # 如果没有提供类别名称映射，使用数值
    if land_cover_class_names is None:
        land_cover_class_names = {cls: f"Class {cls}" for cls in unique_lc_classes}
    
    # 统计每个土地覆盖类别的退化/非退化/无数据面积
    results = []
    
    for lc_class in unique_lc_classes:
        # 创建该土地覆盖类别的掩膜
        lc_mask = (lc_data_reprojected == lc_class)
        
        # 在该类别内统计bins（只统计land_cover有效区域内的bins）
        lc_bins = bins_data_reprojected[lc_mask]
        
        # 统计退化（bin_id == 1）
        degraded_mask = (lc_bins == 1)
        degraded_count = np.count_nonzero(degraded_mask)
        degraded_km2 = float(degraded_count) * float(px_area_km2)
        
        # 统计非退化（bin_id 在 2 到 max_bin_id 之间，包括改善）
        non_degraded_mask = (lc_bins >= 2) & (lc_bins < 99)
        non_degraded_count = np.count_nonzero(non_degraded_mask)
        non_degraded_km2 = float(non_degraded_count) * float(px_area_km2)
        
        # 统计无数据（bin_id == 99）
        no_data_mask = (lc_bins == 99)
        no_data_count = np.count_nonzero(no_data_mask)
        no_data_km2 = float(no_data_count) * float(px_area_km2)
        
        # 获取类别名称
        class_name = land_cover_class_names.get(lc_class, f"Class {lc_class}")
        
        results.append({
            "Land Cover Class": class_name,
            "Degraded (km²)": round(degraded_km2, 2),
            "Non-Degraded (km²)": round(non_degraded_km2, 2),
            "No Data (km²)": round(no_data_km2, 2),
        })
    
    # 创建DataFrame并返回
    df = pd.DataFrame(results)
    return df
#-------------------------------------------------------------------------------
def land_cover_conversion_productivity(
    initial_land_cover_path: Path,
    final_land_cover_path: Path,
    productivity_bins: Path,
    land_cover_class_names: Optional[dict[int, str]] = None,
    output_dir: Optional[Path] = None,
    dst_crs: str = "EPSG:6933",
    dst_res: float = 250.0,
    bin_labels: Optional[dict[int, str]] = None,
) -> pd.DataFrame:
    """分析土地覆盖转换及其对应的生产力动态。

    读取初始和最终土地覆盖分类文件以及生产力bins文件，统计每个土地转换组合
    (From -> To) 在不同生产力动态类别下的面积（km²）。

    Args:
        initial_land_cover_path: 初始土地覆盖分类栅格文件路径。
        final_land_cover_path: 最终土地覆盖分类栅格文件路径。
        productivity_bins: 生产力bins分类栅格文件路径（如 trend/state/performance的bins.tif）。
        land_cover_class_names: 可选的字典，将土地覆盖类别数值映射到名称。
                                如果为None，则使用数值作为类别名称。
                                例如：{1: "Tree-covered areas", 2: "Grasslands", ...}
        output_dir: 输出目录路径，如果提供则保存表格到该目录。
                    默认为 None，不保存文件。
        dst_crs: 目标等面积坐标系，默认为"EPSG:6933"。
        dst_res: 目标分辨率（米），默认为250.0。
        bin_labels: 可选的字典，将bin_id映射到标签名称。
                    如果为None，使用默认标签：
                    {1: "Declining", 2: "Moderate Decline", 3: "Stressed", 
                     4: "Stable", 5: "Increasing"}

    Returns:
        pd.DataFrame: 统计表格，包含以下列：
            - From: 初始土地覆盖类别名称
            - To: 最终土地覆盖类别名称
            - Declining (km²): 退化面积
            - Moderate Decline (km²): 中度退化面积
            - Stressed (km²): 压力面积
            - Stable (km²): 稳定面积
            - Increasing (km²): 改善面积

    Note:
        - bins分类规则（基于5个bins的情况）：
          - bin_id == 1: Declining（退化）
          - bin_id == 2: Moderate Decline（中度退化）
          - bin_id == 3: Stressed（压力）
          - bin_id == 4: Stable（稳定）
          - bin_id == 5: Increasing（改善）
          - bin_id == 99: 无效值（No Data，不统计）
        - 只统计有效的土地转换（From 和 To 都有效）
        - 自动将三个文件都重投影到指定的等面积坐标系和分辨率
        - 使用三个文件的并集边界来确定输出范围

    Example:
        >>> initial_lc = Path("../data/Land_Cover_Each_Year/C3S_2000_reclass_300m.tif")
        >>> final_lc = Path("../data/Land_Cover_Each_Year/C3S_2015_reclass_300m.tif")
        >>> bins = Path("output/trend/2015/bins.tif")
        >>> class_names = {1: "Forest", 2: "Grassland", 3: "Cropland"}
        >>> df = land_cover_conversion_productivity(initial_lc, final_lc, bins, class_names)
        >>> print(df)
    """
    from rasterio.warp import reproject, Resampling, calculate_default_transform, transform_bounds
    
    # 默认 bin 标签
    if bin_labels is None:
        bin_labels = {
            1: "Declining",
            2: "Moderate Decline",
            3: "Stressed",
            4: "Stable",
            5: "Increasing"
        }
    
    # 步骤1: 读取三个文件的元数据，确定并集边界
    with rasterio.open(initial_land_cover_path) as lc1_src, \
         rasterio.open(final_land_cover_path) as lc2_src, \
         rasterio.open(productivity_bins) as bins_src:
        
        # 获取源文件边界并转换为目标坐标系
        lc1_bounds = transform_bounds(lc1_src.crs, dst_crs, *lc1_src.bounds)
        lc2_bounds = transform_bounds(lc2_src.crs, dst_crs, *lc2_src.bounds)
        bins_bounds = transform_bounds(bins_src.crs, dst_crs, *bins_src.bounds)
        
        # 计算并集边界
        union_bounds = (
            min(lc1_bounds[0], lc2_bounds[0], bins_bounds[0]),  # min left
            min(lc1_bounds[1], lc2_bounds[1], bins_bounds[1]),  # min bottom
            max(lc1_bounds[2], lc2_bounds[2], bins_bounds[2]),  # max right
            max(lc1_bounds[3], lc2_bounds[3], bins_bounds[3]),  # max top
        )
        
        # 计算目标transform和尺寸
        try:
            dst_transform, dst_width, dst_height = calculate_default_transform(
                dst_crs, dst_crs, 
                width=None, height=None,
                left=union_bounds[0],
                bottom=union_bounds[1],
                right=union_bounds[2],
                top=union_bounds[3],
                resolution=dst_res
            )
        except Exception as e:
            # 如果calculate_default_transform失败，手动创建transform
            from rasterio.transform import from_bounds
            dst_transform = from_bounds(
                union_bounds[0], union_bounds[1], union_bounds[2], union_bounds[3],
                int((union_bounds[2] - union_bounds[0]) / dst_res),
                int((union_bounds[3] - union_bounds[1]) / dst_res)
            )
            dst_width = int((union_bounds[2] - union_bounds[0]) / dst_res)
            dst_height = int((union_bounds[3] - union_bounds[1]) / dst_res)
        
        if dst_transform is None:
            raise ValueError(f"Failed to create transform for CRS {dst_crs} with resolution {dst_res}")
        
        # 获取nodata值
        lc1_nodata = lc1_src.nodata if lc1_src.nodata is not None else -32768
        lc2_nodata = lc2_src.nodata if lc2_src.nodata is not None else -32768
        bins_nodata = bins_src.nodata if bins_src.nodata is not None else 99
    
    # 步骤2: 重投影初始土地覆盖
    with rasterio.open(initial_land_cover_path) as lc1_src:
        lc1_data = lc1_src.read(1)
        lc1_data_reprojected = np.full((dst_height, dst_width), lc1_nodata, dtype=lc1_src.dtypes[0])
        
        reproject(
            source=lc1_data,
            destination=lc1_data_reprojected,
            src_transform=lc1_src.transform,
            src_crs=lc1_src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,
            src_nodata=lc1_src.nodata,
            dst_nodata=lc1_nodata
        )
    
    # 步骤3: 重投影最终土地覆盖
    with rasterio.open(final_land_cover_path) as lc2_src:
        lc2_data = lc2_src.read(1)
        lc2_data_reprojected = np.full((dst_height, dst_width), lc2_nodata, dtype=lc2_src.dtypes[0])
        
        reproject(
            source=lc2_data,
            destination=lc2_data_reprojected,
            src_transform=lc2_src.transform,
            src_crs=lc2_src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,
            src_nodata=lc2_src.nodata,
            dst_nodata=lc2_nodata
        )
    
    # 步骤4: 重投影生产力bins
    with rasterio.open(productivity_bins) as bins_src:
        bins_data_src = bins_src.read(1).astype(np.int16)
        bins_data_reprojected = np.full((dst_height, dst_width), bins_nodata, dtype=np.int16)
        
        reproject(
            source=bins_data_src,
            destination=bins_data_reprojected,
            src_transform=bins_src.transform,
            src_crs=bins_src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,
            src_nodata=bins_src.nodata,
            dst_nodata=bins_nodata
        )
    
    # 步骤5: 计算像元面积
    px_w = abs(float(dst_transform.a))
    px_h = abs(float(dst_transform.e))
    
    if px_w < 0.1 or px_h < 0.1:
        raise ValueError(
            f"Raster appears to be in degrees (geographic CRS). "
            f"Transform pixel size: {px_w} x {px_h}. "
            f"For accurate area, reproject to an equal-area projection first."
        )
    
    px_area_km2 = (px_w * px_h) / 1e6
    
    if px_area_km2 is None or not np.isfinite(px_area_km2):
        raise ValueError(f"Failed to calculate pixel area.")
    
    # 步骤6: 创建有效像素掩膜（三个数据都有效）
    valid_mask = (
        (lc1_data_reprojected != lc1_nodata) &
        (lc2_data_reprojected != lc2_nodata) &
        (bins_data_reprojected != bins_nodata) &
        (bins_data_reprojected != 99)
    )
    
    # 步骤7: 获取所有土地覆盖类别
    unique_lc1_classes = np.unique(lc1_data_reprojected[valid_mask])
    unique_lc1_classes = unique_lc1_classes[unique_lc1_classes != lc1_nodata]
    
    unique_lc2_classes = np.unique(lc2_data_reprojected[valid_mask])
    unique_lc2_classes = unique_lc2_classes[unique_lc2_classes != lc2_nodata]
    
    # 如果没有提供类别名称映射，使用数值
    if land_cover_class_names is None:
        all_classes = np.unique(np.concatenate([unique_lc1_classes, unique_lc2_classes]))
        land_cover_class_names = {cls: f"Class {cls}" for cls in all_classes}
    
    # 步骤8: 统计每个转换组合在不同bin下的面积
    results = []
    
    # 获取所有有效的bin_id（排除99）
    valid_bins = bins_data_reprojected[valid_mask]
    unique_bins = np.unique(valid_bins)
    unique_bins = unique_bins[unique_bins != bins_nodata]
    unique_bins = unique_bins[unique_bins != 99]
    
    # 对每个初始类别和最终类别的组合进行统计
    for lc1_class in unique_lc1_classes:
        for lc2_class in unique_lc2_classes:
            # 创建该转换组合的掩膜
            conversion_mask = (
                valid_mask &
                (lc1_data_reprojected == lc1_class) &
                (lc2_data_reprojected == lc2_class)
            )
            
            # 如果该转换组合没有像素，跳过
            if not np.any(conversion_mask):
                continue
            
            # 在该转换组合内统计不同bin的面积
            conversion_bins = bins_data_reprojected[conversion_mask]
            
            # 初始化各bin的面积
            bin_areas = {}
            for bin_id in unique_bins:
                bin_mask = (conversion_bins == bin_id)
                bin_count = np.count_nonzero(bin_mask)
                bin_areas[bin_id] = float(bin_count) * float(px_area_km2)
            
            # 获取类别名称
            from_name = land_cover_class_names.get(lc1_class, f"Class {lc1_class}")
            to_name = land_cover_class_names.get(lc2_class, f"Class {lc2_class}")
            
            # 创建结果行（确保列顺序正确）
            result_row = {
                "From": from_name,
                "To": to_name,
            }
            
            # 添加各bin的面积（使用标签名称，按bin_id顺序）
            # 确保所有行都有相同的列
            for bin_id in sorted(unique_bins):
                label = bin_labels.get(bin_id, f"Bin {bin_id}")
                area = round(bin_areas.get(bin_id, 0.0), 1)
                result_row[f"{label} (km²)"] = area
            
            # 只添加至少有一个bin有面积的转换组合
            has_data = any(bin_areas.get(bin_id, 0.0) > 0 for bin_id in unique_bins)
            if has_data:
                results.append(result_row)
    
    # 步骤9: 创建DataFrame
    if len(results) == 0:
        print("Warning: No valid land cover conversions found.")
        return pd.DataFrame()
    
    df = pd.DataFrame(results)
    
    # 确保列的顺序：From, To, 然后是按bin_id排序的列
    column_order = ["From", "To"]
    # 添加所有bin列（按bin_id排序）
    for bin_id in sorted(unique_bins):
        label = bin_labels.get(bin_id, f"Bin {bin_id}")
        column_name = f"{label} (km²)"
        if column_name in df.columns:
            column_order.append(column_name)
    
    # 重新排列列顺序（只保留存在的列）
    existing_columns = [col for col in column_order if col in df.columns]
    df = df[existing_columns]
    
    # 如果没有结果，返回空DataFrame
    if len(df) == 0:
        print("Warning: No valid land cover conversions found.")
        return df
    
    # 步骤10: 保存表格（如果提供了输出目录）
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成输出文件名（基于输入文件名）
        initial_name = initial_land_cover_path.stem
        final_name = final_land_cover_path.stem
        bins_name = productivity_bins.stem
        
        output_file = output_dir / f"land_conversion_{initial_name}_to_{final_name}_{bins_name}.csv"
        df.to_csv(output_file, index=False)
        print(f"Land cover conversion table saved to: {output_file}")
    
    return df
#-------------------------------------------------------------------------------