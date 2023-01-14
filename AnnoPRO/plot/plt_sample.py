import numpy as np
import pandas as pd
from typing import Dict, Tuple, List, Optional
from highcharts import Highchart

colormap= {'Composition': '#65B2CD',
            'Autocorrelation': '#E7A11C',
            'Interaction': '#FF7C00',
            'Physiochemical': '#E898A3',
            'Quasi-sequence-order descriptors': '#2EAC38',
            'PAAC for amino acid index set': '#ECC745',
            'Amphiphilic Pseudo amino acid composition': '#CE86B6'}

def getHeatColor(intensity: float,
                 topColor: str,
                 bottomColor: str = "#ffffff") -> str:
    """
    根据强度值计算热力图颜色值

    Args:
        intensity 强度值,0~1
        topColor intensity为1时的颜色,为十六进制颜色字符串
        bottomColor intensity为0时的颜色,默认为白色

    Returns:
        与强度对应的十六进制颜色字符串
    """
    assert intensity >= 0 and intensity <= 1 + 1e-8, f"Intensity should range from 0 to 1, But got {intensity}"
    bRGBColor = hexToRGB(bottomColor)
    tRGBColor = hexToRGB(topColor)
    color = ((1 - intensity) * bRGBColor[i] + intensity * tRGBColor[i]
             for i in range(3))
    return rgbToHex(tuple(color))


def hexToRGB(hexColor: str) -> Tuple[int]:
    """
    将十六进制颜色字符串转为RGB三色整数元组

    Args:
        hexColor 十六进制颜色字符串

    Returns:
        RGB三色的元组
    """
    assert hexColor.startswith("#"), "Ilegal hex color format"
    try:
        d = int(hexColor.replace("#", ""), base=16)
        b, d = d % 256, d // 256
        g, d = d % 256, d // 256
        r, d = d % 256, d // 256
        return r, g, b
    except Exception:
        raise ValueError("Ilegal hex color format")


def rgbToHex(color: Tuple[int]) -> str:
    """
    将RGB三色整数元组转为十六进制颜色字符串

    Args:
        color RGB三色整数元组

    Returns:
        十六进制颜色字符串
    """
    assert color is not None and len(color) == 3, "Ilegal rgb color format"
    try:
        r = int(color[0] % 256)
        g = int(color[1] % 256)
        b = int(color[2] % 256)
        return "#%02x%02x%02x" % (r, g, b)
    except Exception:
        raise ValueError("Ilegal rgb color format")



def plotSampleGrid(sampleData: np.ndarray,
                   gridData: pd.DataFrame,
                   title: str,
                   subtitle: str,
                   filename: str,
                   shape: Tuple[int]=(39,39),
                   colorMap: Dict[str, str]=colormap,
                   width: int = 1000,
                   height: int = 850):
    """
    绘制样本的ProMap单通道特征图

    Args:
        sampleData  HW单通道的样本数据矩阵
        gridData    包含[x,y,label,subtype]四列的数据框,x指的横向坐标,y指的是纵向坐标
        colorMap    每一个subtype对应的颜色
        shape       网格的形状
        title       标题
        subtitle    子标题
        filename    输出文件名
        width       画布宽度
        height      画布高度

    Returns:
        HighChart对象
    """
    gridData = gridData.copy()
    gridData["v"] = gridData.apply(lambda row: sampleData[row.y][row.x],
                                   axis=1)
    gridData = gridData[gridData.v != 0]
    H = Highchart(width=width, height=height)
    H.set_options('chart', {'type': 'heatmap', 'zoomType': 'xy'})
    H.set_options('title', {'text': title})
    H.set_options('subtitle', {'text': subtitle})
    H.set_options(
        'xAxis', {
            'title': None,
            'min': 0,
            'max': shape[1],
            'startOnTick': False,
            'endOnTick': False,
            'allowDecimals': False,
            'labels': {
                'style': {
                    'fontSize': 20
                }
            }
        })
    H.set_options(
        'yAxis', {
            'title': {
                'text': ' ',
                'style': {
                    'fontSize': 20
                }
            },
            'startOnTick': False,
            'endOnTick': False,
            'gridLineWidth': 0,
            'reversed': True,
            'min': 0,
            'max': shape[0],
            'allowDecimals': False,
            'labels': {
                'style': {
                    'fontSize': 20
                }
            }
        })
    H.set_options(
        'legend', {
            'align': 'right',
            'layout': 'vertical',
            'margin': 1,
            'verticalAlign': 'top',
            'y': 60,
            'symbolHeight': 12,
            'floating': False,
        })
    H.set_options(
        'tooltip', {
            'headerFormat': '<b>{series.name}</b><br>',
            'pointFormat': '{point.label} {point.v}'
        })
    H.set_options('plotOptions', {'series': {'turboThreshold': 50000}})

    for subtype, color in colorMap.items():
        dfi = gridData[gridData['subtype'].apply(
            lambda subtypes: subtype in subtypes.split("||"))].copy()
        if len(dfi) == 0:
            continue
        dfi["color"] = dfi["v"].apply(
            lambda intensity: getHeatColor(intensity, colorMap[subtype]))
        data = dfi.to_dict('records')
        H.add_data_set(data, 'heatmap', name=subtype, color=color)
    if filename:
        H.save_file(filename)
    return H