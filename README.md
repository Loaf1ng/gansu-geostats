# 甘肃地统计实习

甘肃省 53 个台站的降水、蒸发插值与县域林草覆盖率分析。

[课程数据](实习-第一部分.zip) · [ArcMap 工程](gis/gansu.zip) · [模型与预测](results/)

```sh
python -m pip install -r requirements.txt
python analysis/run.py
```

Python 3.14，结果写入 `run/`。查看地图可直接解压工程，打开 `gansu.mxd`。

| 运行选项 | 用途 |
| --- | --- |
| `--figures` | 用 R 重绘 14 幅图，需要 sf、terra 和 UnRAR 或 7-Zip。 |
| `--full` | 重跑空间分块、嵌套验证、变程敏感性分析和台站回归。 |
| `--arcmap "C:\Python27\ArcGIS10.8\python.exe"` | 从原始栅格重算县域统计，需要 Spatial Analyst 和 UnRAR 或 7-Zip。 |
