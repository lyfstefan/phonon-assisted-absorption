# 论文绘图模板：当前设置

本模板只规定通用绘图样式，不包含具体论文的物理模型。

- 每个数据坐标框保持正方形，边长约 2.5 cm；整张画布按面板数量和文字留白调整。
- 全部文字、数字和数学符号使用真正的 LaTeX，字体为 Computer Modern。数学变量、下标和比值写在数学模式中，例如 `$U_1$`、`$U_2$`、`$V/t$`；相字母使用数学斜体。
- 轴标签 10 pt，刻度和普通标注 9 pt，边框 1 pt。
- 图例优先放在坐标框内部。紧凑图例默认 7 pt，可按实际拥挤程度调整为 6.5–8 pt；数值保留两位小数，重复的量放在图例标题中，条目可分成两列。
- 图例不能挡住曲线或数据。必要时只增加显示范围上方的留白；不要改动数据值或归一化。
- 黑色观测点画在下层（`OBSERVATION_ZORDER = 2`），彩色拟合线画在上层（`FIT_ZORDER = 3`）。
- 每个热图使用独立的短细 colorbar：高度为数据框的 60%，宽度为框宽的 4.5%；底端与数据框下沿对齐，标签放在 colorbar 上方。可比面板仍使用相同的数值归一化。
- 导出 PDF、SVG 和 300 dpi PNG，并记录字体、方框和 colorbar 的尺寸。保留完整画布，检查文字裁切和图例遮挡。

LaTeX 设置：

```python
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{textcomp}'
```

使用 `inset_legend(ax, title=..., ncol=2, data_top=...)` 添加框内图例；
`data_top` 可为线性纵轴预留顶部空间。`legend_value(value)` 只把显示数字
格式化为两位小数，保存的数值精度不变。`save()` 会保留这类图例的 7 pt 设置；
其他已有图例可通过 `legend_fontsize=` 明确指定字号。

最小用法：

```python
import matplotlib.pyplot as plt
from publication_style import (
    apply_style, save, inset_legend, C,
    OBSERVATION_ZORDER, FIT_ZORDER,
)

apply_style()
fig, ax = plt.subplots(layout='constrained')
ax.scatter(x, observations, color=C['black'], label='Data',
           zorder=OBSERVATION_ZORDER)
ax.plot(x, prediction, color=C['red'], label='Fit', zorder=FIT_ZORDER)
ax.set(xlabel=r'$x$', ylabel=r'$y$')
inset_legend(ax, ncol=1, data_top=max(max(observations), max(prediction)))
save(fig, 'figure')
```
