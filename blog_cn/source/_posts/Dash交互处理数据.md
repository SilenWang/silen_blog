---
title: Dash实现从图形上取数据点进行计算
categories: Script
date: 2022-03-13 03:09:08
tags: ['数据可视化', '交互式分析', 'Dash', 'Plotly', 'Python']
---

之前用dash实现了一个比较有趣的交互数据处理案例, 这里记录一下
<!-- 摘要部分 -->
<!-- more -->
大概的问题是有很多条的数据, 需要基于数据进行指标计算, 但是甲方并不想使用所有数据, 而是需要将数据条目画成一张点图上, 然后从图形上圈出一部分数据点来, 用这些数据点对应的数据条目完成指标计算, 即交互式的选择部分数据来进行指标计算.

我看到这个要求的第一反应是去翻Dash有没有相关功能, 没想到还真有, [见这个官方文档]()

最终写起来也比较简单, 给函数加上对选择事件的反馈就行了:

```python
# 绘制图形
mport plotly.express as px
data_b1_p1 = pd.read_csv('1st.part1.fmt.csv')

fig_b1_p1 = px.scatter(data_b1_p1, x='FFPE_X', y='FFPE_Y', color = 'Pos', opacity=0.5)
fig_b1_p1 = fig_b1_p1.update_layout(
    width=600,
    height=600
).update_traces(
    mode='markers', marker_size=2
)

# 设定dash app组建的框架结构, 我这里只放一个div用来放前面的图
app.layout = html.Div([
    html.Div([
        dcc.Markdown("demo titile"),
        dcc.Graph(
            id='fig_b1_p1',
            figure=fig_b1_p1
        ),
        html.Pre(id='sel_b1_p1')
    ],
])

# 定义数据处理函数
def show_mean(selectedData, base):
    if selectedData:
        selData = pd.DataFrame(selectedData['points']).rename(columns={'x': 'FFPE_X','y': 'FFPE_Y'})
        selData = selData.merge(base, how='left', on=['FFPE_X', 'FFPE_Y']).drop(
            columns=['curveNumber','pointNumber','pointIndex','FFPE_X','FFPE_Y', 'Pos', 'OID', 'Marker']
        )

        mean = nanmean(selData['col1'])
        return f"{mean}"
    else:
        return '等待选择'

# 设定响应函数
@app.callback(
    Output('sel_b1_p1', 'children'),
    [Input('fig_b1_p1', 'selectedData')]
)
def parse1(selectedData):
    return show_mean(selectedData, data_b1_p1)

```

具体来说就是, 用plotly绘制的交互图形, 本来就有region select和lasso select两个功能, 这两个功能选取数据点后, 图形数据会直接存到`selectedData`这个属性的值中去, 用内置的装饰器可让函数得到这部分数据.

不过需要注意的是, 返回得到的只有用来绘制图形的数据, 比如我这里是点图, 所以通过`selectedData`能拿到的就是x, y坐标, 且数据以列表套字典的形式存储, 我需要将其数据框化后, 将原始数据根据坐标再次合并进来(向左合并), 这样就知道选择的点对应哪些数据条目, 之后计算出指标返回就好了.

原数据不好写示例, 之后找个别的数据写了例子后再更新本blog

以上~
