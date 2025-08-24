---
title: Dash Implementation for Interactive Data Point Selection and Calculation
categories: Coding
date: 2022-03-13 03:09:08
tags: ['ploty', 'dash', 'python']
---

A while ago, I implemented an interesting interactive data processing case using Dash. Here, I'll record it.

<!-- Abstract part -->
<!-- more -->

The problem was that there were many data entries, and the client wanted to perform metric calculations based on these data entries. However, they didn't want to use all the data; instead, they needed to draw the data points on a scatter plot and then select a subset of these points from the graph to perform metric calculations using the corresponding data entries.

My first reaction was to check if Dash had any relevant features. To my surprise, it does! [See this official documentation](https://plotly.com/python/interactive-graphing-with-dash/)

The implementation turned out to be quite simple. I just needed to add feedback for selection events:

```python
# Drawing the graph
import plotly.express as px
data_b1_p1 = pd.read_csv('1st.part1.fmt.csv')

fig_b1_p1 = px.scatter(data_b1_p1, x='FFPE_X', y='FFPE_Y', color='Pos', opacity=0.5)
fig_b1_p1 = fig_b1_p1.update_layout(
    width=600,
    height=600
).update_traces(
    mode='markers', marker_size=2
)

# Setting up the framework structure for Dash app components; here, I only place a div to hold the graph
app.layout = html.Div([
    html.Div([
        dcc.Markdown("demo title"),
        dcc.Graph(
            id='fig_b1_p1',
            figure=fig_b1_p1
        ),
        html.Pre(id='sel_b1_p1')
    ],
])

# Defining the data processing function
def show_mean(selectedData, base):
    if selectedData:
        selData = pd.DataFrame(selectedData['points']).rename(columns={'x': 'FFPE_X', 'y': 'FFPE_Y'})
        selData = selData.merge(base, how='left', on=['FFPE_X', 'FFPE_Y']).drop(
            columns=['curveNumber', 'pointNumber', 'pointIndex', 'FFPE_X', 'FFPE_Y', 'Pos', 'OID', 'Marker']
        )

        mean = nanmean(selData['col1'])
        return f"{mean}"
    else:
        return 'Waiting for selection'

# Defining the callback function
@app.callback(
    Output('sel_b1_p1', 'children'),
    [Input('fig_b1_p1', 'selectedData')]
)
def parse1(selectedData):
    return show_mean(selectedData, data_b1_p1)

```

In summary, using a Plotly interactive graph, which has built-in region select and lasso select functionalities, the selected data points are stored in the `selectedData` attribute of the graph. Using Dash's decorator, you can pass this data to your function.

However, it's important to note that the returned data is only for drawing the graph; for example, if it's a scatter plot, you'll get the x and y coordinates in list-of-dictionaries format. You need to convert these into a DataFrame and merge them with the original data based on the coordinates to know which data entries correspond to the selected points, after which you can perform the calculations and return the results.

The original data is too large to write here; I'll update this blog with an example using another dataset later.
```