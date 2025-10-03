---
title: 使用taipy开发一个数据面板以及遇到的坑
categories: Coding
date: 2025-10-01 22:49:37
tags: ['Taipy', 'Python']
---

Python下快速开发数据或AI相关应用的模块真的很多，我已经使用过的包括Dash、Streamlit、Gradio、NiceGUI，前段时间我甚至又发现了俩，正遇上我需要开发一个展示公司数据的简单数据看板，于是我再次不要命的用了新框架----Taipy。

<!-- more -->

## Taipy的使用

### 页面布局

Taipy设计了多种编写布局的方式，根据个人熟悉的语言，可以从Markdown、HTML以及Python语法中挑选一种，进行页面布局设定。不过感觉官方主推的是Markdown。是的，非常神奇，一个用Python写的库居然主推用Markdown写布局...

我当然是选Python语法了，它的用法跟Gradio基本一致，`with`以及对象用来控制元素的包含关系。

具体来说，通过`taipy.gui.builder`模块，我们可以使用Python语法来定义页面结构。

主要布局组件包括：
- `tgb.layout()`：创建布局容器，可以设置列数（columns）、间距（gap）
- `tgb.part()`：创建内容区块，用于组织相关组件

在我的数据面板中，我采用了三级布局结构：
1. 顶部：公司Logo、标题和主题切换控件
2. 中间：地图、饼图、折线图、柱状图等数据可视化组件
3. 底部：关键指标卡片展示

布局示例：
```python
with tgb.layout(columns="1 2 1", gap="10px"): # 三列，宽度是1:2:1
    # 左侧区块
    with tgb.part(class_name="card"):
        tgb.text("里程碑", class_name="h5")
        tgb.chart(figure="{timeline}", plot_config={"staticPlot": True})
    # 中间和右侧区块类似
```

### 样式设定

Taipy支持通过CSS文件和Stylekit来自定义样式。

- **CSS样式**：CSS的使用跟一般写网页一样，创建了`style.css`文件来定义全局样式和组件样式：
```css
.card {
    --element-padding: 0.4rem;
}

.scrollable-part {
    overflow-y: auto; /* 垂直滚动 */
}
```

然后在项目创建时，指定要加载的额外CSS就行了。

```python
gui = Gui(page, css_file="style.css")
gui.run()
```


- **Stylekit主题**：通过Python字典定义主题颜色，然后在启动GUI时传入：
```python
STYLEKIT = {
    "color-primary": "#E72410",
    "color-secondary": "#F28E2A",
    "color-background-light": "#EDEDED",
    "color-paper-light": "#FFFFFF",
}

gui.run(stylekit=STYLEKIT, dark_mode=False)
```
注意，Stylekit能控制的内容远少于经典的设置class然后CSS文件内调整，可调整的内容见[Taipy的文档](https://docs.taipy.io/en/latest/userman/gui/styling/)。

除上述两种方式之外，Taipy的不同控件还可以通过`class_name`参数来设定一些模块预设的样式，每个控件能用的东西不同，具体需要参考各个控件的文档


### 数值绑定

Taipy提供一种名为[数值绑定](https://docs.taipy.io/en/latest/userman/gui/binding/)的方式，来让界面内容的更新更为简单。

具体使用花括号`{}`语法将Python变量绑定到界面组件。举最简单的文字例子，我可以如下设置显示一段文字：

```python
text = 'My Text'
tgb.text("{text}")
```

在上面的设置中，先复制变量`text`，再在文字控件的内容部分用`"{text}"`将变量于内容绑定看上去有些多此一举，但是这在Taipy中非常重要。Taipy中，会将命名控件中定义过的所有对象/变量，放在一个名为`state`的特殊对象中。在上面的代码中，实际上同时创造了`state.text`这个属性，而任何对`state.text`进行的改变，也会反映在`tgb.text("{text}")`这个文字控件上。于是，如果我想更新这里显示的文字，只要在能获取state的部分中，直接修改其内容，文字就会被修改，并显示在页面上了，比如下面这样：

```python
state.assign(text, "Changed Text")
```

这样相较于常用的Dash/Steamlit/Gradio，其实就少写了更新控件内容部分的代码，使开发更快。尤其对于数据展示的App，很多时候其实图和表的绘制代码是固定的，我们只是重新整理数据框或计算指标，Taipy的这种设计能让开发者将时间着重放在更新数据上，而不用在展示部分画太多时间。

### 回调函数

如其他框架一样，Taipy的各种控件也是可以设置回调函数的，比如最基础的按钮有`on_action`回调，就是按钮被按下时触发的回调。不同控件可用的回调不同，需要参考官方文档。

与其他框架不同的是，Taipy还可以设置全局的默认回调函数，我们可以直接定义名为`on_action`的函数，放在命名控件内，Taipy的逻辑是，如果一个控件可以触发`on_action`回调，且有名为`on_action`的函数在当前命名控件，则会直接触发回调。因此在Taipy中，用这种全局回调需要做好逻辑判断，误触回调。


### 后台任务

Taipy虽然内置了定时调度方法，但却放在了付费版本中。不过好在，社区版本中可以通过`invoke_long_callback`这种高负载的任务回调来时间定时刷新数据的目的。我这边参考了[官方给的一个例子](https://github.com/Avaiga/demo-realtime-pollution)


```python
from random import randint
from time import sleep


def countdown():
    while True:
        sleep(10)


def update_value(value):
    value = value + randint(1, 4)
    return value


def update_display_value(state):
    state.total_client = update_value(state.total_client)
    state.total_visits = update_value(state.total_visits)
    state.total_reads = update_value(state.total_reads)


def on_init(state):
    set_icon(state)
    invoke_long_callback(
        state,
        user_function = countdown,
        user_status_function = update_display_value,
        period = 3000
    )
```

上面的代码中，`invoke_long_callback`是Taipy提供的一个函数，我们通过它来设定一项高负载任务，也就是`user_function`，这个函数会单开线程执行，然后设定一个检查高负载任务结果的线程，调用`user_status_function`中传入的函数，来更新state中的数据。这样一旦数据更新，就能触发Taipy的数值绑定逻辑，更新页面上显示的内容了。


## 踩坑记录

原本想，Taipy都有明确的商业化方案了，这个模块的基本功能应该做得相当不错了，但是在实际使用的过程中，还是碰到让我卡了一天以上的坑...

### 4.1.0版本的Toggle控件在Theme模式不触发

在pixi/conda能获取的4.1.0版本中，我发现Toggle控件在设置`theme=True`，即设置为主题转换按钮时，数值绑定和函数回调都是实效的... 这个问题苦恼了我一天半... 我试着让AI给了我很多解决方法，但是Taipy似乎并不是一个那么流行的模块，不能像其他之前的框架那样，直接注入Javascript来解决。因此我最后还是跑去看了源代码... 发现在4.1.0版本中，`ThemeToggle.tsx`中定义的主题切换控件，就没有数值绑定和回调的逻辑... 只有在最新的`main`分支中才有... 最后，靠从源代码安装Taipy解决了这个问题...

### Stylekit并不能一次性设定所有的颜色

官方虽然在文档中说，使用Stylekit这个设计可以让开发者灵活和简洁的控制主题配色，但是实际使用后发现，它似乎并不能控制页面的总背景色，需要通过css额外解决：

```css
/* Setting theme color */
body {
    background-color: var(--color-background);
}
```

### plotly 没有绘制时间线图的函数

这个... 纯属对plotly吐槽以下了，虽然Python在数据科学方面应用甚广，但是不论是Plotly、matplotlib、seaborn还是altair，其图形的丰富程度，还是比不过R下的生态... 当然R下有人写了，但是包管理太混乱，我用不了别人写的代码，又是另一回事了...

不过还好里程碑时间线这种东西也就是点、线加个文字，最后还是靠AI手措出来了...

### 对Plotly的设定方式说明不足

这个问题也卡了我很久，Taipy的文档中，有说明如何不使用Taipy提供的接口，而是直接传递Ploty对象，将图形绘制在页面上。但是如何对Plotly对象进行配置，文档就说参考Plotly文档了。然后我就体会到了Plotly这个文档、看上去写得很细，实际上内容很贫乏的问题了。在他们官方的文档中，很多的配置项，居然只有`fig.show()`的时候才能传入... 而不是直接对Plotly对象去配置... 然后我绕了一大圈，发现Taipy的`chart`控件有一个`plot_config`参数... 这里可以传入`fig.show()`时的参数...