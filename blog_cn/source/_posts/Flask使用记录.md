---
title: Flask使用记录
categories: Script
date: 2020-08-09 01:25:45
tags: ['Flask', 'Python', 'Web开发', '后端', 'API']
---

因项目需要, 最近终于是实际看了下Flask怎么用, 照例来个记录.
<!-- 摘要部分 -->
<!-- more -->

由于我这次的任务比较简单, 只要用Flask做个服务器, 能接受请求然后调用Python代码计算结果就行, 因此我这次只涉及最基本的使用. 当然由于我本身并不熟HTML和js, 因此这里记录的东西纯属我自己的理解...可能包含大量不准确和有误导性的描述...

## 创建一个Flask服务(App)

所有的操作从创建个Flask App开始, 创建也比较简单, 把`Flask`对象实例化就行了.

```python
from flask import Flask
app = Flask(__name__)

if __name__ == '__main__':
    app.run(host='0.0.0.0')
```

上述代码保存为Python文件, 然后`python app.py`就能启动服务了. `host`参数设到`0.0.0.0`是为了所有ip都能访问这个服务.

## 将访问与函数绑定

Flask可以将对页面的请求转换成对Python函数的调用. 具体在对应函数前用`@app.route`装饰器就行. `'/'`的话对应的是`index.html`, 如果是`'/login'`这种的话则对应`login.html`, 当然我实测下来对应的html文件不需要存在, 在返回值里直接写html字串就行.

我这里是准备好了一个说明API使用方法的页面, 读进来直接返回就好.

```python
@app.route('/')
def index():
    with open('index.html') as f:
        content = f.read()
    return content
```

## 处理POST请求

我这次的程序需要接收一些参数和上传的文件内容以进行进一步处理, 文件的内容自然不能走`GET`, 所以只处理`POST`.

Flask内置的`request`对象可以方便的获取请求中带有的所有信息.

```python
from io import BytesIO
from flask import request
from iNeo_MS import sheet_parse

@app.route('/parse', methods=['POST'])
def parse():
    args = request.args
    method = args['method']
    byteData = request.files.get('file')
    data = sheet_parse(BytesIO(byteData.read()), method)
    return str(data.to_dict('records'))
```

## 用Python进行API测试

进行`POST`服务测试是不能直接靠浏览器的, 得安装专门的插件... 然而家里网络堪忧... 于是还是靠Python了.

具体的方法当然是用`requests`模块的`post`方法了, 将参数放到`params`参数中, 文件则读入内容后通过`files`参数送出.

```python
import requests

params = {
    'method': 'pNovo',
}

with open("data/results.res", "rb") as f:
    content = f.read()

files = {
    "file" :('result.res', content),
}

response = requests.post(
    'http://127.0.0.1:5000/parse',
    params=params,
    files=files
```
