---
title: 基于dash与flask搭建一个任务监控web-app
categories: Bioinfomatic
date: 2021-06-21 01:24:23
tags: ['python', 'flask', 'dash', 'mongoDB', 'dash', 'pymongo', 'flasgger']
---

鸽了都快大半年了, 自从接锅以来, 事情一直都多得做不完... 也就没那个兴致写东西, 不过这一次算是将半年来新学的一些东西做了综合应用, 写了我人生中第一个带数据库操作的Web-App, 所以还是应该纪念一下, 不过时间有限, 先写文本, 图片以后再补...
<!-- 摘要部分 -->
<!-- more -->

这个app实现的是类似NCBI Blast那种分析服务, 目的是让用户能通过Web界面提交分析任务, 并从监控列表中了解已提交任务的分状态, app提供了一个极其简陋的分析界面, 用户可在界面选择必要的输入文件和配置文件(事先登录到数据库)后, 开始分析任务, 并通过监控栏了解已提交任务的运行状态. 当然这个Web见面其实不是给用户用的, 而是作为与前端IT交流的Demo, 并在日后供我自己测试使用, 真正派上用场的是内建的Web-API, 应用通过响应前端发来的Web请求来完成实际的业务工作

## 模块介绍
本次的Web-App功能虽然简陋, 但是实现起来其实用了不少东西, 同时也花了不少时间在一些内容的组合实现上, 本次只着重于实际实现过程中一些问题的解决方式, 一些新模块的基本操作(dash, python操作mongoDB之类)之后可能会单独再水一次

### Web API
API基于Flask构建, 之所以不采用Django是因为介绍里说Flask更轻量简洁, 毕竟我们也不是搭什么巨型高并发高计算量的API, 我觉得用这个更简单

### Web操作界面
Web界面基于Dash构建, Dash与R语言中的Shiny是完全同性质的东西, 只不过这东西在中文生信群体中似乎完全没有存在感, 中文材料基本都是英文机翻. 当然英文材料似乎也不如Shiny来的多, 主要还是以官方的文档和官方的Q&A为主. 不过我既然现在主语言是Python, 使用基于Python的东西肯定比和R混编简单得多. 另外虽然我的App很简单, 但是相对过去写的生信小脚本那还是复杂了好几个量级的... 因此Python的错误追踪还是及其救命的, 不至于让我在Debug上浪费巨量的时间

另外值得一提的是, Dash是从Plotly这个绘图包发展而来的, 同时Dash提供的Web-App也是在Flask基础上实现的, 所以Dash界面+Flask处理Web请求在我当前的应用场景下极佳的组合, 因为我可以在一个主程序文件内同时写Dash函数和Flask函数, 减少重复工作, 同时也保证我测试展现的东西和同事通过API获取的东西是一致的

### 后端任务控制
任务控制基于Snakemake进行, 这个就不多说了, 值得一提的是用了这么长时间之后, 我对Snakemake文档里面的那句, Snakemake is python code这句话有了更好的理解, 因为这个特性的存在, 我能在Snakemake里面通过自编写的一些函数实现Snakemake当前没有实现的特性. 比如在这个应用中的任务状态自动更新

### 文件记录与任务记录(数据库)
任务和文件记录使用了mongoDB进行数据存储, 其实如果单单只是要构建这个App的话, 并不需要用到任何数据库, 使用简单的文件做存储, 或者直接将内容存储在设定好结构的目录里也可以, 使用数据库纯粹是为了进行简单的数据库操作练习, 因为后面有别的工作会用.

而选择mongoDB, 不使用之前已经用过的SQLite, 也是出于后面工作的考虑, 因此提前拿简单的内容进行mongoDB练手

### Api文档
Api文档这个其实最开始并不在应用规划内容里面, 是因为近期搭建了非常多简单的Web服务以实现公司内部的分析上系统, Api一多了, 一个Api准备一份文档存和查询起来相当不方便, 同时因为预计不同的Api会交给不同人, 另外维护一份总的文档也比较烦, 因此我想是否有方式能把Api文档直接生成在单个Web-App下面, 最后就选定了`flasgger`, 这个模块的作用是在flask内引入一系列API, 方便作者直接在Python代码内编写Api文档, 编写的文档可以调用Swagger UI直接在应用路径下形成一份Api文档

## 应用总览
本App其实逻辑很简单, 就是前端接受请求, 将请求转发到后端启动分析进程, 分析进程在运行过程中不断将状态更新到数据库, 然后前端定时或者根据请求从数据库获取程序运行状态就好, 大概的信息流见下图

![这里咕了一张图]()

## 实现要点
虽然这个应用简单, 但是实际实施过程中还是有问题是要额外想办法解决的:
1. 保证Dash-App处理和Flask-API处理的一致性, 即使用Dash-App完成测试后不需要额外对Flask-API部分进行改动
2. Dash-App/Flask-App/API文档的混用, 即让所有这些东西都在同一个App下, 而不是搭建多个Web-App然后相互跳转
3. 利用Snakemake控制进程本身将运行状态更新到数据库, 而不另外设立单独的监控进程

### Dash-App和Flask-API处理的一致性
这一点其实很好实现了... 就是把实际的业务处理函数跟Dash和Flask的函数分开就好, 也就是Dash和Flask的响应函数将数据处理成统一的格式, 然后交由一个实际的业务函数来执行数据处理, 两个响应函数得到处理的数据后再分别格式化为Dash需要的内容或API中设定好的格式:

```python
@app.callback(
                [
                    Output('info-div', 'children')
                ],
                [
                    Input('upload-info', 'contents')
                ],
                [
                    State('upload-info', 'filename')
                ]
            )
def step2_dash(content, filename):
    '''
    dash响应函数
    '''
    # 需要从文件解析目标路径
    if content:
        _, content_stream = content.split(',')
        info_data = step2(content_stream, filename)
        json_fmt = dumps(info_data, sort_keys=True, indent=4, separators=(",", ":"), ensure_ascii=False)
        return [html.Pre(f'解析成功, 从信息表中提取的信息:\n{json_fmt}')]
    else:
        return [html.Pre('无解析结果')]


@app.server.route('/step2', methods=['POST'])
def step2_flask():
    '''
    flask响应函数
    '''
     data = loads(flask.request.form.get('data'))
     fileUrl = f'{data["dowloadFileUrl"]}?type={data["type"]}&fileName={data["fileName"]}'
     filename = re.findall(r'fileName=(.+)', fileUrl)[0]
     r = requests.get(fileUrl)
     info_data = step2(r.content, filename, decode=False)
     res_data = {'errCode':0, 'data': info_data, 'msg': ''}
    json_fmt = dumps(res_data, sort_keys=True, ensure_ascii=False)
    return json_fmt.encode(encoding="UTF-8")


def step2(content_stream, filename, decode=True):
    '''
    step2的实际处理部分, 拆分出来分别接收不同的请求
    '''
    if decode:
        decoded = base64.b64decode(content_stream)
    else:
        decoded = content_stream
    info_data = data_format(decoded, filename)
    return info_data
```

### Dash-App/Flask-App/API文档混用
Flask和Api文档的混用其实比较简单, 因为flasgger本来就是个flask加API文档的. 而Dash和Flask的混用其实也不难, Dash对象有一个`server`属性, 这里面存的就是Flask App对象, 所以其实需要用到该对象的部分直接`app.server`就好了:

```python
import dash

app = dash.Dash("Demo")

app.server.config['SWAGGER'] = {
    'title': '质谱信息录入/报告生成服务文档',
    'uiversion': 2
}
swagger = Swagger(app.server)

app.layout = html.Div(
    [
        html.H1('Demo'),
        html.Div([
            dcc.Upload(
                id='upload-info',
                children=html.Div([
                    'upload here',
                ]),
            ),
            html.Div([
                html.Pre(''),
            ], id='info-div')
        ])
    ],
)

@app.callback(
                [
                    Output('info-div', 'children')
                ],
                [
                    Input('upload-info', 'contents')
                ],
                [
                    State('upload-info', 'filename')
                ]
            )
def step2_dash(content, filename):
    '''
    dash响应函数
    '''
    return [html.Pre('parsed')]


@app.server.route('/step2', methods=['POST'])
def step2_flask():
    '''
    flask响应函数
    ---
    tags:
        - Step2
    parameters:
        - name: dowloadFileUrl
          in: query
          type: string
          required: true
          description: 文件下载地址
    responses:
        200:
            description: 返回解析完成的数据
            schema:
                properties:
                    errorCode:
                        type: number
                        description: 状态码
                        default: 0
                    data:
                        type: string
                        description: 投递任务的key
                        default: parsed
                    msg:
                        type: string
                        description: 正常时为空字串, 状态码不为0时会附带错误信息
                        default: ""
    '''
    return 'parsed'
```


### 利用Snakemake控制进程本身更新运行状态

最开始我是考虑用Snakemake自带的`onstart`和`onerror`特性来实现的, 毕竟我要的其实就是任务开始和结束的时候向数据库提交最新的运行信息.
但实际操作以后发现这两个特性是针对整个Snakefile起效的, 也就是总任务开始前执行一次, 然后总任务失败的时候执行一次, 并不提供单个任务的响应, 那就只能自己动手了...

这时候就要利用Snakemake is python code的特性了. 在一条规则(rule)内, 指定`shell`关键字时其实是在调用snakemake内部Api中的`shell`函数执行系统级的命令, 所以其实写`run: shell('BASH_SCRIPT')`的方式是一样的(程序表现上还有一点额外差异, 这里忽略), 但是改用这种方式之后, 我就可以在系统命令的前后来额外增加其他的python代码来实现我要的功能了, 比如执行前后都更新一次任务状态:

```smk
rule step1:
    input:
        done = '{rid}/analysis.done'
    output:
        txt = '{rid}/done.txt',
    run:
        # update_status是更新状态的函数, rid是任务编号
        update_status(rid, 'step1', 'running')
        shell(
            '''
            cp {input.done} {output.txt}
            '''
        )
        update_status(rid, 'step1', 'done')
```

如此就方便的实现了目的了. 另外, 对于错误时的处理不需要这么进行, 因为错误状态在本应用中不需要精确到规则, 利用`onerror`就可以了

不过还需要注意一点, 就是Sanakemake其实支持一项名为papernote的服务, 该服务可专门收集主控流程的任务状态的, 也就是说Snakemake其实内建有类似的东西, 不过这就需要去读一读相关部分的代码来了解怎么实现的了, 也许内部API有更好用的方式呢?

## 后续
去年的时候我面试过一个生信工程师, 他在疫情期间回不了公司, 就自己在家自学代码然后写了个分析提交/监控服务. 现在我总算是也会了... 回看一下的话, 其实真的不难, 只是需要时间和一定的外部动力推着我学... 另外就是实践还是重要的, 只有实践才能发现实际的问题, 然后解决了问题, 能力也就随之提高了, 我以后还是要多加油...
