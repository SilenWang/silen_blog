---
title: 又踩坑了--plotly的scattermap不能对大部分marker图标设置颜色
categories: Coding
date: 2025-12-02 15:20:39
tags: ['taipy', 'scattermap']
---

啊, 没想到plotly这么多年了, 还是稍微一用就能触到它的能力边界, 之前是画不了时间轴图, 这次又发现地图的Marker自定义程度不够...

<!-- more -->

我需要用plotly带的`Scattermap`对象绘制一个地图, 地图上需要设置一些标点. 然后需求方看了第一版, 觉得想用地图上常见的反水滴形图标来代替默认的圆形. 正好`Scattermap`支持, 于是我这么写:

```python
map_fig = (
    go.Figure()
    .add_trace(
        go.Scattermap(
            lat = Lat,
            lon = Lon,
            mode = 'markers',
            marker = dict(
                size = 12,
                color = AIMINGMED_COLOR['AIMING_RED'],
                opacity = 0.5,
                symbol = 'marker',
                allowoverlap = True
            ),
            hoverinfo = 'none'
        )
    )
)
```

结果发现, 图标虽然能自定义, 但是并不能变色... 于是又是一番查找, 最后又发现了一个去年的[Issue](https://github.com/plotly/plotly.py/issues/5153), 这里面给了很绝望的答案, `plotly`到目前为止就是不支持这个改变, 因为它调用的是`maki-icons`, 这里面的图标是固定黑色的, 也没有提供通过请求改颜色的接口. 

这个问题的等级是P3, 也就是... 排上日程, 但是最不重要的那一级别... 又是个猴年马月的东西了...

但是我可不能等, 我又折磨了Copilot一上午, 最后得到的了一个更绝望的答案, `plotly.py`只是`plotly.js`的一个包装, 只负责将数据和参数按照要求整理好, 实际的绘制完全由`plotly.js`实现, 也就是说, 如果我要修复这个问题, 我又要看JS代码, 并且不是修改了就好, 我还要把JS部分的接口写出来, 然后再跟Python的部分对接上... 这... 不是不行, 但是我时间不够...

于是我又一排脑门想邪门办法了, 根据观察, `plotly.js`实际上会通过https请求图标的svg文件, 然后用这个文件进行绘图, 那... 有没有可能, 我靠油猴脚本, 拦截这个请求, 塞我准备好的svg文件给他呢?

这次是chatGPT救我狗命了, deepseek完败... chatGPT确实给我了我一个能够拦截并替换svg文件的脚本:

```javascript
// ==UserScript==
// @name         Intercept Mapbox Maki Icon via Image.src
// @namespace    http://tampermonkey.net/
// @version      1.0
// @description  拦截 Mapbox 图标请求并替换成自定义 SVG
// @match        *://*/*
// @run-at       document-start
// @grant        none
// ==/UserScript==

(function() {
    'use strict';

    const TARGET_PREFIX = "https://unpkg.com/maki@2.1.0/icons/";
    const CUSTOM_SVG_CONTENT = `
<svg version="1.1"
	 id="svg4619" inkscape:version="0.91 r13725" sodipodi:docname="marker-15.svg" xmlns:cc="http://creativecommons.org/ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" xmlns:svg="http://www.w3.org/2000/svg"
	 xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" width="15px" height="15px"
	 viewBox="0 0 15 15" style="enable-background:new 0 0 15 15;" xml:space="preserve">
<path id="path4133" inkscape:connector-curvature="0" d="M7.5,0C5.0676,0,2.2297,1.4865,2.2297,5.2703
	C2.2297,7.8378,6.2838,13.5135,7.5,15c1.0811-1.4865,5.2703-7.027,5.2703-9.7297C12.7703,1.4865,9.9324,0,7.5,0z" fill="#E72410"/>
</svg>
`;
    const CUSTOM_SVG_URL =
        "data:image/svg+xml;base64," + btoa(CUSTOM_SVG_CONTENT);

    // 劫持 Image.src
    Object.defineProperty(Image.prototype, "src", {
        set(url) {
            if (typeof url === "string" && url.startsWith(TARGET_PREFIX)) {
                console.log("[Tampermonkey] 拦截 Mapbox 图标: ", url);
                url = CUSTOM_SVG_URL;
            }
            this.setAttribute("src", url);
        },
        get() {
            return this.getAttribute("src");
        }
    });

})();
```

总算是... 先实现了, 之后的问题, 之后再说吧哎, 继续加班了...