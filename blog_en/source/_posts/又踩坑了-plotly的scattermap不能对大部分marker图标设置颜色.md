---
title: Another Pitfall, Plotly's Scattermap Cannot Set Colors for Most Marker Icons
categories: Coding
date: 2025-12-02 15:20:39
tags: ['taipy', 'scattermap', 'Tampermonkey', 'script']
---

Ah, I didn't expect that after so many years, Plotly still hits its limits with just a little use. Previously it couldn't draw timeline charts, and this time I found that the customization of map markers is insufficient...

<!-- more -->

I needed to use Plotly's `Scattermap` object to draw a map with some markers on it. After seeing the first version, the client wanted to use the common inverted‑teardrop icon (like a map pin) instead of the default circle. Fortunately `Scattermap` supports that, so I wrote:

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

It turned out that although the icon could be customized, its color couldn't be changed... After some searching, I found a [last‑year's issue](https://github.com/plotly/plotly.py/issues/5153) that gave a rather hopeless answer: up to now, `plotly` simply does not support this change because it uses `maki‑icons`, where the icons are fixed black and there is no interface to change the color via requests.

The issue is labeled P3, meaning… it's on the schedule, but at the lowest priority… something that will be done who‑knows‑when…

But I couldn't wait. I tortured Copilot for another morning and got an even more desperate answer: `plotly.py` is just a wrapper for `plotly.js`; it only organizes the data and parameters as required, while the actual drawing is completely handled by `plotly.js`. In other words, if I wanted to fix this problem, I would have to look at the JS code, and not just modify it—I would also need to expose the interface on the JS side and then connect it with the Python part… That… isn't impossible, but I didn't have enough time…

So I slapped my forehead and thought of a hack. From observation, `plotly.js` actually fetches the icon's SVG file via an HTTPS request and uses that file for drawing. Well… is it possible to intercept that request with a userscript (like Tampermonkey) and feed it my prepared SVG file?

This time ChatGPT saved my ass, Deepseek completely failed… ChatGPT indeed gave me a script that could intercept and replace the SVG file:

```javascript
// ==UserScript==
// @name         Intercept Mapbox Maki Icon via Image.src
// @namespace    http://tampermonkey.net/
// @version      1.0
// @description  Intercept Mapbox icon requests and replace with custom SVG
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

    // Hijack Image.src
    Object.defineProperty(Image.prototype, "src", {
        set(url) {
            if (typeof url === "string" && url.startsWith(TARGET_PREFIX)) {
                console.log("[Tampermonkey] Intercepted Mapbox icon: ", url);
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

Finally… I got it working for now. The remaining problems can be dealt with later… Sigh, back to overtime work…
