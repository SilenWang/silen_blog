---
title: Graphviz used with draw.io
categories: Others
date: 2021-10-20 01:48:35
tags: ['graphviz', 'drawio']
---

Since I learned Graphviz, I often use it to create flowcharts or simple diagrams. Most of the time, using it is quite convenient because the flowcharts are small and only used for illustration, not for display. However, when the flowchart is complex and there's a need for both beautification and presentation, it becomes more troublesome.

Therefore, I found an alternative approach: first use Graphviz to create the rough sketch, then use [graphviz2drawio](https://github.com/hbmartin/graphviz2drawio) to convert the written Graphviz into a format that can be imported into [drawio](https://github.com/jgraph/drawio) for manual arrangement.

In practice, this solution isn't particularly convenient either... The flowchart framework is there after import, but the original shapes from Graphviz (like circles and squares) are not retained. Additionally, some objects in Graphviz do not exist in drawio, so using Graphviz to create sketches can only retain basic skeletons; a lot of adjustments are still needed. Therefore, in the future, either learn another language for using drawio or find another similar tool that can import Graphviz and further edit it.
