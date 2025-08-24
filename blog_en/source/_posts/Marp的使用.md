---
title: 'Using Marp'
categories: Others
date: 2022-06-22 15:40:28
tags: ['Marp', 'Python']
---

In daily work, there are always times when you need to create a presentation document to show project progress or results. However, creating a PPT is quite time-consuming (my obsessive-compulsive disorder always makes me adjust the font and position). I've never used any of the fancy elements and features of PPTs, so I found Marp to save some effort.


<!-- Summary -->
<!-- more -->


Marp is a project that generates presentation documents using markdown syntax. With basic markdown syntax, you can quickly create a presentable document (it's actually quite ugly), supporting generation into pptx, pdf, and web pages. I've been using it for about a month, and the experience has been good (the audience might find it eye-straining).

## Basic Usage of Marp in VScode
The Marp project has multiple available programs, but the most convenient one is definitely the VScode plugin. No special settings are required; just install and use it.


![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_20220622_1655884337.jpg)

After installation, when you open a markdown file, an additional button will appear in the upper right corner of the editor. You can switch the current file to Marp mode.

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Set_20220622_1655884577.jpg)


Then you can start writing... Since it's all basic markdown syntax, there's not much to pay attention to except for using `---` to separate the content of each page. For more details, see the [official documentation](https://marpit.marp.app/markdown).

## Use Cases for Reference

### Inserting Images Left/Right

In my presentations, I often need to include an image and then provide some explanations. This side-by-side arrangement looks more reasonable, but basic markdown syntax doesn't support this... You can achieve it using Marp's background image syntax.

Here, `contain` automatically adjusts the size of the image; otherwise, it usually exceeds the slide range. The example uses right placement, but you can also use `left` to place it on the left side. For vertical alignment, `contain` doesn't seem to work.

```markdown
---
marp: true
---

## A statistics chart is on the right
- The image is randomly found

![bg contain right](https://seaborn.pydata.org/_images/grouped_barplot.png)
```

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Slide1_20220622_1655885449.jpg)

### Footnotes and Page Numbers

By adding two settings at the beginning, you can enable footnotes and page numbers. However, the position seems to be unchangeable using markdown syntax.

```markdown
---
marp: true
footer: 20220610

---
## Displaying page numbers
```

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Slide2_20220622_1655885767.jpg)
