---
title: 'Marp\'s Usage'
categories: Script
date: 2022-06-22 15:40:28
tags: ['Marp', 'Python']
---

Daily work often requires creating a demonstration document to show project progress or results, but making a PPT is actually quite time-consuming (perfectionism always controls my hands when adjusting fonts and positions), and the content I usually present consists of simple text and charts. I've never used the fancy elements and features that PPT offers, so I found Marp to help me procrastinate...

<!-- Abstract part -->
<!-- more -->

Marp is a project that generates presentation documents using markdown syntax. With basic markdown syntax, you can quickly generate a document that looks decent (although it might look a bit ugly). It supports generating pptx, pdf, and web pages. I've been using it for a month, and the experience has been quite good (the viewers might think my eyes are blind).

## Marp_VScode Basic Usage
Marp offers multiple available programs, but the most convenient one is undoubtedly the VScode plugin. No special settings are required; just install and use it.

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_20220622_1655884337.jpg)

After installation, when you open a markdown file, an additional button will appear in the upper right corner of the editor, allowing you to switch the current file to Marp mode.

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Set_20220622_1655884577.jpg)

You can now start writing... because it's all basic markdown syntax. The only thing to note is that you use `---` to separate the content of each page, and there's not much else to worry about. Check out the [official documentation](https://marpit.marp.app/markdown) if needed.

## Use Case Reference

### Inserting Images Left/Right Align

In my presentations, I often need to place a chart next to some text for explanation purposes. This layout is more reasonable than having the image on one side and the text on the other. However, markdown's basic syntax doesn't support this directly. You can achieve it using Marp's background image syntax.

Here, `contain` makes the image automatically adjust its size; otherwise, it would usually exceed the slide range. The example uses right alignment; you can also use `left` for left alignment. For top and bottom positions, `contain` doesn't seem to work...

```markdown
---
marp: true
---

## Right-Sided Chart
- The chart is randomly found.

![bg contain right](https://seaborn.pydata.org/_images/grouped_barplot.png)
```

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Slide1_20220622_1655885449.jpg)

### Footers and Page Numbers

Adding two settings at the beginning will enable footers and page numbers. However, their positions can't be adjusted using markdown syntax.

```markdown
---
marp: true
footer: 20220610

---
## Slide with Page Number
```

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Slide2_20220622_1655885767.jpg)
```