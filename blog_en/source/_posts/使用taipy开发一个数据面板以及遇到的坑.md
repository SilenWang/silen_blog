---
title: Developing a Data Dashboard with Taipy and the Pitfalls Encountered
categories: Coding
date: 2025-10-01 22:49:37
tags: ['Taipy', 'Python']
---

There are really many modules for rapidly developing data or AI-related applications in Python. I've already used Dash, Streamlit, Gradio, NiceGUI, and recently I discovered two more. Just when I needed to develop a simple data dashboard to display company data, I once again recklessly decided to use a new framework - Taipy.

<!-- more -->

## Using Taipy

### Page Layout

Taipy designs multiple ways to write layouts. Depending on the language you're familiar with, you can choose from Markdown, HTML, or Python syntax for page layout settings. However, it seems the official recommendation is Markdown. Yes, it's quite surprising - a library written in Python primarily promotes using Markdown for layouts...

I naturally chose Python syntax. Its usage is basically consistent with Gradio, using `with` statements and objects to control element containment relationships.

Specifically, through the `taipy.gui.builder` module, we can use Python syntax to define page structure.

Main layout components include:
- `tgb.layout()`: Creates a layout container, where you can set the number of columns and gap
- `tgb.part()`: Creates content blocks for organizing related components

In my data dashboard, I adopted a three-level layout structure:
1. Top: Company logo, title, and theme toggle control
2. Middle: Data visualization components like maps, pie charts, line charts, bar charts
3. Bottom: Key metric card displays

Layout example:
```python
with tgb.layout(columns="1 2 1", gap="10px"): # Three columns with width ratio 1:2:1
    # Left block
    with tgb.part(class_name="card"):
        tgb.text("Milestones", class_name="h5")
        tgb.chart(figure="{timeline}", plot_config={"staticPlot": True})
    # Middle and right blocks are similar
```

### Style Configuration

Taipy supports custom styling through CSS files and Stylekit.

- **CSS Styles**: Using CSS is the same as writing regular web pages. I created a `style.css` file to define global styles and component styles:
```css
.card {
    --element-padding: 0.4rem;
}

.scrollable-part {
    overflow-y: auto; /* Vertical scrolling */
}
```

Then when creating the project, specify the additional CSS to load:

```python
gui = Gui(page, css_file="style.css")
gui.run()
```

- **Stylekit Themes**: Define theme colors through Python dictionaries, then pass them when starting the GUI:
```python
STYLEKIT = {
    "color-primary": "#E72410",
    "color-secondary": "#F28E2A",
    "color-background-light": "#EDEDED",
    "color-paper-light": "#FFFFFF",
}

gui.run(stylekit=STYLEKIT, dark_mode=False)
```
Note that Stylekit controls far fewer elements than the classic approach of setting classes and adjusting them in CSS files. For what can be adjusted, refer to [Taipy's documentation](https://docs.taipy.io/en/latest/userman/gui/styling/).

Besides these two methods, different Taipy controls can also use the `class_name` parameter to set some module-preset styles. What's available varies by control, so you need to refer to each control's documentation.

### Value Binding

Taipy provides a method called [Value Binding](https://docs.taipy.io/en/latest/userman/gui/binding/) to make updating interface content simpler.

Specifically, it uses curly brace `{}` syntax to bind Python variables to interface components. Taking the simplest text example, I can set text display as follows:

```python
text = 'My Text'
tgb.text("{text}")
```

In the above step, first assigning the variable `text`, then using `"{text}"` in the text control's content section to bind the variable to the content might seem redundant, but this is very important in Taipy. In Taipy, all objects/variables defined in the namespace are placed in a special object called `state`. In the code above, `state.text` is actually created, and any changes to `state.text` will be reflected in the `tgb.text("{text}")` text control. So if I want to update the displayed text, I just need to modify its content where I can access state, and the text will be modified and displayed on the page, like this:

```python
state.assign(text, "Changed Text")
```

Compared to commonly used Dash/Streamlit/Gradio, this saves writing code to update control content, making development faster. Especially for data display apps, often the chart and table drawing code is fixed, and we're just reorganizing data frames or calculating metrics. Taipy's design allows developers to focus on updating data rather than spending too much time on display aspects.

### Callback Functions

Like other frameworks, various controls in Taipy can also set callback functions. For example, the most basic button has an `on_action` callback, which is triggered when the button is pressed. Different controls have different available callbacks, so you need to refer to the official documentation.

Unlike other frameworks, Taipy can also set global default callback functions. We can directly define a function named `on_action` and place it in the namespace. Taipy's logic is that if a control can trigger an `on_action` callback and there's a function named `on_action` in the current namespace, the callback will be directly triggered. Therefore, in Taipy, using this global callback requires good logic judgment to avoid accidental callback triggers.

### Background Tasks

Although Taipy has built-in scheduling methods, they're placed in the paid version. However, fortunately, in the community version, we can use `invoke_long_callback` for high-load task callbacks to achieve timed data refresh. I referenced [an official example](https://github.com/Avaiga/demo-realtime-pollution):

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

In the code above, `invoke_long_callback` is a function provided by Taipy. We use it to set up a high-load task, which is `user_function`. This function will be executed in a separate thread, then a thread is set up to check the high-load task results, calling the function passed in `user_status_function` to update the data in state. This way, once the data is updated, it triggers Taipy's value binding logic, updating the content displayed on the page.

## Pitfall Records

Originally, I thought that since Taipy has a clear commercialization plan, the basic functions of this module should be quite good. However, during actual use, I encountered pitfalls that stumped me for over a day...

### Toggle Control in Version 4.1.0 Doesn't Trigger in Theme Mode

In version 4.1.0 available through pixi/conda, I found that the Toggle control with `theme=True` (set as theme toggle button) had ineffective value binding and function callbacks... This problem troubled me for a day and a half... I tried having AI give me many solutions, but Taipy doesn't seem to be that popular a module, so I couldn't directly inject JavaScript to solve it like with other frameworks. So I eventually went to look at the source code... I found that in version 4.1.0, the theme toggle control defined in `ThemeToggle.tsx` simply didn't have value binding and callback logic... Only the latest `main` branch had it... Finally, I solved this problem by installing Taipy from source code...

### Stylekit Can't Set All Colors at Once

Although the official documentation says that using Stylekit design allows developers to flexibly and concisely control theme colors, in practice I found that it can't control the page's overall background color, which needs to be solved with additional CSS:

```css
/* Setting theme color */
body {
    background-color: var(--color-background);
}
```

### Plotly Doesn't Have a Function for Drawing Timeline Charts

This... is purely a complaint about Plotly. Although Python is widely used in data science, whether it's Plotly, matplotlib, seaborn, or altair, the richness of their graphics still can't match the ecosystem under R... Of course, someone has written it in R, but package management is too chaotic, and I can't use others' code, which is another matter...

But fortunately, milestone timeline things are just points, lines, and some text, so I finally managed to create it manually with AI's help...

### Insufficient Documentation on Plotly Configuration Methods

This problem also stumped me for a long time. Taipy's documentation explains how to not use Taipy's interface but directly pass Plotly objects to draw graphics on the page. But how to configure Plotly objects, the documentation just says to refer to Plotly documentation. Then I experienced the problem with Plotly's documentation - it seems very detailed but is actually quite poor in content. In their official documentation, many configuration items can only be passed when using `fig.show()`... rather than directly configuring the Plotly object... Then after going in circles, I found that Taipy's `chart` control has a `plot_config` parameter... This can pass parameters equivalent to those in `fig.show()`...
