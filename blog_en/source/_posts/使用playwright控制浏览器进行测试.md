---
title: Using Playwright for Browser Testing
categories: Script
date: 2025-06-09 00:44:33
tags: ['python', 'pytest', 'playwright']
---

For websites requiring interactive operations, using pytest for API testing alone is insufficient as it doesn't cover frontend interactions. Previously I used Selenium for browser control in web scraping tasks, but now I've found that Microsoft's Playwright can also control browsers for automated testing, and it's even simpler to use than Selenium since it supports automatic test case recording.

<!-- more -->

## Installation

Using pixi makes Playwright installation quick, but after installation you still need to call Playwright to install browser dependencies. This differs from Selenium which directly uses browsers supported by the system.

```bash
# Install firefox and chromium
playwright install --with-deps firefox chromium
```

## Using with Pytest

Playwright supports multiple languages and also works with Python's Pytest framework. Using the command `playwright codegen -b firefox http://127.0.0.1:8000` to open the browser, then selecting Pytest in the recording window allows you to record Pytest code. A recorded test would look like this:

```python
BASE_URL = 'http://127.0.0.1:8000'

def test_add_sample(lab_mem_page):
    try:
        lab_mem_page.goto(f"{BASE_URL}/page", wait_until="load")
    except TimeoutError: # Retry once if failed
        lab_mem_page.reload(wait_until="load")

    lab_mem_page.get_by_role("button", name="添加/修改样本").click()
    lab_mem_page.get_by_role("menuitem", name="添加样本").click()
    lab_mem_page.get_by_role("textbox", name="*内部编号").fill("CODE")
    lab_mem_page.get_by_role("textbox", name="*住院号").fill("6666")
    lab_mem_page.get_by_role("textbox", name="*姓名").fill("2222")
    lab_mem_page.get_by_role("textbox", name="*性别").click()
    lab_mem_page.get_by_role("listitem").filter(has_text="男").click()
    lab_mem_page.get_by_role("button", name="提交").click()
    expect(lab_mem_page.get_by_text("添加成功")).to_be_visible()
```

After creating test cases, you can run tests and see if any content needs trimming or adjustment.

## Leveraging Pytest Features

When using Pytest for testing, Pytest's fixture system and special `conftest.py` files are all available, allowing information sharing between tests just like regular Pytest testing to make tests more efficient and stable.

## Issues with Recording System

While the recording feature significantly speeds up test writing, automatic recording doesn't always perfectly select page elements, especially when interactive elements are dynamically generated. Therefore, to write stable test cases, some frontend knowledge is still needed to uniquely identify target elements in the resource tree.
