---
title: 使用playwright控制浏览器进行测试
categories: Script
date: 2025-06-09 00:44:33
tags: ['python', 'pytest', 'playwright']
---


对于要交互操作的网站，光使用pytest进行API测试是不够的，因为这不涉及前端。以前我使用过Selenium控制浏览器来完成爬虫任务，现在发现测试领域，微软的playwright也能控制浏览器以进行自动化的测试，且使用起来比Selenium还简单，因为它支持自动录制测试用例。

<!-- more -->

## 安装

使用pixi可以快速完成playwright的安装，不过安装之后还需要调用playwright来安装浏览器依赖，这一点与Selenium不同，Selenium是直接调用系统中支持的浏览器的。


```bash
# 安装 firfox和 chromium
playwright install --with-deps firefox chromium
```

## 结合Pytest使用

Playwright支持多种语言，同时还支持Python下的Pytest框架，在命令行中使用`playwright codegen -b firefox http://127.0.0.1:8000` 打开浏览器后，在录制窗口中选择Pytest，就可以录制Pytest代码了。录制后一个测试会长这样：

```python
BASE_URL = 'http://127.0.0.1:8000'

def test_add_sample(lab_mem_page):
    try:
        lab_mem_page.goto(f"{BASE_URL}/page", wait_until="load")
    except TimeoutError: # 如果失败，进行一次重试
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

完成测试用例后，可以进行测试，然后看是否有内容需要删减和调整。

## 利用Pytest的特性

在结合Pytest使用进行测试时，Pytest的fixture系统以及特殊的`conftest.py`文件都是可用的，可以像普通Pytest测试一样进行某些信息的共享，以让测试更高效和稳定。

## 录制系统的问题

录制功能确实很大程度上加速了测试的编写，但是自动录制对页面元素的选取不是那么一劳永逸的，尤其当页面上要交互的元素是动态生成的时候。所以要写出稳定的测试用例，还是要懂一些前端知识，能自己找出目标元素在资源树中的唯一定位方式才行。