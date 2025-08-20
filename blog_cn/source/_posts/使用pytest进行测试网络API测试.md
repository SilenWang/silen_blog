---
title: 使用pytest进行测试网络API测试
categories: Coding
date: 2025-06-09 00:25:03
tags: ['测试', 'API测试', 'pytest', 'python', '自动化测试']
---

这半年有相当一部分工作是维护公司的网站，由于是半路接手的代码，难保维护的时候不带来新的bug，因此根据领导的要求，学习了如何进行自动测试。

<!-- more -->

测试使用的是python下使用最多的测试框架，pytest。

框架使用非常简单，创建以`test_`开头的文件，然后在其中编写以`test_`开头的函数之后，在函数内写测试代码，并加上`assert`断言就能完成最基本的测试了：

```python
BASE_URL = 'http://127.0.0.1:4567'

def test_login_success():
    login_data = {"username": "USER", "password": "123456"}
    response = requests.post(f"{BASE_URL}/login", json=login_data)
    
    # 断言关键结果
    assert response.status_code == 200
    assert "accessToken" in response.json().get("data")
```

在实际测试的过程中，我们可能会需要设置一些在不同测试间共用的信息，比如我要测试的API全都是要登录验证才能使用的，因此正好使用pytest的`fixture`系统，在测试开始前先获取登录token，然后每个测试函数都用同样的token来访问API：

```python
@pytest.fixture(scope='session')
def prod_admin_auth_header():
    return auth_header("admin", "123456")

def test_get_dictionary_id(prod_admin_auth_header):
    response = requests.get(f"{BASE_URL}/GetDictionaryid", headers=lab_mem_auth_header)
    assert response.status_code == 200
```

这里，被设定为`fixture`的函数，其函数名可以作为测试函数的参数直接传到测试函数中去，不需要做额外的操作。

另外还可以在测试文件夹创建一个特殊文件`conftest.py`，将整个测试需要共享的`fixture`都放在里面，这样在其他测试文件中，不需要手动导入，也可以直接使用设定过的`fixture`。
