---
title: Testing Network APIs with pytest
categories: Others
date: 2025-06-09 00:25:03
tags: ['python', 'pytest']
---

For the past six months, a significant part of my work has been maintaining the company's website. Since I took over the code midway, there's always a risk of introducing new bugs during maintenance. Therefore, as required by my supervisor, I learned how to implement automated testing.

<!-- more -->

The testing framework used is pytest, the most popular testing framework in Python.

The framework is very simple to use. Create files starting with `test_`, then write functions starting with `test_` inside them. Within these functions, write test code and add `assert` statements to complete the most basic tests:

```python
BASE_URL = 'http://127.0.0.1:4567'

def test_login_success():
    login_data = {"username": "USER", "password": "123456"}
    response = requests.post(f"{BASE_URL}/login", json=login_data)
    
    # Assert key results
    assert response.status_code == 200
    assert "accessToken" in response.json().get("data")
```

In actual testing, we often need to set up some shared information between different tests. For example, all the APIs I need to test require authentication, so it's perfect to use pytest's `fixture` system. We can first obtain a login token before testing begins, then use the same token to access APIs in each test function:

```python
@pytest.fixture(scope='session')
def prod_admin_auth_header():
    return auth_header("admin", "123456")

def test_get_dictionary_id(prod_admin_auth_header):
    response = requests.get(f"{BASE_URL}/GetDictionaryid", headers=lab_mem_auth_header)
    assert response.status_code == 200
```

Here, functions marked as `fixture` can be directly passed into test functions as parameters by their function names, without requiring additional operations.

Additionally, you can create a special file `conftest.py` in the test folder to store all the `fixture`s that need to be shared across the entire test suite. This way, in other test files, you can directly use the configured `fixture`s without manual imports.
