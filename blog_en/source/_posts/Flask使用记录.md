---
title: Flask Usage Record
categories: Script
date: 2020-08-09 01:25:45
tags: ['Flask', 'Python']
---

Due to project requirements, I finally took the time to actually look into how to use Flask. Here's a record of what I learned.
<!-- Abstract part -->
<!-- more -->

Since my task this time was relatively simple - just create a Flask server that can accept requests and call Python code to compute results - I only dealt with the most basic usage. Of course, since I'm not familiar with HTML and JavaScript, my understanding here is purely based on my own experience... it may contain many inaccuracies and misleading descriptions...

## Creating a Flask Service (App)

All operations start with creating a Flask App. It's quite simple - just instantiate the `Flask` object.

```python
from flask import Flask
app = Flask(__name__)

if __name__ == '__main__':
    app.run(host='0.0.0.0')
```

Save this as a Python file, and then run `python app.py` to start the service. Setting the `host` parameter to `'0.0.0.0'` allows all IPs to access this service.

## Binding Requests to Functions

Flask can convert requests for pages into calls to Python functions. Specifically, use the `@app.route` decorator before the corresponding function. For `'/'`, it corresponds to `index.html`, and for something like `'/login'`, it corresponds to `login.html`. Of course, I found that the HTML files don't necessarily need to exist; you can directly return HTML strings in the response.

I prepared a page explaining how to use the API, which I read in and returned as is.

```python
@app.route('/')
def index():
    with open('index.html') as f:
        content = f.read()
    return content
```

## Handling POST Requests

For my program, I needed to receive some parameters and uploaded file contents for further processing. Since the file content can't be sent via `GET`, I only handled `POST`.

Flask's built-in `request` object makes it convenient to get all information from the request.

```python
from io import BytesIO
from flask import request
from iNeo_MS import sheet_parse

@app.route('/parse', methods=['POST'])
def parse():
    args = request.args
    method = args['method']
    byteData = request.files.get('file')
    data = sheet_parse(BytesIO(byteData.read()), method)
    return str(data.to_dict('records'))
```

## Testing the API with Python

Testing a `POST` service can't be done directly through a browser, so you need to install specialized plugins... but my home network is unstable... so I still rely on Python.

The specific method is using the `post` method from the `requests` module. Place the parameters in the `params` parameter and read the file content before sending it via the `files` parameter.

```python
import requests

params = {
    'method': 'pNovo',
}

with open("data/results.res", "rb") as f:
    content = f.read()

files = {
    "file" :('result.res', content),
}

response = requests.post(
    'http://127.0.0.1:5000/parse',
    params=params,
    files=files
