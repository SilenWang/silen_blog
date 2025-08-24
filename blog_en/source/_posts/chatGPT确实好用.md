---
title: chatGPT is indeed useful
tags:
  - chatGPT
  - python
  - flask
categories: Others
date: 2023-02-23 18:24:06
toc: False
---

chatGPT has been popular for so long, but I only created an account yesterday when I wanted to use it to generate some content. As many people online have said, while it may be unreliable on professional questions, for writing code - since I ultimately need to run and test it myself - practice yields true knowledge. Having a roughly correct direction is already helpful!

<!-- Excerpt -->
<!-- more -->

The background of this problem is:

- I used `mkdocs` to generate a technical documentation for company use (similar to my manual), which is a static webpage that can be viewed directly by adding an HTTP server
- I want to embed this documentation into an existing system using `iframe`
- Since this documentation is important, it has certain confidentiality requirements, so it cannot be accessed directly and needs some verification

After discussing with IT colleagues, the final solution was that they would give me a token on the first request, I would use this token to verify against their API, and return the content only if the verification passes. I only know basic flask and fastAPI, and flask is simpler to use, so in summary, I need to write a program that:

- Uses flask to build an HTTP server with the same functionality as Python's simple http.server, allowing users to access the static website created by mkdocs
- This server needs verification functionality - receiving a token from the request to verify whether the user is legitimate
- After the first verification, requests for internal jumps within the static website are handled by this service itself, so to maintain continuous verification, I need to save the token myself and complete the verification

Although I've used flask and dash before, I actually don't have much concept of web-related technologies. Previous experience testing caddy reverse proxy made me think this program was certainly implementable, but I wasn't sure exactly how. The colleague's mention of saving tokens myself made sense conceptually but I didn't know how to implement it.

This is where chatGPT came into play. I described my requirements roughly and asked how to write the code I wanted. It didn't give me a perfect solution matching all my requirements, but provided a code example that let me test and figure out how the program should be written. For parts I didn't understand, I could ask follow-up questions or request specific explanations. While its answers weren't completely correct (and sometimes due to my unclear questions), through discussion I gradually found the implementation I wanted.

In the end, I completed this simple program in one morning. While the program is simple, in the past I would have needed to:
1. Find several Chinese web pages to see if there were ready examples
2. Guess the English description methods
3. Browse about 10 English posts on stackoverflow
4. Look for answers to specific parts of my problem
5. Try to assemble a workable solution...

This process wouldn't be painful - it's a learning process after all - but certainly not efficient...

Now, as long as I can basically describe what I want, chatGPT can directly give me a code example to debug and learn from, saving most of the time searching Chinese and English materials and jumping straight to debugging with examples. The efficiency improvement is transformative...

Suddenly I feel that with its help, learning some frontend tricks... writing my own websites and apps is no longer a dream...

#### update20220227

Today chatGPT has too many requests and is temporarily unavailable (not sure about paid users). Seems like early adoption brings benefits... play with it while you can before free access disappears...


{% noteblock:: Appendix: Complete conversation process about this issue %}

#### You:
How to use flask to build an HTTP server that uses tokens for verification?

#### ChatGPT:
To use Flask to build an HTTP server that uses tokens for verification, you can use Flask's `Flask-JWT` extension. This extension provides JSON Web Token (JWT) authentication mechanism for Flask applications.

First install required extensions:
```
pip install Flask
pip install Flask-JWT
```

Then write Flask application code as follows:
```python
from flask import Flask, jsonify, request
from flask_jwt import JWT, jwt_required, current_identity
from werkzeug.security import safe_str_cmp

class User(object):
    def __init__(self, id, username, password):
        self.id = id
        self.username = username
        self.password = password

users = [
    User(1, 'user1', 'password1'),
    User(2, 'user2', 'password2')
]

def authenticate(username, password):
    for user in users:
        if user.username == username and safe_str_cmp(user.password.encode('utf-8'), password.encode('utf-8')):
            return user

def identity(payload):
    user_id = payload['identity']
    return next((user for user in users if user.id == user_id), None)

app = Flask(__name__)
app.config['SECRET_KEY'] = 'super-secret'
jwt = JWT(app, authenticate, identity)

@app.route('/protected')
@jwt_required()
def protected():
    return jsonify({'user_id': current_identity.id, 'username': current_identity.username})

if __name__ == '__main__':
    app.run()
```

#### You:
This flask server will receive a token when processing requests. I need to use this token to verify against another server, and only return static content after verification passes. How to implement this?

#### ChatGPT:
To use the received token to verify against another server in Flask, you can use Python's `requests` library to conveniently send HTTP requests to other servers and get responses.

Example:
```python
import requests
from flask import Flask, request

app = Flask(__name__)

def verify_token(token):
    response = requests.get('https://example.com/verify-token', headers={'Authorization': 'Bearer ' + token})
    return response.status_code == 200

@app.route('/')
def index():
    token = request.headers.get('Authorization').split(' ')[1]
    return 'Hello, world!' if verify_token(token) else 'Invalid token'

if __name__ == '__main__':
    app.run()
```

#### You:
If I want to return not a simple website but pre-generated static pages like those from mkdocs or gitbook, how to implement?

#### ChatGPT:
To return pre-generated static pages, use Flask's `send_from_directory` method to serve files from a specific directory.

Example assuming static files are in `static/docs`:
```python
from flask import Flask, send_from_directory

app = Flask(__name__)

@app.route('/docs/<path:filename>')
def serve_docs(filename):
    return send_from_directory('static/docs', filename)

if __name__ == '__main__':
    app.run()
```

#### You:
How does flask get parameters from GET requests?

#### ChatGPT:
In Flask, use `request.args` to access GET parameters. It's a dict containing all parameters and their values.

Example for URL `http://localhost:5000/search?q=python&page=1`:
```python
from flask import Flask, request

app = Flask(__name__)

@app.route('/search')
def search():
    query = request.args.get('q')
    page = request.args.get('page')
    return f'Searching for "{query}" on page {page}'

if __name__ == '__main__':
    app.run()
```

#### You:
The above code can only access first-level files under static/docs/. How to access files in subdirectories?

#### ChatGPT:
To access files in subdirectories, use `os.path.join` to combine paths:

```python
import os
from flask import Flask, send_from_directory

app = Flask(__name__)

@app.route('/<path:filename>')
def serve_docs(filename):
    if filename == '':
        filename = 'index.html'
    path = os.path.join('static', filename)
    return send_from_directory(path, filename)

if __name__ == '__main__':
    app.run()
```

#### You:
In flask, how to redirect '/' to '/index.html'?

#### ChatGPT:
Use Flask's `redirect` function:

```python
from flask import Flask, redirect

app = Flask(__name__)

@app.route('/')
def redirect_to_index():
    return redirect('/index.html')

if __name__ == '__main__':
    app.run()
```

#### You:
How to generate a recommended secret key for flask?

#### ChatGPT:
Use Python's `secrets` module to generate secure random keys:

```python
import secrets

app = Flask(__name__)
app.secret_key = secrets.token_hex(16)  # 16-byte random string
```

#### You:
How to implement token expiration in flask?

#### ChatGPT:
For token expiration, consider using JWT (JSON Web Tokens) with the `cryptography` library:

```python
import datetime
import jwt

app = Flask(__name__)
app.config['JWT_SECRET_KEY'] = 'mysecretkey'
app.config['JWT_EXPIRATION_DELTA'] = datetime.timedelta(seconds=3600)

@app.route('/login')
def login():
    payload = {
        'sub': '1234567890', 
        'exp': datetime.datetime.utcnow() + app.config['JWT_EXPIRATION_DELTA']
    }
    token = jwt.encode(payload, app.config['JWT_SECRET_KEY'], algorithm='HS256')
    return {'access_token': token.decode('UTF-8'), 'token_type': 'Bearer'}

@app.route('/protected')
@jwt_required()
def protected():
    return {'data': 'Protected resource'}
```

{% endnoteblock %}
