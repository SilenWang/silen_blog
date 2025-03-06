---
title: Building a Task Monitoring Web App with Dash and Flask
categories: Bioinfomatic
date: 2021-06-21 01:24:23
tags: ['python', 'flask', 'dash', 'mongoDB', 'pymongo', 'flasgger']
---

I've been extremely busy lately, with too much work to finish... so I haven't had much motivation to write. However, this time I've applied what I've learned in the past six months to build my first web app with database operations, which is worth commemorating. Although I'm short on time, I'll start with the text and add images later.
<!-- Summary -->
<!-- more -->

This app implements functionality similar to NCBI Blast, allowing users to submit analysis tasks through a web interface and monitor their status from a task list. The app provides a very basic analysis interface where users can select necessary input files and configuration files (previously logged into the database) to start the analysis task and track its progress from the monitoring bar. Although this web interface is not intended for end-users, it serves as a demo for front-end IT communication and will be used by me for testing in the future. The real utility lies in the built-in Web API, which responds to frontend requests to perform actual business tasks.

## Module Introduction
Although the app's functionality is basic, it involves several components and took considerable time to implement. This post focuses on problem-solving during the implementation process and will cover basic operations of new modules (like Dash and Python MongoDB operations) in future posts.

### Web API
The API is built using Flask because it is described as lightweight and simple compared to Django. Since we are not building a large, high-concurrency, high-computation API, I believe Flask is simpler for our needs.

### Web Interface
The web interface is built using Dash, which is similar to Shiny in R. However, there seems to be little presence of Dash in the Chinese bioinformatics community, with most materials being English machine translations. While English materials are also limited, they mainly rely on official documentation and Q&A. Since I am now primarily using Python, using a Python-based tool like Dash is simpler than mixing it with R. Additionally, although my app is simple, it is several orders of magnitude more complex than the bioinformatics scripts I used to write. Therefore, Python's error tracking is incredibly helpful, saving me a lot of time during debugging.

Another point worth mentioning is that Dash evolved from the Plotly plotting package and provides web apps built on Flask. This makes Dash + Flask an excellent combination for my current application scenario because I can write both Dash functions and Flask functions in one main program file, reducing redundancy while ensuring consistency between what I test and demonstrate and what colleagues obtain through the API.

### Task Control
Task control is based on Snakemake. While I won't go into detail here, after using it for a long time, I have a better understanding of the statement "Snakemake is Python code" in its documentation. This feature allows me to implement custom functionalities not currently available in Snakemake by writing my own functions within Snakemake. For example, automatic task status updates in this application.

### File and Task Records (Database)
Task and file records are stored using MongoDB for data storage. If we were only building this app, we wouldn't need any database; simple files or structured directories could suffice. Using a database was purely for practice, as I have other tasks that will require it later.

Choosing MongoDB over SQLite was also due to future work considerations, so I took the opportunity to practice with MongoDB now.

### API Documentation
API documentation was not initially part of the app's planning. Recently, I've built many simple web services for internal analysis systems, and maintaining separate documents for each API is inconvenient. Additionally, different APIs might be assigned to different people, making a single comprehensive document difficult to maintain. Therefore, I considered whether there was a way to generate API documentation directly within a single Web-App. After some research, I chose `flasgger`, which introduces a series of APIs into Flask, allowing authors to write API documentation directly in Python code. The generated documentation can be accessed via Swagger UI directly from the app path.

## App Overview
This app's logic is quite simple: it accepts requests from the frontend, forwards them to the backend to start analysis processes, and updates task statuses in the database during process execution. The frontend then retrieves the program run status either periodically or based on requests. A rough information flow diagram is shown below.

![Here I've skipped a picture]()

## Implementation Highlights
Although this app is simple, there are issues that need additional solutions during implementation:
1. Ensuring consistency between Dash-App and Flask-API processing, meaning that after testing with Dash-App, no additional changes are needed for the Flask-API part.
2. Mixing Dash-App/Flask-App/API documentation within a single App, rather than building multiple Web-Apps and redirecting between them.
3. Utilizing Snakemake to control processes and update their status in the database without setting up a separate monitoring process.

### Consistency Between Dash-App and Flask-API Processing
This is relatively easy to implement... by separating actual business processing functions from Dash and Flask functions, meaning that both Dash and Flask response functions format data into a unified format before passing it to an actual business function. The processed data is then reformatted for Dash or API-specified formats:

```python
@app.callback(
                [
                    Output('info-div', 'children')
                ],
                [
                    Input('upload-info', 'contents')
                ],
                [
                    State('upload-info', 'filename')
                ]
            )
def step2_dash(content, filename):
    '''
    Dash response function
    '''
    # Need to parse target path from file
    if content:
        _, content_stream = content.split(',')
        info_data = step2(content_stream, filename)
        json_fmt = dumps(info_data, sort_keys=True, indent=4, separators=(",", ":"), ensure_ascii=False)
        return [html.Pre(f'Parsing successful, extracted information from the information table:\n{json_fmt}')]
    else:
        return [html.Pre('No parsing result')]


@app.server.route('/step2', methods=['POST'])
def step2_flask():
    '''
    Flask response function
    '''
     data = loads(flask.request.form.get('data'))
     fileUrl = f'{data["dowloadFileUrl"]}?type={data["type"]}&fileName={data["fileName"]}'
     filename = re.findall(r'fileName=(.+)', fileUrl)[0]
     r = requests.get(fileUrl)
     info_data = step2(r.content, filename, decode=False)
     res_data = {'errCode':0, 'data': info_data, 'msg': ''}
    json_fmt = dumps(res_data, sort_keys=True, ensure_ascii=False)
    return json_fmt.encode(encoding="UTF-8")


def step2(content_stream, filename, decode=True):
    '''
    Actual processing part of step2, split to receive different requests
    '''
    if decode:
        decoded = base64.b64decode(content_stream)
    else:
        decoded = content_stream
    info_data = data_format(decoded, filename)
    return info_data
```

### Mixing Dash-App/Flask-App/API Documentation
Mixing Flask and API documentation is relatively simple because Flasgger is essentially a Flask + API documentation combination. Mixing Dash with Flask is also not difficult; the Dash object has a `server` attribute that stores the Flask App object, so any part that needs to use this object can directly use `app.server`:

```python
import dash

app = dash.Dash("Demo")

app.server.config['SWAGGER'] = {
    'title': 'Mass Spectrometry Information Entry/Report Generation Service Documentation',
    'uiversion': 2
}
swagger = Swagger(app.server)

app.layout = html.Div(
    [
        html.H1('Demo'),
        html.Div([
            dcc.Upload(
                id='upload-info',
                children=html.Div([
                    'Upload here',
                ]),
            ),
            html.Div([
                html.Pre(''),
            ], id='info-div')
        ])
    ],
)

@app.callback(
                [
                    Output('info-div', 'children')
                ],
                [
                    Input('upload-info', 'contents')
                ],
                [
                    State('upload-info', 'filename')
                ]
            )
def step2_dash(content, filename):
    '''
    Dash response function
    '''
    return [html.Pre('parsed')]


@app.server.route('/step2', methods=['POST'])
def step2_flask():
    '''
    Flask response function
    ---
    tags:
        - Step2
    parameters:
        - name: dowloadFileUrl
          in: query
          type: string
          required: true
          description: File download URL
    responses:
        200:
            description: Returns parsed data
            schema:
                properties:
                    errorCode:
                        type: number
                        description: Status code
                        default: 0
                    data:
                        type: string
                        description: Task submission key
                        default: parsed
                    msg:
                        type: string
                        description: Empty string when normal, error message if status code is not 0
                        default: ""
    '''
    return 'parsed'
```


### Utilizing Snakemake to Control Processes and Update Status

Initially, I considered using Snakemake's built-in `onstart` and `onerror` features to implement this, as I essentially needed task start and end information submitted to the database. However, after actual implementation, I found that these features are effective for the entire Snakefile, executing once at the start of the total task and once when it fails, without providing responses for individual tasks. Therefore, I had to do it myself...

At this point, we can leverage the "Snakemake is Python code" feature. When specifying the `shell` keyword in a rule, it actually calls Snakemake's internal API `shell` function to execute system-level commands, so writing `run: shell('BASH_SCRIPT')` is essentially the same (there are some minor differences in program behavior that we can ignore). However, by using this approach, we can add additional Python code before and after the system command to implement our desired functionality, such as updating task status before and after execution:

```smk
rule step1:
    input:
        done = '{rid}/analysis.done'
    output:
        txt = '{rid}/done.txt',
    run:
        # update_status is the function to update status, rid is the task ID
        update_status(rid, 'step1', 'running')
        shell(
            '''
            cp {input.done} {output.txt}
            '''
        )
        update_status(rid, 'step1', 'done')
```

This conveniently achieves our goal. Additionally, for error handling, we don't need to do this because the error status in this application does not require granularity down to individual rules; `onerror` is sufficient.

However, one thing to note is that Snakemake actually supports a service called papernote, which specifically collects task statuses from the main control flow. This means Snakemake has something similar built-in, but we would need to read the relevant code to understand how it works and whether there are better internal APIs available.

## Follow-up
Last year, I interviewed a bioinformatics engineer who was unable to return to the company during the pandemic. He learned coding at home and built an analysis submission/monitoring service. Now that I've also done it... looking back, it's not difficult, but it requires time and external motivation to push me to learn. Additionally, practice is important; only through practice can we discover actual problems and solve them, thereby improving our skills. I need to keep working hard...
