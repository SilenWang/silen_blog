---
title: Building an Offline‑Capable PWA with Pyodide
date: 2026-05-28 07:30:00
tags:
  - Pyodide
  - PWA
  - WebAssembly
  - Python
  - Service Worker
  - JavaScript
categories: Coding
---

I recently helped a friend build a small tool that needs to run Python code in the browser for data processing, and the user’s environment may not have a network connection. A few years ago this requirement would have been nearly impossible – you either had to set up a backend server to run Python, or rewrite the entire logic in pure JavaScript. But with the combination of Pyodide + PWA, it becomes feasible.

<!-- more -->

## Background First

My friend’s requirement was this: there is a fairly specific data processing workflow, the code is written in Python, and it depends on some common scientific computing libraries (pandas, numpy, etc.). Users need to upload files in the browser, process them, and download the results, and often they are in an intranet environment without a network, or in airplane mode.

This scenario is actually quite typical. Many small tools in bioinformatics and data analysis have their core logic written in Python, but the distribution method is either to have everyone install a Python environment and run scripts, or to set up a web service. The former is not friendly to non‑technical users, and the latter requires server costs and maintenance.

## The Former Barriers to Building a PWA

The concept of PWA (Progressive Web App) has been around for several years, and its core capability is offline access, background sync, etc. through Service Workers. But in the past, building a practical PWA tool faced two major problems:

### 1. Business Logic Had to Be Rewritten in JS

If your existing core logic is written in Python and you want to run it in the browser, there were basically only two paths:

- **Rewrite as JavaScript**: Sounds simple, but Python’s scientific computing ecosystem (pandas, numpy, scipy, etc.) has no exact equivalent in JS. Even if there are libraries with similar functionality, the differences in behavior details are enough to keep you debugging for a while.
- **Set up a backend service**: The frontend is only responsible for display, the backend runs the Python logic, and they communicate via API. This approach is feasible, but it introduces server costs, network dependency, and operational complexity. Moreover, if users need offline use, this approach won’t work.

### 2. The Hard Requirement of HTTPS

The core feature of PWA (Service Worker) requires an HTTPS environment to register (except localhost). This means that if you want to make a PWA tool for others to use, you first have to solve the HTTPS problem.

For a production environment with a fixed domain name, applying for an SSL certificate is not a problem. But if you just want to make a temporary tool to share with friends or colleagues, it becomes a bit awkward:

- Applying for a domain name and certificate for a temporary tool is too much trouble
- Using a self‑signed certificate will cause the browser to report it as unsafe
- Using tools like ngrok can quickly expose an HTTPS port, but they rely on third‑party services

In practice, during development and testing, the best solution is to directly use `localhost` (the browser exempts localhost from the HTTPS requirement), and then use the `--host` parameter to allow others on the same LAN to access it. But this method only works on a LAN and is limited to development and testing.

## Pyodide: Running Python in the Browser

[Pyodide](https://pyodide.org/) is a project that compiles CPython to WebAssembly, allowing you to run Python code directly in the browser. It doesn’t just move the Python interpreter into the browser; it also does a lot of adaptation work through Emscripten, enabling commonly used scientific computing libraries such as numpy, pandas, scipy, and matplotlib to run in the browser as well.

Usage is simple: include it in your HTML:

```html
<script src="https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js"></script>
<script>
  async function main() {
    let pyodide = await loadPyodide();
    await pyodide.loadPackage(['numpy', 'pandas']);
    // Now you can run Python directly in the browser
    let result = await pyodide.runPython(`
      import numpy as np
      import pandas as pd
      data = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
      data.describe()
    `);
    console.log(result);
  }
  main();
</script>
```

Of course, in actual use you wouldn’t directly write large strings to pass Python code. A more reasonable approach is to treat the Python script as a separate file and load and execute it through the file system API provided by Pyodide.

```javascript
// Write the Python file into Pyodide's virtual file system
pyodide.FS.writeFile('/main.py', pythonCode);

// Execute
let result = pyodide.runPython(`
  import sys
  sys.path.insert(0, '/')
  import main
  main.run()
`);
```

## Combining with PWA for Offline Capability

With Pyodide solving the “run Python in the browser” problem, combining it with PWA’s offline capability makes it possible to have a complete toolchain that does not depend on the network at all.

The overall architecture is roughly like this:

1. **Core logic**: Python script, executed on the browser side via Pyodide
2. **Frontend interface**: HTML + CSS + JS, providing file upload, parameter configuration, and result display
3. **Service Worker**: Caches all static resources (HTML, JS, CSS, Pyodide’s wasm files, Python packages, etc.) to ensure offline availability
4. **IndexedDB**: If you need to persist user data or cache processing results, you can store them in IndexedDB

There is a key point here: Pyodide’s wasm files and Python packages (.whl or .data files) are not small in size and require a network connection for the first load. But once they are cached via the Service Worker, subsequent use works normally even when completely offline.

```javascript
// Caching strategy in the Service Worker
const CACHE_NAME = 'pyodide-app-v1';
const PRECACHE_URLS = [
  '/',
  '/index.html',
  '/app.js',
  '/style.css',
  '/main.py',
  // Pyodide related resources
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.asm.wasm',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.asm.data',
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME).then(cache => cache.addAll(PRECACHE_URLS))
  );
});
```

## Solutions for Temporary HTTPS Access

As mentioned earlier, PWA requires HTTPS, and temporary tools sometimes make it inconvenient to deal with domain certificates. Here are a few solutions I have tried:

### Option 1: localhost + LAN Sharing

If it’s for development testing or small‑scale use, simply start an HTTPS Server locally (using a self‑signed certificate or mkcert to generate a trusted local certificate), then use `--host 0.0.0.0` to listen on all network interfaces. Devices on the LAN can then access it via `https://192.168.x.x:port`.

```bash
# Use mkcert to generate a locally trusted certificate
mkcert -install
mkcert localhost 192.168.1.100 127.0.0.1

# Use Node's https module or Caddy to start the server
caddy file-server --domain localhost
```

The certificate generated by mkcert will be trusted by the system, the browser will not report it as unsafe, and the experience is close to a production environment.

### Option 2: Cloudflare Tunnel

If you need a publicly accessible HTTPS address for a demo (for example, to let friends in different cities test it), Cloudflare Tunnel is a good choice. It establishes a tunnel to Cloudflare’s edge nodes and automatically provides HTTPS.

```bash
cloudflared tunnel --url http://localhost:8080
```

After execution, you will get an address like `https://xxx.trycloudflare.com`. Just open it in the browser, and both HTTPS and PWA will work normally. No domain name is needed, no certificate configuration is required, and you can just shut it down when you are done.

### Option 3: GitHub Pages

If the tool itself is purely frontend (no backend API involved), you can also deploy it directly to GitHub Pages. It automatically provides HTTPS and has a fixed URL, making it suitable for long‑term tools. The downside is that updates require going through the git workflow, which is less flexible than local development.

## Actual Results and Experience

With this combination, the final effect is: the first time a user visits, they need to be online to load resources (mainly Pyodide’s wasm files and Python packages). After that, even if they are completely disconnected from the network, refreshing the page still allows them to use all features normally.

In terms of processing speed, Pyodide runs Python through WebAssembly. Although its performance is not as good as native CPython, it is sufficient for most data processing tasks. Especially in scenarios with small amounts of data (a few MB to tens of MB of tabular data), the processing time is basically on the order of seconds, and users do not perceive any noticeable delay.

## Summary

The Pyodide + PWA combination turns the idea of “running Python code offline in the browser” from something that sounds cool but unrealistic into a truly deployable solution.

Its core advantages are:

- **No backend needed**: Pure frontend architecture eliminates server costs and maintenance
- **Offline capable**: After Service Worker caching, it works even when completely disconnected
- **Reuse the Python ecosystem**: No need to rewrite existing Python code into JS; libraries like numpy and pandas can be used directly
- **Zero deployment cost**: Static files can be run anywhere – GitHub Pages, a local server, or even opened from a USB drive

Of course, it also has limitations: it is not suitable for computationally heavy tasks (Pyodide has its performance ceiling), it is not suitable for scenarios that require calling system APIs (browser sandbox restrictions), and the initial resource load size is not small.

But for the scenario of “having a Python‑written tool that you want non‑technical users to be able to use offline in the browser,” this is probably the best solution available today.
