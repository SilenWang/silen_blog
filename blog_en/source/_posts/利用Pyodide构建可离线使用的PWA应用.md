---
title: Building an Offline-Capable PWA with Pyodide
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

I recently helped a friend build a small tool that processes sequences from an input file and generates XML in a specific format. My friend isn't comfortable with the command line, so a GUI was needed. The input data also had confidentiality requirements, meaning all computation had to happen locally. Normally I'd build a desktop app, but considering long-term maintenance, I wanted to try an offline-capable PWA. I checked and found that Pyodide already ships pandas as a wasm package, so there was nothing extra to do on the runtime side — let's go!

<!-- more -->

## The Former Barrier to Building a PWA

The concept of PWA (Progressive Web App) has been around for years. Its core capability is offline access and background sync through Service Workers. But building an offline PWA tool used to be quite challenging for me — I can't write JS or TS from scratch.

## Pyodide: Running Python in the Browser

[Pyodide](https://pyodide.org/) compiles CPython to WebAssembly, letting us run Python code directly in the browser. It doesn't just port the Python interpreter — it also does extensive adaptation through Emscripten so that numpy, pandas, scipy, matplotlib and other scientific computing libraries work in the browser too. These libraries aren't written in Python under the hood, so the interpreter alone wouldn't be enough.

## Combining with PWA for Offline Capability

Pyodide solves the "run Python in the browser" problem. Combined with PWA's offline capabilities, we get a complete toolchain that doesn't depend on the network.

The architecture looks like this:

1. **Core logic**: Python scripts executed via Pyodide in the browser
2. **Frontend**: HTML + CSS + JS for file upload, parameter config, and result display
3. **Service Worker**: Caches all static assets (HTML, JS, CSS, Pyodide's wasm files, Python packages, etc.) for offline use

Pyodide's wasm files and Python packages (.whl or .data) aren't small — the first load requires a network. But once cached by the Service Worker, everything works fine even when fully offline.

## AI Did Most of the Coding

The architecture made the idea feasible. But in the past, actually writing this program would still have been painful — multi-language coding, and using a Python framework for the web UI wasn't worth it (it would cache more Python libraries and might introduce libraries that are hard to compile to wasm).

Now there's DeepSeek. The UI for this tool is just two buttons and a log box — I could leave everything to DeepSeek.

## Actual Results and Experience

The final result: on first visit, users need to be online to load resources (mainly Pyodide's wasm files and Python packages). After that, even with no network at all, refreshing the page works and all features function normally.

As for speed, running Python through WebAssembly is slower than native CPython. But this tool processes fewer than 10,000 records — performance isn't a concern.

I can see myself using this combo to build more small tools in the future.
