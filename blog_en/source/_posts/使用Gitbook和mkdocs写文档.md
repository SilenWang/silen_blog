---
title: 'Using Gitbook and mkdocs for Documentation'
categories: Others
date: 2019-07-07 22:13:31
tags: ['documentation', 'gitbook', 'mkdocs']
---

My blog is generated using Hexo, a static blog framework based on Node.js. It helps people like me who know nothing about frontend development to quickly set up a decent-looking blog, while I only need to know how to write posts in Markdown. Similarly, during development there's always a need to write technical documentation, hence tools like Gitbook and mkdocs exist.

<!-- more -->

Gitbook and mkdocs actually have different positioning. Gitbook was originally designed for quickly writing books (including documentation) using Markdown, while mkdocs helps technical personnel quickly write and maintain easily modifiable documentation. But in the hands of beginners like me...their basic functionalities are quite similar, so there's not much difference...

## Gitbook

`gitbook-cli` is a Node.js-based tool that can be installed via npm: `npm install gitbook-cli`. Note that due to changes in the Gitbook project's direction, this CLI tool is no longer being updated. However, the basic features are sufficient and it still works fine...

After installation, you can use the `gitbook` command to create and write documentation. The creation command is `gitbook init`, which generates `README.md` and `SUMMARY.md` files in the directory. `SUMMARY.md` serves as both the table of contents and structure file - `gitbook-cli` uses this file to generate the final documentation. `SUMMARY.md` itself is also a Markdown file:

```markdown
# Chapter 1
* [Documentation](README.md)
* [Part 1](Part1/README.md)
    * [Content 1](Part1/content1.md)
        * [Sub-content](Part1/content1.md)
# Chapter 2
* [XXXX](NNNN.md)
```

Book chapters are controlled by first-level heading syntax (`#`). Different chapters are separated by dividers in the table of contents tree. If auto-numbering is enabled, the numbering will differ between chapters.

Chapter contents are stored as unordered lists, with each list item using link syntax (`[]()`) to point to specific Markdown files. Within these Markdown files, first-level headings serve as anchors. Currently only first-level headings work as anchors - if you need to anchor to second-level headings or other content, you'll need to manually add HTML tags as anchors.

Gitbook supports all basic Markdown syntax as well as HTML, so if you need to include complex tables, you can write them in HTML format and paste them directly.

With this knowledge, you can basically start writing books. If you need to customize the appearance or use plugins, you can create a `book.json` configuration file and fill in the settings accordingly.

After writing, you can use `gitbook serve` to preview the book locally (default port 8000). The preview service can remain running, and modifications to Markdown source files will be reflected in real-time (may require refresh). When you're satisfied with the book, use `gitbook build` to generate HTML files in the `_book/` directory.

## mkdocs

`mkdocs` is a Python-based tool that can be easily installed via `pip` or `conda`. After installation, the `mkdocs` command becomes available. To create documentation, use `mkdocs new proj_name`, which creates an `mkdocs.yml` file and `docs/` folder in the `proj_name/` directory. `mkdocs.yml` serves as both configuration and table of contents file, using YAML format:

```yaml
site_name: "My Documentation"

pages:
- Table of Contents: index.md
- Chapter 1: chapter1.md
- Chapter 2: chapter2.md

theme: readthedocs
```

mkdocs generates a table of contents similar to common Markdown TOCs, without needing additional anchor settings like Gitbook. Since my use case was relatively simple, I didn't explore more advanced features.

Completed documentation can be previewed using `mkdocs serve`, which also updates in real-time. Alternatively, use `mkdocs build` to generate static files in the `site/` directory.

## Comparison of the Two Tools

I've used both tools while preparing handover documentation and have some comparative observations.

In terms of basic functionality, there's almost no difference for my needs - both use Markdown for writing documentation, have simple commands, and support local preview. However, since gitbook-cli is no longer updated, there might be future compatibility issues.

Regarding extensibility...they're also quite similar, with many plugins available for both. One notable difference is that Gitbook seems to lack third-party themes? I searched many times but couldn't find ready-to-use themes, while mkdocs comes with several built-in ones.

Also, since this was handover documentation for work, my colleagues weren't accustomed to reading HTML, so I tried PDF export plugins for both tools...they work...but not very elegantly...compared to the HTML output, both PDFs look rather crude...though if I had to choose, mkdocs' PDF output is slightly better.