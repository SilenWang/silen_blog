---
title: Basic Usage of GitHub Actions
categories: Script
date: 2025-04-13 16:40:47
tags: ['github', 'github actions']
---

In modern software development, Continuous Integration (CI) and Continuous Deployment (CD) are key practices for improving development efficiency and quality. [GitHub Actions](https://github.com/features/actions) is a powerful automation tool provided by GitHub that helps developers easily implement CI/CD and broader automation workflows.

<!-- Excerpt section -->
<!-- more -->

## The Role of GitHub Actions
GitHub Actions is a powerful automation platform that allows developers to define and run workflows in GitHub repositories. It can automatically execute a series of tasks (such as code building, testing, deployment, etc.) based on specific events (like code pushes, pull requests, tag releases, etc.). Through GitHub Actions, developers can achieve full automation from code submission to deployment, saving time and effort while reducing human errors.

## GitHub Actions Workflow Structure
GitHub Actions workflows consist of several components, with the most important being jobs and steps.

### Jobs
A workflow can contain multiple jobs, each being an independent execution unit that can run in different environments.

```yaml
jobs:
  build: # Job name
    runs-on: [ubuntu-latest] # Runtime environment, here specifying GitHub-hosted runners
```

`runs-on` specifies the environment for running the job - here we're using GitHub's latest Ubuntu environment.

### Steps
Each job can contain multiple steps, where each step is an independent task that can execute a series of commands or call predefined actions. Predefined actions can be found in the [GitHub Marketplace](https://github.com/marketplace).

```yaml
- name: Get current tag # Step name
  id: get_tag # Step ID
  uses: zingimmick/github-action-get-current-tag@v1 # Using a specific pre-written step
```

## Variable Passing

Just like in bioinformatics analysis pipelines I've written before, different steps often need to pass information to complete a task. As a simple configuration file, YAML can't implement variable passing directly, so we need to rely on some environment variables set in the GitHub Actions framework.

### Passing Information Between Steps

When executing code in steps, you can pass information between steps by appending content to the special environment variable `$GITHUB_OUTPUT`:

```yaml
- name: Save info
    id: save_info
    run: |
        echo "tag=3" >> $GITHUB_OUTPUT
- name: Use info
    id: use_info
    run: |
        echo ${{ steps.save_info.outputs.tag }}
```

### Passing Information Between Jobs

In jobs, you can specify to pass outputs from specific steps to job inputs as follows:

```yaml
jobs:
    get_info:
        steps:
            - name: Save info
                id: save_info
                run: |
                    echo "tag=3" >> $GITHUB_OUTPUT
    outputs:
      tag: ${{ steps.save_info.outputs.tag }}

    use_info:
        needs: [get_info] # Set dependencies to get information from other jobs
        steps:
            - name: Save info
                id: save_info
                run: |
                    echo ${{ needs.get_info.outputs.tag }}
```

## Example: Building a Project and Publishing

```yaml
name: BUILD

on:
  push:
    tags:
      - "v*" # Triggered by version tags

jobs:
  build:
    runs-on: [ubuntu-latest]
    permissions: # Need write permission for releases
      contents: write

    steps: 

      - name: Get current tag
        id: get_tag
        uses: zingimmick/github-action-get-current-tag@v1
  
      - name: Checkout code
        uses: actions/checkout@v4
        with:
          fetch-depth: 0
          ref: ${{ steps.get_tag.outputs.tag }}

      - name: Build frontend
        run: pixi run front_prod
        
      - name: Build backend
        run: pixi run backend

      - name: Pack to zip
        run: pixi run release

      - name: Release
        uses: ncipollo/release-action@v1.15.0
        with:
          draft: false
          generateReleaseNotes: true  # Automatically generate release notes
          artifacts: '${{ github.workspace }}/release/*.zip'
          tag: ${{ steps.get_tag.outputs.tag }}
          name: Release ${{ steps.get_tag.outputs.tag }}
```
