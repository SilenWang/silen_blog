---
title: Hands on Writing an Automated Deployment GitHub Action
categories: Coding
date: 2025-12-17 01:44:13
tags: ['github action', 'CD', 'automated deployment']
---

As the number of projects I need to maintain keeps growing (3 official websites, 2 systems, and 1 mini‑program, all with separate front‑end/back‑end and independent databases), I’ve encountered many tasks that are trivial to do once but become chaotic when repeated many times across different contexts. I had already experimented with CI on GitHub, and this time I learned and practiced CD (Continuous Deployment).

<!-- more -->

## The Pain Points of Maintaining Multiple Projects Alone

As a solo developer who has to handle everything, I need to maintain several projects with different technology stacks. Each project has its own repository, build process, and deployment environment. My previous workflow was: modify code locally → run build scripts manually → upload artifacts to the server via SCP → log into the server to replace files and restart services. This process is highly repetitive and prone to mistakes due to oversight or human error, which could lead to production issues.

Especially after fixing an urgent bug and needing to switch between multiple projects, the time and mental overhead of manual deployment increase significantly. Therefore, I wanted to hand over these repetitive operations to a machine, allowing me to focus more on the actual code logic.

## Why CD (Continuous Deployment) is Needed

CI (Continuous Integration) handles automated building and testing of code, ensuring every commit is runnable. CD (Continuous Deployment) goes a step further by automatically deploying the tested code to a production or staging environment. For teams that require frequent updates and work on multiple projects in parallel, CD can greatly reduce repetitive manual work, lower the chance of human error, and let developers concentrate on feature development.

Of course, I’m not such a team—I’m just one person. But as mentioned above, a usable CD workflow can still free me from mechanical repetition and let me devote my energy to fixing the actual code problems.

## Writing a GitHub Actions CD Workflow

Below is the deployment workflow I designed for the official website projects I maintain. It triggers automatically when code is pushed to the `main` branch and changes occur under `app/Official_Site/**` or `app/Official_Site_EN/**`. It also supports manual triggering via the GitHub interface (`workflow_dispatch`).

```yaml
name: DEPLOY

on:
  push:
    branches:
      - main
    paths:
      # Only trigger when these directories change
      - 'app/Official_Site/**' 
      - 'app/Official_Site_EN/**'
  workflow_dispatch:

jobs:
  build:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        environment: [website]

    steps:
      - name: Checkout repository
        uses: actions/checkout@v3
        with:
          ref: main

      - name: Setup Pixi
        uses: prefix-dev/setup-pixi@v0.9.3
        with:
          pixi-version: v0.59.0
          environments: ${{ matrix.environment }}

      - name: Build Website for ZH/EN
        run: |
          pixi run web_build
          mv app/Official_Site/dist front
          mv app/Official_Site_EN/dist front-en

      - name: Copy files ZH
        uses: appleboy/scp-action@v1
        with:
          host: ${{ vars.PROD_IP }}
          username: root
          password: ${{ secrets.PROD_PASSWD }}
          port: 22
          source: "front"
          target: /opt/project/uploadPath

      - name: Copy files EN
        uses: appleboy/scp-action@v1
        with:
          host: ${{ vars.PROD_IP }}
          username: root
          password: ${{ secrets.PROD_PASSWD }}
          port: 22
          source: "front-en"
          target: /opt/project/uploadPath

      - name: Deploy
        uses: appleboy/ssh-action@v1
        with:
          host: ${{ vars.PROD_IP }}
          username: root
          password: ${{ secrets.PROD_PASSWD }}
          port: 22
          script: |
            cd /opt/project
            rm -rf front.bak front-en.bak
            mv front front.bak
            mv front-en front-en.bak
            mv uploadPath/front front
            mv uploadPath/front-en front-en
            systemctl restart nginx
            
```

### Key Steps Explained

1. **Trigger Conditions (`on`)**
   - `push` to the `main` branch only when the specified directories change, avoiding builds triggered by unrelated commits.
   - `workflow_dispatch` provides a manual trigger for immediate deployment when needed.

2. **Code Checkout**
   - Uses the official `actions/checkout` action to pull the code, specifying `ref: main` to ensure the build is based on the latest commit.

3. **Dependency Installation & Build Artifacts**
   - Employs [Pixi](https://pixi.sh/) as a cross‑language package manager. Here it sets up the environment and builds the static websites with its built‑in command.

4. **File Transfer (SCP)**
   - Uses `appleboy/scp-action` to upload the built artifacts to a temporary directory `/opt/project/uploadPath` on the server.
   - The server IP is passed via `vars.PROD_IP` (a project variable), and the password is stored in `secrets.PROD_PASSWD` (a repository secret) to avoid exposing sensitive information in the code.

5. **Deployment & Rollback (SSH)**
   - Logs into the server via `appleboy/ssh-action` to perform the replacement.
   - Key logic: first back up the currently running `front` and `front‑en` directories with a `.bak` suffix, then move the new versions from the upload directory into place, and finally restart Nginx.
   - The backup operation provides a quick rollback capability—simply move the backup directories back to restore the previous version.

## Conclusion

With this simple GitHub Actions workflow, I’ve successfully reduced the deployment time for static official websites from about 10+ minutes of manual work per deployment to 2‑3 minutes of automated execution, and completely eliminated problems caused by operational oversights.

CD isn’t just for large teams. Individual projects can also achieve significant efficiency gains through automation. If you’re also maintaining multiple projects that require frequent updates, consider spending an hour setting up your own continuous deployment pipeline. I believe it will bring a qualitative leap to your development experience.
