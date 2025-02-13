---
title: Installing R Packages with Conda
categories: Others
date: 2018-08-28 22:37:57
tags: ['R', 'conda', 'jupyter']
---

Due to work requirements, I needed to install a package mentioned in a literature. Initially, I thought it would be done with two commands, but who knew there were so many dependencies to resolve... and R packages on Linux need to be compiled, which is quite time-consuming... Fortunately, there's Conda!

<!-- more -->

Speaking of which, when reading wechat public accounts before, I didn't know why they recommended using Conda to install necessary tools... This experience has made me realize how annoying the dependencies in Linux can be... And as a normal user without full permissions, many paths are written forbidden... Conda's method is very useful, even though it wastes space.

# Conda Installation

- Conda comes in two versions... I don't know if they can be called distributions. One is Anaconda, and the other is Miniconda. You can think of them as full installation and minimal installation. Anaconda comes with a complete Python environment, even Jupyter is included, while Miniconda only includes the basic Python and not much else. Since I mainly use it to install R, I chose [Miniconda](https://conda.io/miniconda.html).
- The downloaded Linux installation file is an exceptionally large shell script. Running it will start the installation.
- After installation, you can use the `conda` command to install what you want. For example, I can start installing R.

```bash
conda install r
```

- Everything installed using Conda is uniformly placed in the Conda installation directory. The folder structure inside seems a bit different from the root directory.
- When initially completing the Conda installation, if no special instructions are given, Conda will modify the original `~/.bash` file and add the Conda bin directory to the `PATH` variable. However, its priority is after the system default directories (`/bin`). Since I want to use Conda's configuration to override the system defaults, I manually added this line to `~/.bash_profile` and placed Conda's bin directory at the front of the system defaults.
- Conda also has mirror sites in China. Edit `~/.condarc` (create it if it doesn't exist) and add the following content to use the USTC mirror (Tsinghua's reportedly has issues, so I didn't use it). Of course, this is just an example; under `pkg`, there are other sources as well. Add them if needed.

```yaml
channels:
  - https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
  - defaults
show_channel_urls: true
```

- After that, excute R and install the relevant packages inside. The installed R packages will not be placed in `~/R` but will also be entirely in Conda's directory. Refer to the installation output for specifics.
- Since R packages need to be compiled, I used Miniconda, so there were unsatisfied dependencies. Those packages can be installed using `conda install`. Additionally, you can check all available versions on Anaconda's website and then execute the corresponding commands.

This time, with Conda's assistance, most dependency issues were resolved. The only one that couldn't be resolved was compiled locally and uploaded... It took two days to finally package it up... I hope there won't be any problems when using it.

# Jupyter Installation and Configuration

R is on the cluster, which makes debugging somewhat inconvenient since I'm used to writing and testing in RStudio or Jupyter Notebook. Fortunately, Conda can also install Jupyter, so I searched for remote access configuration methods.

- First, use Conda to install Jupyter.

```sh
conda install jupyter
```

- After completion, configure Jupyter.

```sh
# Generate the configuration file
jupyter notebook --generate-config
# Set a password. Follow the prompts and enter twice.
jupyter notebook password
```

- Open the generated configuration file with an editor (for regular users, it's usually `~/.jupyter/jupyter_notebook_config.py`, as indicated by the prompt when generating the file).
- Find the following items, uncomment them, and change the corresponding values. The password is in the JSON file mentioned during generation; there will be a long string inside. Copy and paste it all over.

```txt
c.NotebookApp.ip='*'
c.NotebookApp.password = u'your_key_str'
c.NotebookApp.open_browser = False # Represents starting the notebook service without opening the browser and accessing it.
c.NotebookApp.port = 8888 # This can be omitted; a port will be automatically assigned.
```

- Then, you can open your browser on any computer that can connect to the cluster IP with `0.0.0.0:8888` and access it by entering the password.
