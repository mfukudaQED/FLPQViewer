# About
FLPQViewer is a code to visualize local physical quantities.

# Installation
Depending on the Python version, it may not function correctly.
Therefore, we recommend setting up a conda environment, as described later, for the installation.

```
git clone https://github.com/mfukudaQED/FLPQViewer.git
cd FLPQViewer
pip install -e .
```


# Usage

```
python ${PATH_FLPQV}/flpqv/main_flpqv.py input.toml
```


The default.toml file is available at `${PATH_FLPQV}/flpqv/default.toml`.

# Installation using conda environment
## Setting of miniconda
### Installation of miniconda
- Download the installer from the official website and install it.
```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ bash ./Miniconda3-latest-Linux-x86_64.sh
```

- During the installation, the following initial conda configuration is automatically added to the `.bashrc` file.
```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/fukuda/program/miniconda/exit/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/fukuda/program/miniconda/exit/etc/profile.d/conda.sh" ]; then
        . "/home/fukuda/program/miniconda/exit/etc/profile.d/conda.sh"
    else
        export PATH="/home/fukuda/program/miniconda/exit/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

### Set up for conda environment
- Creating a clean virtual environment for FLPQViewer,
```
$ conda create -n flpqviewer_env python=3.11
```

- When using FLPQViewer, always activate this environment with conda activate.
```
$ conda activate flpqviewer_env
```

- If you want to leave this flpqviewer_env environment,
```
$ conda deactivate
```

## Installation of FLPQViewer
- Install FLPQViewer using the pip command while the conda environment is activated.
```
$ git clone https://github.com/mfukudaQED/FLPQViewer.git
$ cd FLPQViewer
$ pip install .
```

### When packages such as Python conflict with each other
- Spack or similar tools may install their own version of Python, and in such cases,
  activating a Conda environment with conda activate can still result in the Python path being overridden by the Spack-managed one.
  In this situation, it is recommended to use conda run to ensure that the Python executable from the intended Conda environment is used.
```
cd FLPQViewer
conda run -n flpqv_env pip install .
conda deaactivate
conda run -n flpqv_env python3 FLPQViewer/flpqv/main_flpqv.py chempot.toml
```

# Contact
Masahiro FUKUDA (ISSP, Univ. of Tokyo)  
masahiro.fukuda__at__issp.u-tokyo.ac.jp  
Please replace `__at__` by @.
