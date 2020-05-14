# tmprss2
Code and documentation for QSAR modeling and virtual screening of small molecule TMPRSS2 inhibitors

# Installation
1. Install Miniconda for Python 3 from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)
2. Open a terminal, and navigate to the directory you want this repo to live in.
3. Clone and enter this repo by running
```shell script
git clone --recurse-submodules https://github.com/srensi/tmprss2.git
# (advanced) If you are authenticating GitHub with ssh, use this line instead
git clone --recurse-submodules git@github.com:srensi/tmprss2.git
cd tmprss2
```
4. Install dependencies by running these from the tmprss2 directory:
```shell script
# Install dependencies:
conda env create -f environment.yml
conda activate tmprss2

# Install chemprop as package
pip install -e chemprop_fitting/chemprop_repo

# Install pubchempy
pip install pubchempy

```

Now, you should have the dependencies needed by our code and by chemprop, and you should
be able to do things like ```from chemprop import ...``` in python scripts.
