# tmprss2
Code and documentation for QSAR modeling and virtual screening of small molecule TMPRSS2 inhibitors

# Installation
1. Clone or download this repo.
2. Install Miniconda for Python 3 from [https://conda.io/miniconda.html](https://conda.io/miniconda.html)
3. Run these in terminal
```shell script
#Install dependencies:
cd /path/to/tmprss2
conda env create -f environment.yml
conda env update -n tmprss2 -f chemprop_fitting/chemprop_repo/environment.yml
conda activate tmprss2
```

```shell script
# Install chemprop as package
pip install -e chemprop_fitting/chemprop_repo
```

Now, you should have the dependencies needed by our code and by chemprop, and you should
be able to do things like ```from chemprop import ...``` in python scripts.
