# tmprss2
Code and documentation for QSAR modeling and virtual screening of small molecule TMPRSS2 inhibitors

# Installation
1. Clone or download this repo.
2. Install ```conda``` for python3, either through miniconda or anaconda.
3. Run these in terminal
```shell script
#Install dependencies:
cd /path/to/tmprss2
conda env create -f environment.yml
conda env update -n tmprss2 -f chemprop/environment.yml
conda activate tmprss2

#Install chemprop as package
cd ./chemprop_repo
pip install -e .
```

Now, you should have the dependencies needed by our code and by chemprop, and you should
be able to do things like ```from chemprop import ...``` in python scripts.
