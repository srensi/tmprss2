# tmprss2
Code and documentation for QSAR modeling and virtual screening of small molecule TMPRSS2 inhibitors

# Installation
Install conda, either through miniconda or anaconda.

Clone or download this repo.

Navigate to the tmprss2 directory.

Install dependencies:
```python
conda env create -f environment.yml
conda env update -n tmprss2 -f chemprop/environment.yml
conda activate tmprss2
```

Install chemprop as a package:
```
cd /path/to/tmprss2/chemprop
pip install -e .
```

Hopefully, good to go.
