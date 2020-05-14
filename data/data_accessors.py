"""
Houses scripts to convert raw data into more malleable objects
"""

import glob
import os
import pandas as pd

data_folder = os.path.dirname(__file__)


def tmprss2_to_pandas():
    """collects the tmprss2 assay data into a pandas dataframe.
        currently incomplete, only uses csvs, doesn't validate between csv and xlsx datasets"""
    csvs = glob.glob(os.path.join(data_folder, 'tmprss2_meyer_et_al', '*.csv'))
    return pd.concat((pd.read_csv(f) for f in csvs), ignore_index=True)
