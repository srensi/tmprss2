import pandas as pd
import glob
import os

csvs = glob.glob(os.path.join('data', '*.csv'))
df = pd.concat((pd.read_csv(f) for f in csvs))

