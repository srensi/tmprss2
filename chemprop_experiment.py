import pandas as pd
import glob
import os

csvs = glob.glob(os.path.join('data', '*.csv'))
raw_data = pd.concat((pd.read_csv(f) for f in csvs))

chemprop_data = raw_data[['SMILES', 'Activity']]
chemprop_data.to_csv('chemprop_in.csv', index=False)

os.system('python chemprop/train.py --data_path chemprop_in.csv --dataset_type regression --save_dir models')