import pandas as pd
import glob
import os

from chemprop.args import TrainArgs
from chemprop.train import cross_validate
from chemprop.utils import create_logger

csvs = glob.glob(os.path.join('data', '*.csv'))
raw_data = pd.concat((pd.read_csv(f) for f in csvs))

chemprop_data = raw_data[['SMILES', 'Activity']]
chemprop_data.to_csv('chemprop_in.csv', index=False)

# argument passing pretty janky but it's set up to use command line
args = TrainArgs().parse_args(['--data_path', 'chemprop_in.csv', '--dataset_type', 'regression',
                               '--save_dir', 'models'])
logger = create_logger(name='train', save_dir=args.save_dir, quiet=args.quiet)
cross_validate(args, logger)
