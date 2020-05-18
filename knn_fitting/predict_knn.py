from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

from data.data_accessors import tmprss2_to_pandas
import numpy as np
from rdkit import Chem, DataStructs
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors, KNeighborsRegressor

# collect dataset
dataset = tmprss2_to_pandas()

# visualize dataset
plt.figure()
hist, bins, _ = plt.hist(dataset.Activity, bins=8)
plt.clf()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
plt.hist(dataset.Activity, bins=logbins)
plt.xscale('log')
plt.xlabel('Activity (log scale)')
plt.ylabel('# molecules')
plt.show()

# calculate fingerprints
dataset['fingerprints'] = dataset.SMILES.apply(lambda s: Chem.RDKFingerprint(Chem.MolFromSmiles(s)))
# convert fingerprints to arrays
X = np.zeros([len(dataset), 2048])
for i in range(len(dataset)):
    DataStructs.ConvertToNumpyArray(dataset.fingerprints.iloc[i], X[i])

# run the test pieces through knn
y = dataset.Activity.values

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size=0.2)


def tanimoto_similarity(fp1, fp2):
    """from http://infochim.u-strasbg.fr/CS3_2014/Slides/CS3_2014_Willett.pdf"""
    c = np.sum((fp1 == fp2) & (fp1 == 1))  # bits set in common
    a = np.sum(fp1)  # bits set in fp1
    b = np.sum(fp2)  # bits set in fp2
    if (a + b - c) == 0:
        print(a, b, c)
    return c / (a + b - c)


rmse = []
ks = np.arange(1, 10)
for k in ks:
    nbrs = KNeighborsRegressor(n_neighbors=k, algorithm='auto',
                               metric=tanimoto_similarity)
    nbrs.fit(train_X, train_y)
    pred = nbrs.predict(train_X)

    rmse.append(mean_squared_error(train_y, pred))

plt.figure()
plt.plot(ks, rmse)
