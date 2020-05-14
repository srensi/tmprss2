from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

from data.data_accessors import tmprss2_to_pandas
import numpy as np
from rdkit import Chem, DataStructs
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier


# collect dataset
dataset = tmprss2_to_pandas()

# visualize dataset
hist, bins, _ = plt.hist(dataset.Activity, bins=8)
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
plt.figure()
plt.hist(dataset.Activity, bins=logbins)
plt.xscale('log')
plt.xlabel('Activity (log scale)')
plt.ylabel('# molecules')
plt.show()

# calculate fingerprints
dataset['fingerprints'] = dataset.SMILES.apply(lambda s: Chem.RDKFingerprint(Chem.MolFromSmiles(s)))

# run the test pieces through knn
X = np.array(dataset.fingerprints.data)
y = np.array(dataset.Activity)

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size=0.2)

# instantiate learning model with k=3
knn = KNeighborsClassifier(n_neighbors=3)
# fitting model
knn.fit(train_X, train_y)
# predict
pred = knn.predict(test_X)
print("Accuracy:{}".format(accuracy_score(test_y, pred)))

# currently euclidean; we should use tanimoto similarity instead
