import pdb
import numpy as np
import pandas as pd
import pdb
from sklearn import decomposition
from sklearn import mixture


nikpay = pd.read_csv("super_diff.nikpay.to_plot.gz", delimiter = ' ', header = None)
christophersen =  pd.read_csv("super_diff.christophersen.to_plot.gz", delimiter = ' ', header = None)
mahajan =  pd.read_csv("super_diff.mahajan.to_plot.gz", delimiter = ' ', header = None)
michailidou =  pd.read_csv("super_diff.michailidou.to_plot.gz", delimiter = ' ', header = None)
liu =  pd.read_csv("super_diff.liu.to_plot.gz", delimiter = ' ', header = None)
lambert =  pd.read_csv("super_diff.lambert.to_plot.gz", delimiter = ' ', header = None)
lopezisac =  pd.read_csv("super_diff.lopez-isac.to_plot.gz", delimiter = ' ', header = None)
malik = pd.read_csv("super_diff.malik.to_plot.gz", delimiter = ' ', header = None)

nikpay = nikpay.to_numpy()
christophersen = christophersen.to_numpy()
mahajan = mahajan.to_numpy()
michailidou = michailidou.to_numpy()
liu = liu.to_numpy()
lambert = lambert.to_numpy()
lopezisac = lopezisac.to_numpy()
malik = malik.to_numpy()

superm = np.hstack((nikpay, christophersen,mahajan,michailidou,liu,lambert,lopezisac,malik))
super_mtx = np.transpose(superm)



pca = decomposition.PCA(n_components=5)
pca.fit(super_mtx)
super_PCA = pca.transform(super_mtx)


print("PCA explained variance \n", pca.explained_variance_ratio_)

np.savetxt('PCA_components.super.to_plot.gz', super_PCA)
