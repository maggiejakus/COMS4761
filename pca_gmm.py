import pdb
import numpy as np
import pandas as pd
import pdb
from sklearn import decomposition
from sklearn import mixture


nikpay = pd.read_csv("gmm_diff.nikpay.to_plot.gz", delimiter = ' ', header = None)
christophersen =  pd.read_csv("gmm_diff.christophersen.to_plot.gz", delimiter = ' ', header = None) 
mahajan =  pd.read_csv("gmm_diff.mahajan.to_plot.gz", delimiter = ' ', header = None)
michailidou =  pd.read_csv("gmm_diff.michailidou.to_plot.gz", delimiter = ' ', header = None)
liu =  pd.read_csv("gmm_diff.liu.to_plot.gz", delimiter = ' ', header = None)
lambert =  pd.read_csv("gmm_diff.lambert.to_plot.gz", delimiter = ' ', header = None)
lopezisac =  pd.read_csv("gmm_diff.lopez-isac.to_plot.gz", delimiter = ' ', header = None)
malik = pd.read_csv("gmm_diff.malik.to_plot.gz", delimiter = ' ', header = None)

nikpay = nikpay.to_numpy()
christophersen = christophersen.to_numpy()
mahajan = mahajan.to_numpy()
michailidou = michailidou.to_numpy()
liu = liu.to_numpy()
lambert = lambert.to_numpy()
lopezisac = lopezisac.to_numpy()
malik = malik.to_numpy()

gmm = np.hstack((nikpay, christophersen,mahajan,michailidou,liu,lambert,lopezisac,malik))
gmm_mtx = np.transpose(gmm)



pca = decomposition.PCA(n_components=5)
pca.fit(gmm_mtx)
gmm_PCA = pca.transform(gmm_mtx)


print("PCA explained variance \n", pca.explained_variance_ratio_)

np.savetxt('PCA_components.GMM.to_plot.gz', gmm_PCA)




