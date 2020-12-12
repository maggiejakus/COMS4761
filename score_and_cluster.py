import pdb
import numpy as np
import pandas as pd
import pdb
from sklearn import decomposition
from sklearn import mixture
#import mca

author = "malik"


#annotations file - converting alkes grch38 to grch37
#pdb.set_trace()
mich = pd.read_csv("converted_annot.txt", delimiter = ' ')

#rsIDs of superclump SNPs - this time, everything is bpPOS
super_snps = pd.read_csv("superclump.snps."+author+".to_plot")


#confused about functionality..
#import PGS

## MALIK PGS doesn't have BP POS - using RSID instead
PGS = pd.read_csv("PGS000038.txt", delimiter = '\t', comment='#')
PGS_snps = PGS.iloc[:, 0].to_numpy() #PGS BP POS LOCATION

#pdb.set_trace()
#mich is annotations for cols 14 onward
annot_names = np.array(list(mich)[14:])
annot = mich.iloc[:,14:].to_numpy()

#only take annotations with "extend" in name
annot = annot[:,["extend" in x for x in annot_names]]
annot_names = annot_names[["extend" in x for x in annot_names]]

#identify values in mich
#THESE VALUES CHANGE W DIFF GWAS
bp_pos = mich.iloc[:,1].to_numpy()
beta = (mich.iloc[:,6]).to_numpy() #NOTE - used LOG OF ODDS RATIO FOR PHELAN
pval_weight = np.log10(mich.iloc[:,7].to_numpy()) * -1.0
mich = mich.iloc[:,2].to_numpy() #this is rsID

#scale annotations by beta vals
annot = annot*np.abs(beta[:,np.newaxis])
#annot = annot*np.abs(beta[:,np.newaxis])*pval_weight[:,np.newaxis]
#pdb.set_trace()

#if doing R logistic PCA - didn't give much
#np.savetxt('annot_for_r.txt', annot)
#logpca = pd.read_csv("logpca_pcs.txt", delimiter = ' ', header = None)
#ss_PCA = logpca.to_numpy()

#if doing normal PCA - fit using annotations scoes
pca = decomposition.PCA(n_components=10)
pca.fit(annot)
ss_PCA = pca.transform(annot)
print(pca.explained_variance_ratio_)

#noticed a few extreme outliers - remove them
bp_pos = bp_pos[ss_PCA[:,0] < np.quantile(ss_PCA[:,0], [0.999])]
mich = mich[ss_PCA[:,0] < np.quantile(ss_PCA[:,0], [0.999])]
annot = annot[ss_PCA[:,0] < np.quantile(ss_PCA[:,0], [0.999]),:]
#annot_weight = annot_weight[ss_PCA[:,0] < np.quantile(ss_PCA[:,0],[0.999]),:]
ss_PCA = ss_PCA[ss_PCA[:,0] < np.quantile(ss_PCA[:,0], [0.999]), :]

#calculate variance explained by PCA
explained_variance = pca.explained_variance_ratio_

#calculate functional annotation scores for PGS

#pdb.set_trace()
### ~~~!!!!
bool_PGS = np.in1d(mich, PGS_snps) #make sure this makes sense
PGS_annot = annot[bool_PGS, :]
PGS_weights = np.sum(PGS_annot, axis = 0)


PGS_weights = PGS_weights / np.sum(PGS_weights)
PGS_vec = np.vstack((annot_names, PGS_weights))
PGS_not = [not(x) for x in bool_PGS]
PGS_nonrisky = annot[PGS_not, :]
PGS_nonrisky_weights = np.sum(PGS_nonrisky, axis = 0)
PGS_nonrisky_weights = PGS_nonrisky_weights / np.sum(PGS_nonrisky_weights)
PGS_nonrisky_vec = np.vstack((annot_names, PGS_nonrisky_weights))
PGS_diff = [float(PGS_nonrisky_vec[1][i]) - float(PGS_vec[1][i]) for i in range(len(PGS_vec[1]))]

#calculate functional annotation scores for superclump

#pdb.set_trace()
super_snps = super_snps.to_numpy()[:,0]
bool_super = np.in1d(mich, super_snps) #CHANGE HERE
super_annot = annot[bool_super, :]
super_weights = np.sum(super_annot, axis = 0)
super_weights = super_weights / np.sum(super_weights)
super_vec = np.vstack((annot_names, super_weights))

super_nonrisky_idx = [not(x) for x in bool_super]
super_nonrisky = annot[super_nonrisky_idx, :]
super_nonrisky_weights = np.sum(super_nonrisky, axis = 0)
super_nonrisky_weights = super_nonrisky_weights / np.sum(super_nonrisky_weights)
super_nonrisky_vec = np.vstack((annot_names, super_nonrisky_weights))
super_diff = [float(super_nonrisky_vec[1][i]) - float(super_vec[1][i]) for i in range(len(super_nonrisky_vec[1]))]

#initialize GMM with superclump
guess_means = np.vstack((np.mean(ss_PCA[np.invert(bool_super),:], axis=0), np.mean(ss_PCA[bool_super,:], axis=0)))

#use GMMs to find risky vs nonrisky SNPs
gmm = mixture.GaussianMixture(n_components=2, covariance_type='full', means_init = guess_means)
gmm.fit(ss_PCA)

probabilities = gmm.predict_proba(ss_PCA)


np.savetxt('PGS_annot_sums.'+author+'.to_plot.gz', PGS_vec, '%s')
np.savetxt('PCA_components.'+author+'.to_plot.gz', ss_PCA)
np.savetxt('probabilities.'+author+'.to_plot.gz', probabilities)
np.savetxt('means.'+author+'.to_plot.gz', gmm.means_)
np.savetxt('weights.'+author+'.to_plot.gz', gmm.weights_)
np.savetxt('pca_rsids.'+author+'.to_plot.gz', mich, '%s')
np.savetxt('annot.'+author+'.to_plot.gz', annot)
np.savetxt('annot_names.'+author+'.to_plot.gz', annot_names, '%s')
np.savetxt('super_sums.'+author+'.to_plot.gz', super_vec, '%s')
np.savetxt('super_nonrisky.'+author+'.to_plot.gz', super_nonrisky_vec, '%s')

np.savetxt('super_risky.ttest.'+author+'.to_plot.gz', super_annot, '%s')
np.savetxt('super_nonrisky.ttest.'+author+'.to_plot.gz', super_nonrisky, '%s')

np.savetxt('PGS_risky.ttest.'+author+'.to_plot.gz', PGS_annot, '%s')
np.savetxt('PGS_nonrisky.ttest.'+author+'.to_plot.gz', PGS_nonrisky, '%s')


#super_probs = probabilities[bool_super, :]

#calculate GMM functional annotation scores
nonrisky = annot[probabilities[:,0] > 0.5, :]
risky = annot[probabilities[:,0] <= 0.5, :]
nonrisky_weights = np.sum(nonrisky, axis = 0)
nonrisky_weights = nonrisky_weights / np.sum(nonrisky_weights)
risky_weights = np.sum(risky, axis = 0)
risky_weights = risky_weights / np.sum(risky_weights)


total_weights = np.vstack((nonrisky_weights, risky_weights))
total_vec = np.vstack((annot_names, total_weights))
gmm_diff = [(nonrisky_weights[i]) - (risky_weights[i]) for i in range(len(risky_weights))]
np.savetxt('gmm_diff.'+author+'.to_plot.gz', gmm_diff)
np.savetxt('super_diff.'+author+'.to_plot.gz', super_diff)
np.savetxt('PGS_diff.'+author+'.to_plot.gz', PGS_diff)
np.savetxt('total_weights.'+author+'.to_plot.gz', total_vec, '%s')

#export nonrisky and risky GMM annotations in whole:
np.savetxt('gmm_nonrisky.ttest.'+author+'.to_plot.gz', nonrisky)
np.savetxt('gmm_risky.ttest.'+author+'.to_plot.gz', risky)


print("done")



