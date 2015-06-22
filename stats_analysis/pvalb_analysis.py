# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:44:05 2015

@author: xiaoxiaoliu
"""


########################  PVALB vis
pv_tag = merged.IsPVALB
idxs = np.nonzero(pv_tag== 'Y')[0]
pv_features = glFeatures[idxs,:]
no_inds = np.nonzero(pv_tag== 'N')[0]
nonpv_features = glFeatures[no_inds,:]


glFeatures_n= normalizeFeatures(glFeatures)
idx1=np.nonzero(glFeatures_n>3)
glFeatures_n[idx1]=3
idx2=np.nonzero(glFeatures_n<-3)
glFeatures_n[idx2]=-3



#### pca
pca = PCA(n_components=2)
scores =pca.fit(glFeatures_n).transform(glFeatures_n)

pl.scatter(scores[idxs,0],scores[idxs,1],c='y',label='Pvlab',s=50)
pl.scatter(scores[no_inds,0],scores[no_inds,1],c='r',label='other aspiny',s=50)

spiny_tag=merged.dendrite_type_tag
spiny_idxs = np.nonzero(spiny_tag== 'spiny')[0]
aspiny_idxs =np.nonzero(spiny_tag== 'aspiny')[0]
pl.scatter(scores[spiny_idxs,0],scores[spiny_idxs,1],c='b',label='spiny',s=50)
#pl.scatter(scores[aspiny_idxs,0],scores[aspiny_idxs,1],c='c',label='aspiny')
pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()




##################  scatter plot for top two features
pl.figure()
pl.scatter(glFeatures_n[idxs,13],glFeatures_n[idxs,15],c='y',label='Pvlab',s=50 )
pl.scatter(glFeatures_n[no_inds,13],glFeatures_n[no_inds,15],c='r',label='other aspiny',s=50)
pl.scatter(glFeatures_n[spiny_idxs,13],glFeatures_n[spiny_idxs,15],c='b',label='spiny',s=50 )
#pl.scatter(glFeatures_n[aspiny_idxs,13],glFeatures_n[aspiny_idxs,15],c='c',label='aspiny',s=50 )
pl.xlabel('Max Euclidean Path Dist')
pl.ylabel('Max Branch Order')
pl.legend()



#normalized = normalizeFeatures(glFeatures)
#selectFeatures_MRMR(zscores)  # save zscores into a csv file
#plotFeatureVector(glFeatures, "normalized global features")
#plotScoreMatrix(glFeatures, "Similarity Matrix based on global features")
#plotScoreMatrix(gmiFeatures, "Similarity Matrix based on GMI features")


########################  Layer?
pca = PCA(n_components=5)
pv_features=glFeatures_n[idxs,:]
pvalb_all=merged.loc[idxs,:]
pv_scores =pca.fit(pv_features).transform(pv_features)

layer_group = pvalb_all.layer
l23 = np.nonzero(layer_group== 'L2/3')[0]
l4 =np.nonzero(layer_group== 'L4')[0]
l5 =np.nonzero(layer_group== 'L5')[0]
l6 =np.nonzero(layer_group== 'L6')[0]
pl.figure()
pl.scatter(pv_scores[l23,0],pv_scores[l23,1],c='b',label='L2/3',s=50)
pl.scatter(pv_scores[l4,0],pv_scores[l4,1],c='y',label='L4',s=50)
pl.scatter(pv_scores[l5,0],pv_scores[l5,1],c='c',label='L5',s=50)
pl.scatter(pv_scores[l6,0],pv_scores[l6,1],c='r',label='L6',s=50)


pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()

############  region
pca = PCA(n_components=5)
pv_features=glFeatures_n[idxs,:]
pvalb_all=merged.loc[idxs,:]
pv_scores =pca.fit(pv_features).transform(pv_features)

region_group = pvalb_all.region
region1 = np.nonzero(region_group== 'AUDd')[0]
region2 =np.nonzero(region_group== 'VISl')[0]
region3 =np.nonzero(region_group== 'VISp')[0]
region4 =np.nonzero(region_group== 'VISal')[0]
pl.figure()
pl.scatter(pv_scores[region1,0],pv_scores[region1,1],c='b',label='AUDd',s=50)
pl.scatter(pv_scores[region2,0],pv_scores[region2,1],c='y',label='VIS1',s=50)
pl.scatter(pv_scores[region3,0],pv_scores[region3,1],c='c',label='VISp',s=50)
pl.scatter(pv_scores[region4,0],pv_scores[region4,1],c='r',label='VISa1',s=50)


pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()

############  cell shape
pca = PCA(n_components=5)
pv_features=glFeatures_n[idxs,:]
pvalb_all=merged.loc[idxs,:]
pv_scores =pca.fit(pv_features).transform(pv_features)

region_group = pvalb_all.cell_shape
region1 = np.nonzero(region_group== 'asymm')[0]
region2 =np.nonzero(region_group== 'symm')[0]
region3 =np.nonzero(region_group== 'bitufted')[0]

pl.figure()
pl.scatter(pv_scores[region1,0],pv_scores[region1,1],c='b',label='asymm',s=50)
pl.scatter(pv_scores[region2,0],pv_scores[region2,1],c='y',label='symm',s=50)
pl.scatter(pv_scores[region3,0],pv_scores[region3,1],c='c',label='bifuted',s=50)



pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()

