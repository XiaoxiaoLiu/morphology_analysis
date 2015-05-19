import numpy as np
import matplotlib.pylab as pl
from scipy.stats import norm
import pandas as pd
from sklearn.lda import LDA
import scipy.cluster.hierarchy as hac
from mpl_toolkits.mplot3d import Axes3D


#===================== DATA PATH ==============================================

data_DIR = "/Users/xiaoxiaoliu/work/data/lims2/manual"
fn = data_DIR+'/zscore_glFeatures_subtypes.csv'
#fn = data_DIR+'/zscore_ei_glFeatures.csv'

#==============================================================================

def compareHist(data1, data2, _title):
    pl.figure()
    pl.show()
    pl.hist(data1, bins=20,normed=True, alpha=0.6, color='b')
    pl.hist(data2, bins=20,normed=True, alpha=0.6, color='r')

    # Fit a normal distribution to the data:
    mu1, std1 = norm.fit(data1)
    xmin, xmax = pl.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu1, std1)
    pl.plot(x, p, 'k', linewidth=2, color='b')

    # Fit a normal distribution to the data:
    mu2, std2 = norm.fit(data2)
    xmin, xmax = pl.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu2, std2)
    pl.plot(x, p, 'k', linewidth=2, color='r')
    pl.title(_title)
    pl.savefig(data_DIR + '/'+ _title + '.png',bbox_inches='tight')
    #pl.close()
    return



def singleFeatureComparison(featureId, label1, label2):
    pos = (np.array(np.where(label==label1)) -1 ).flatten() # inhibitory
    inh = raw_features[pos, featureId-1]

    pos = (np.array(np.where(label==label2)) -1 ).flatten()
    exc = raw_features[pos, featureId-1]

    compareHist(inh, exc,featureNames[featureId])
    return



def twoFeatureScatterPlotMultiClusters(featureId1, featureId2, class_id_list,  fig_title):
    class_f1=[]
    class_f2=[]
    colors=['blue','red','green','yellow','black','cyan','magenta','gray']
    fig = pl.figure()
    pl.show()
    i=0
    handles=[]
    for  LABEL in class_id_list:
	    pos = (np.array(np.where(label== LABEL)) ).flatten() # inhibitory
	    class_f1 = raw_features[pos, featureId1-1]
	    class_f2 = raw_features[pos, featureId2-1]
	    pl.scatter( class_f1, class_f2 ,s=50, c=colors[i] )
	    i=i+1



    pl.title(fig_title)
    pl.xlabel(featureNames[featureId1])
    pl.ylabel(featureNames[featureId2])
    fig.legend(class_id_list, 'upper left')
    pl.savefig(data_DIR + '/'+ fig_title + '.png',bbox_inches='tight')

    #pl.close()
    return




def twoFeatureScatterPlot(featureId1,featureId2, fig_title):
    pos = (np.array(np.where(label==1)) -1 ).flatten() # inhibitory
    inh1 = raw_features[pos, featureId1-1]
    inh2 = raw_features[pos, featureId2-1]

    pos = (np.array(np.where(label==2)) -1 ).flatten()
    exc1 = raw_features[pos, featureId1-1]
    exc2 = raw_features[pos, featureId2-1]

    pl.figure()
    pl.show()
    pl.scatter(inh1,inh2,s=50,c=u'b')
    pl.scatter(exc1,exc2,s=50,c=u'r')
    pl.title(fig_title)
    pl.xlabel(featureNames[featureId1])
    pl.ylabel(featureNames[featureId2])
    pl.savefig(data_DIR + '/'+ fig_title + '.png',bbox_inches='tight')
    #pl.close()
    return


#==============================================================================

def lda(X,y):
    lda = LDA(n_components=3)
    X_r2 = lda.fit(X,y).transform(X)

    plt.figure()
    for c, i, target_name in zip("gbr", [0,1, 2], ['others','inhibitory','excitatory']):
        plt.scatter(X_r2[y == i, 0], X_r2[y == i, 1], c=c, label=target_name)
        plt.legend()
        plt.title('LDA')


def pca(X):
    pca = PCA(n_components=2)
    X_r = pca.fit(X).transform(X)

    plt.figure()
    for c, i, target_name in zip("gbr", [0,1, 2], ['others','inhibitory','excitatory']):
        plt.scatter(X_r[y == i, 0], X_r[y == i, 1], c=c, label=target_name)
        plt.legend()
        plt.title('PCA')


def hcluster(X):
    z = hac.linkage(X)
    d = hac.dendrogram(z)



#=====================       MAIN   ==========================================
df=pd.read_csv(fn, sep=',',header=0)
featuresTB = df.values
r,c = featuresTB.shape


label = featuresTB[:,0] # 0- N/A  1-inhibitory  2-excitatory
featureNames = df.columns
raw_features = featuresTB[0:r,1:c]
print raw_features.shape
pl.figure()
pl.boxplot(raw_features)
pl.xlabel('feature id')
pl.ylabel('zscores')
pl.show()



#=======================  Feature Sele tion ===========================================

# top 5 features
singleFeatureComparison(14,1,7)
singleFeatureComparison(3,2,4)

#twoFeatureScatterPlotMultiClusters(14,3,[1,2,4,6,7],'Top 2 Features')

featureId1 = 14
featureId2 = 3
featureId3 = 4
class_id_list=[1,2,4,6,7,3,5,8]
label_names=['Tall-tufted','Star-pyramid','Simple-tufted','Bi-tufted','Basket-cell','CorticoCortical','Martinotti','CorticalThamalic']
colors=['b', 'c', 'y', 'm', 'r','yellow','black','red']
markers=['x','o','x','o','x','*','*','*']
class_f1=[]
class_f2=[]
pl.figure()
pl.show()
i=0
ax = pl.subplot(111, projection='3d')

for  LABEL in class_id_list:
	    pos = (np.array(np.where(label== LABEL)) ).flatten() # inhibitory
	    class_f1 = raw_features[pos, featureId1-1]
	    class_f2 = raw_features[pos, featureId2-1]
	    class_f3 = raw_features[pos, featureId3-1]
	    ax.plot(class_f1, class_f2,class_f3, markers[i], color=colors[i], label= label_names[i])
	    i=i+1

pl.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8, bbox_to_anchor=(0, 0))

pl.show()



#twoFeatureScatterPlot(14,7,'Top two features scatter plot')

# featureId1 =7
# featureId2 =8
# pos = (np.array(np.where(label==1)) -1 ).flatten() # inhibitory
# inh1 = raw_features[pos, featureId1-1]
# inh2 = raw_features[pos, featureId2-1]
# pos = (np.array(np.where(label==2)) -1 ).flatten()
# exc1 = raw_features[pos, featureId1-1]
# exc2 = raw_features[pos, featureId2-1]
#
# compareHist(inh2/inh1, exc2/exc1, 'ratio')



#=====================   CLUSTERING ============================================
X=featuresTB[0:r,1:c]
y=(label).astype(int)

#lda(X,y)
#pca(X)
#hcluster(X)


