__author__ = 'xiaoxiaol'



import pandas as pd
import matplotlib.pylab as plt
from sklearn.metrics import confusion_matrix

def plot_confusion_matrix(cm, title='Confusion matrix', cmap=plt.cm.Blues):
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()

    # tick_marks = np.arange(len(iris.target_names))
    # plt.xticks(tick_marks, iris.target_names, rotation=45)
    # plt.yticks(tick_marks, iris.target_names)
    plt.tight_layout()
    plt.ylabel('x')
    plt.xlabel('y')
    plt.show()




data_DIR = "/data/mat/xiaoxiaol/data/lims2/ivscc_0519"
# output_dir = data_DIR+'/ephys_overlap_clustering_result_pca_aligned'
# df_e=pd.read_csv(output_dir+"/ward_ol_removed/cluster_id.csv")
# df_m=pd.read_csv(data_DIR+"/Ephys.ClusterID.03072016.csv")
#
# df_merge = pd.merge(df_e,df_m, on="specimen_id")
#
# df_merge.to_csv(data_DIR+"/em_types.csv")
df_merge=pd.read_csv(data_DIR+"/em_types_curated.csv")
df_merge['M-cluster']=df_merge['M-cluster'].astype(str)



cm = confusion_matrix(df_merge['M-cluster'], df_merge['E-cluster'])
plot_confusion_matrix(cm)
