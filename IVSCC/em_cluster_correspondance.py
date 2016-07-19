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


import matplotlib.cbook as cbook

# Load a numpy record array from yahoo csv data with fields date,
# open, close, volume, adj_close from the mpl-data/example directory.
# The record array stores python datetime.date as an object array in
# the date column
datafile = cbook.get_sample_data('goog.npy')
try:
    # Python3 cannot load python2 .npy files with datetime(object) arrays
    # unless the encoding is set to bytes. Hovever this option was
    # not added until numpy 1.10 so this example will only work with
    # python 2 or with numpy 1.10 and later
    price_data = np.load(datafile, encoding='bytes').view(np.recarray)
except TypeError:
    price_data = np.load(datafile).view(np.recarray)
price_data = price_data[-250:]  # get the most recent 250 trading days

delta1 = np.diff(price_data.adj_close)/price_data.adj_close[:-1]

# Marker size in units of points^2
volume = (15 * price_data.volume[:-2] / price_data.volume[0])**2
close = 0.003 * price_data.close[:-2] / 0.003 * price_data.open[:-2]

fig, ax = plt.subplots()
ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)

ax.set_xlabel(r'$\Delta_i$', fontsize=20)
ax.set_ylabel(r'$\Delta_{i+1}$', fontsize=20)
ax.set_title('Volume and percent change')

ax.grid(True)
fig.tight_layout()

plt.show()
