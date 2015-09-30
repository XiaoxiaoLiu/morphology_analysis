
######  iterative PCA ##############
from rpy2.robjects.packages import importr
R_base     = importr('base')
R_stats    = importr('stats')
R_graphics = importr('graphics')
R_gplots =  importr('gplots')

from rpy2.robjects import pandas2ri
pandas2ri.activate()

def binary_partition(df_data, ):
    df_left = pd.DataFrame()
    df_right = pd.DataFrame()

    r_df_data = pandas2ri(df_data)
    pca = R_stats.princomp(r_df_data)

    return df_left, df_right

def iterative_pca(df_all,feature_names,num_clusters):
    df_zscores = get_zscore_features(df_all, feature_names, None, 1)
    df_left, df_right = binary_partition(df_zscores)

    #iterative_pca(df_left, features_names, )


    return
############################################################################################
#############################################################################################