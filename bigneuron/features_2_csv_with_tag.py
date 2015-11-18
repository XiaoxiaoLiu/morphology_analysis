# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:00 2015

@author: xiaoxiaol
"""



import utilities.morph_nfb_2_csv as nfb

def generateALLFeatureCSV_gold166(feature_file, feature_csv_file):
    swc_file_nameList, glFeatures, gmiFeatures = readDBFeatures(feature_file)

    allFeatures = np.append(glFeatures, gmiFeatures, 1)
    allColumns = np.append(GL_FEATURE_TAGS, GMI_FEATURE_TAGS, 0)

    df = pd.DataFrame(allFeatures, columns=allColumns)

    df['swc_file'] = pd.Series(swc_file_nameList, index=df.index)

    algorithmList = []
    imageList = []
    for swc_file in swc_file_nameList:
        fn = swc_file.split('/')[-1]
        algorithm = fn.split('v3dpbd_')[-1]
        algorithm = algorithm.split('.')[0]
        if "app1" in algorithm:   # for patterns like *x245_y234_z234_app1.swc
              algorithm = "app1"
        if "app2" in algorithm:
              algorithm = "app2"
        if  "spanningtree" in algorithm: # fastmarching_spanningtree is too long
              algorithm = "spanningtree"
        image = fn.split('.v3dpbd')[0]
        image = image.split('sorted_')[-1]  # for sorted_* swc_files
        algorithmList.append(algorithm)
        imageList.append(image)

    df['algorithm'] = pd.Series(algorithmList, index=df.index)
    df['image'] = pd.Series(imageList, index=df.index)

    allColumns = np.append(np.array(['image', 'algorithm', 'swc_file']), allColumns, 0)

    df = df[allColumns]

    df.to_csv(feature_csv_file, index=False)

    print 'output all feature csv file to :', feature_csv_file
    return




rerun_Gold_standard = 0
if rerun_Gold_standard == 1:

    ########################################## Gold Standard
    data_DIR= "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed"
    #########################################################
    FEATURE_FILE = data_DIR + '/features.nfb'
    nfb.generateALLFeatureCSV_gold166(FEATURE_FILE, data_DIR +'/features_with_tags.csv')




# ##########################################  recon results
data_DIR= "/data/mat/xiaoxiaol/data/20151030_rhea_reconstructions_for_allen300_silver_set"
#/data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted"
#########################################################

FEATURE_FILE = data_DIR + '/features.nfb'

# generate features with meta tags
generateALLFeatureCSV_gold166(FEATURE_FILE, data_DIR +'/features_with_tags.csv')

# generate ano files
nfb.generateLinkerFileFromCSV(data_DIR, data_DIR +'/features_with_tags.csv', 'image')
#




