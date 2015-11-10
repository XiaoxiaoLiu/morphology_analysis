# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:00 2015

@author: xiaoxiaol
"""



import utilities.morph_nfb_2_csv as nfb



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
nfb.generateALLFeatureCSV_gold166(FEATURE_FILE, data_DIR +'/features_with_tags.csv')

# generate ano files
nfb.generateLinkerFileFromCSV(data_DIR, data_DIR +'/features_with_tags.csv', 'image')
#




