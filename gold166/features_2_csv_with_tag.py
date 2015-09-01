# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:00 2015

@author: xiaoxiaol
"""




import numpy as np
import pandas as pd
from os import sys, path

p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
sys.path.append(p+'/utilities')

import morph_nfb_2_csv as nfb





########################################## Gold Standard
data_DIR= "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed"
#########################################################
FEATURE_FILE = data_DIR + '/features.nfb'

nfb.generateALLFeatureCSV(FEATURE_FILE, data_DIR +'/features_with_tags.csv')



#
# ##########################################  recon results
# data_DIR= "/data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted"
# #########################################################
#
# FEATURE_FILE = data_DIR + '/features.nfb'
#
# nfb.generateALLFeatureCSV(FEATURE_FILE, data_DIR +'/features_with_tags.csv')
# nfb.generateLinkerFileFromCSV(data_DIR, data_DIR +'/features_with_tags.csv', 'image')
#




