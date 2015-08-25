# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:00 2015

@author: xiaoxiaol
"""

import numpy as np
import pandas as pd
from os import sys, path




GL_FEATURE_TAGS = np.array(['num_nodes', 'soma_surface', 'num_stems','num_bifurcations', 'num_branches', 'num_of_tips',  'overall_width', 'overall_height',  'overall_depth', 'average_diameter',    'total_length', 'total_surface', 'total_volume', 'max_euclidean_distance',       'max_path_distance', 'max_branch_order',  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])
GMI_FEATURE_TAGS = np.array(['moment1', 'moment2', 'moment3','moment4', 'moment5', 'moment6',  'moment7', 'moment8',  'moment9', 'moment10',    'moment11', 'moment12', 'moment13', 'avgR'])




#===================================================================
def readDBFeatures(FEATURE_FILE):
    # TODO: detect nan values
    glf_featureList = []  # each row is a feature vector
    gmi_featureList = []
    swc_file_nameList = []
    with open (FEATURE_FILE,'r') as  f:
        for fn_line in f: # ignore the SWCFILE=* line
                         
            swc_file= fn_line[8:].strip()
            swc_file_nameList.append(swc_file)
            
            line_globalFeature = (f.next()).strip()
            glf = map(float,line_globalFeature.split('\t'))
            glf_featureList.append(glf)

            line_GMI = (f.next()).strip()
            gmi = map(float,line_GMI.split('\t'))
            gmi_featureList.append(gmi)

    return  swc_file_nameList, np.array(glf_featureList), np.array(gmi_featureList)




def  generateALLFeatureCSV(feature_file, feature_csv_file):
    
      swc_file_nameList, glFeatures, gmiFeatures = readDBFeatures(feature_file)

    # attacheh specimen id, nrrd id to the tabe
  
    allFeatures = np.append(glFeatures,gmiFeatures,1)
    allColums = np.append(gl_feature_names,gmi_feature_names,0)
    df = pd.DataFrame(allFeatures,  columns = allColums )
    
    df['swc_file'] = pd.Series(swc_file_nameList, index=df.index)

    df.to_csv(feature_csv_file, index=False)

 
    print 'output all feature csv file to :',feature_csv_file
    return





#===================================================================
def zscore(features):
    zscores = stats.zscore(features,0)
    return zscores


def normalizeFeatures(features):
    meanFeatures = np.mean(features,0)
    stdFeatures = np.std(features, 0)
    if np.count_nonzero(stdFeatures)< len(stdFeatures):
          print "zero detected"
          print stdFeatures
    normalized = (features - meanFeatures)/stdFeatures
    return normalized








def concatCSVs(csv1, csv2, outcsv):
    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)

    #out_df = pd.merge(df1, df2)
    out_df = pd.concat([df1,df2], axis=1)
    out_df.to_csv(outcsv, index=False)
    return

def saveScoreMatrix(featureArray,scoreMatrix_file, REMOVE_OUTLIER=1):
    scoreMatrix =[]
    normalized = zscore(featureArray)

    # remove outliers!!!
    normalized[normalized < -3]  =-3
    normalized[normalized > 3] = 3

    for i in range(len(normalized)):
        queryFeature = normalized[i] # each row is a feature vecto
        #scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)
        scores = np.sum(np.abs(normalized-queryFeature)**2,1)
        scoreMatrix.append(scores)

    df = pd.DataFrame(scoreMatrix)
    df.to_csv(scoreMatrix_file)


def  generateGlFeatureCSV(new_csv_file):
    gl_feature_csv_file	= data_DIR+"/glFeatures.csv"  # mid result
    # attacheh specimen id, nrrd id,
    glfFeatures, gmiFeatures = readDBFeatures(FEATURE_FILE)


    #attach  specimen_name
    list_df = pd.read_csv(LIST_CSV_FILE)
    sp_ids = list_df.specimen_id
    sp_names = np.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    df = pd.DataFrame(glfFeatures, columns = gl_feature_names)
    df['specimen_name'] = pd.Series(sp_names, index=df.index)
    #df = df['specimen_name',feature_names]
    df.to_csv(gl_feature_csv_file, index=False)

    concatCSVs(LIST_CSV_FILE,gl_feature_csv_file, new_csv_file)
    return

def  generateGmiFeatureCSV(new_csv_file):
    gmi_feature_csv_file= data_DIR+"/gmiFeatures.csv" # mid result

    # attacheh specimen id, nrrd id,
    glFeatures, gmiFeatures = readDBFeatures(FEATURE_FILE)

    #attach  specimen_name
    list_df = pd.read_csv(LIST_CSV_FILE)
    sp_ids = list_df.specimen_id
    sp_names = np.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    df = pd.DataFrame(gmiFeatures, columns = gmi_feature_names)
    df['specimen_name'] = pd.Series(sp_names, index=df.index)
    #df = df['specimen_name',feature_names]
    df.to_csv(gmi_feature_csv_file, index=False)

    concatCSVs(LIST_CSV_FILE,gmi_feature_csv_file, new_csv_file)
    return


def generateLinkerFileFromCSV(result_dir, csvfile, column_name):
	df = pd.read_csv(csvfile)
	types = df[column_name]
	for atype in np.unique(types):
		idxs = np.nonzero(types==atype)[0]
		swc_files = df['orca_path']
		with open(result_dir+'/'+atype+'.ano','w') as outf:
        	   for afile in swc_files[idxs]:
                       filename = afile.split('/')[-1]
                       line='SWCFILE='+filename+'\n'
                       outf.write(line)
                   outf.close()


def generateFeatureMergedCSV(outFile):
    all_feature_csv_with_id_file = data_DIR+"/allFeatures_withid.csv"
    generateALLFeatureCSV( all_feature_csv_with_id_file)

    df_complete = pd.read_csv(all_feature_csv_with_id_file)
    mycolumns = np.array(['specimen_id','specimen_name','nrid','orca_path'])
    mycolumns = np.append(mycolumns,gl_feature_names,0)
    mycolumns = np.append(mycolumns,gmi_feature_names,0)
    df_complete = df_complete.reindex(columns=mycolumns)
    print df_complete.columns

        # merge all info
    df_type = pd.read_csv(data_DIR+'/custom_report-IVSCC_classification-April_2015.csv')
    merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
    merged.to_csv(outFile)






#==================================================================================================
def main():
   
    WORK_PATH="/Users/xiaoxiaoliu/work"
     
    ########################################## data dir
    #data_DIR= WORK_PATH+"/data/lims2/neuron_recon_2"
    data_DIR= "/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned"
    LIST_CSV_FILE =  data_DIR+'/list.csv'
    #########################################################
    
    
    data_linker_file =  data_DIR+'/original/mylinker.ano'
    FEATURE_FILE = data_DIR + '/preprocessed/prep_features.nfb'
    readDBFeatures(FEATURE_FILE)
    
    generateALLFeatureCSV(FEATURE_FILE)
    
    


if __name__ == "__main__":
      main()