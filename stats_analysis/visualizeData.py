import os
import scipy.stats
import numpy
import matplotlib.pylab as pl
import pandas as pd

from os import sys, path
p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
sys.path.append(p+'/utilities')
import lims_orca_utils as lu

# program path on this machine
#===================================================================
MRMR= "/Users/xiaoxiaoliu/work/src/mrmr_c_src/mrmr"

V3D="/local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d"

# data dir
data_DIR= "/home/xiaoxiaol/work/data/lims2/neuron_recon_2"

list_csv_file =  data_DIR+'/list.csv'
data_linker_file =  data_DIR+'/original/mylinker.ano' 
preprocessed_data_linker_file = data_DIR+'/preprocessed/mylinker.ano'
feature_file = data_DIR + '/preprocessed/prep_features.nfb'



#===================================================================
def zscore(features):
    zscores=scipy.stats.mstats.zscore(features,0)
    return zscores


def normalizeFeatures(features):
    meanFeatures = numpy.mean(features,0)
    stdFeatures = numpy.std(features, 0)
    normalized = (features - meanFeatures)/stdFeatures

     ## hackish :rearrange data
#    tmp8= normalized[8]
#    tmp9= normalized[9]
#    tmp31= normalized[31]
#
#    normalized=numpy.delete(normalized,8,0)
#    normalized=numpy.delete(normalized,9-1,0)
#    normalized=numpy.delete(normalized,31-2,0)
#
#    normalized=numpy.append(normalized,[tmp31],axis=0)
#    normalized=numpy.append(normalized,[tmp9],axis=0)
#    normalized=numpy.append(normalized,[tmp8],axis=0)
#
    return normalized


def plotFeatureVector(featureArray,fig_title):
    normalized = normalizeFeatures(featureArray)
    pl.figure()
    pl.imshow(normalized,interpolation='none')
    pl.colorbar()
    pl.title(fig_title)
    pl.xlabel('feature ID')
    pl.ylabel('neuron ID')
    pl.show()
    plt.savefig(data_DIR+'/'+fig_title+'.png', bbox_inches="tight")


def plotScoreMatrix(featureArray, fig_title):
    scoreMatrix =[]
    normalized = normalizeFeatures(featureArray)

    for i in range(len(normalized)):
        queryFeature = normalized[i] # each row is a feature vecto
        scores = numpy.exp(-numpy.sum(abs(normalized-queryFeature)**2,1)/100)
        #scores = numpy.sum(abs(normalized-queryFeature)**2,1)
        scoreMatrix.append(scores)

    pl.figure()
    pl.imshow(scoreMatrix,interpolation='none')
    pl.colorbar()
    pl.title(fig_title)
    pl.xlabel('neuron ID')
    pl.ylabel('neuron ID')
    pl.show()
    plt.savefig(data_DIR+'/'+fig_title+'.png', bbox_inches="tight")


def selectFeatures_MRMR(featureArray, threshold=0, number_of_features=5, selection_method='MID'):
    #write out feature array into a csv file, then execute MRMR
    csvfile = data_DIR+"/zscore_glFeatures.csv"

    numpy.savetxt(csvfile, featureArray, delimiter=",")

    # call MRMR
    # cmd = MRMR +  " -i "+ csvfile + " -t "+ threshold + " -n " + number_of_features
    # print cmd
    #os.system(cmd)
    return

def readDBFeatures(feature_file):
    # TODO: detect nan values
    glf_featureList = []  # each row is a feature vector
    gmi_featureList = []
    with open (feature_file,'r') as  f:
        for fn_line in f: # ignore the SWCFILE=* line
            line_globalFeature = (f.next()).strip()
            glf = map(float,line_globalFeature.split('\t'))
            glf_featureList.append(glf)

            line_GMI = (f.next()).strip()
            gmi = map(float,line_GMI.split('\t'))
            gmi_featureList.append(gmi)

    return  numpy.array(glf_featureList), numpy.array(gmi_featureList)



def concatCSVs(csv1, csv2, outcsv):

    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)

    #out_df = pd.merge(df1, df2)
    out_df = pd.concat([df1,df2], axis=1)

    out_df.to_csv(outcsv, index=False)

    return



def  generateCompleFeatureCSV(new_csv_file):
    # attacheh specimen id, nrrd id,
    glfFeatures, gmiFeatures = readDBFeatures(feature_file)
    feature_names=['num_nodes',	'soma_surface',	'num_stems','num_bifurcations',	'num_branches',	'num_of_tips',	'overall_width',	'overall_height',	'overall_depth',
    'average_diameter',	'total_length',	'total_surface',	'total_volume',	'max_euclidean_distance',	'max_path_distance',	'max_branch_order',	'average_contraction',
    'average fragmentation',	'parent_daughter_ratio',	'bifurcation_angle_local',	'bifurcation_angle_remote']

    feature_csv_file = data_DIR+"/glFeatures.csv"

    #attach  specimen_name
    list_df = pd.read_csv(list_csv_file)
    sp_ids = list_df.specimen_id
    sp_names = numpy.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    df = pd.DataFrame(glfFeatures, columns = feature_names)
    df['specimen_name'] = pd.Series(sp_names, index=df.index)
    #df = df['specimen_name',feature_names]
    df.to_csv(feature_csv_file, index=False)

    concatCSVs(list_csv_file,feature_csv_file, new_csv_file)
    return

def main():

    new_csv_file = data_DIR+"/glFeatures_withid.csv"
    #generateCompleFeatureCSV(new_csv_file)

    df1 = pd.read_csv(new_csv_file)

    if df1.IsPval


    #normalized = normalizeFeatures(glfFeatures)
    #zscores = zscore(glfFeatures)

    #selectFeatures_MRMR(zscores)  # save zscores into a csv file

    #plotFeatureVector(glfFeatures, "normalized global features")
    #plotScoreMatrix(glfFeatures, "Similarity Matrix based on global features")
    #plotScoreMatrix(gmiFeatures, "Similarity Matrix based on GMI features")

    return


if __name__ == "__main__":
      main()
