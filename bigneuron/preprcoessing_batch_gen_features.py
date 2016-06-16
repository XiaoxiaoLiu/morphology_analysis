import os
from os import sys, path
import platform


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)


import utilities.morph_nfb_2_csv as nfb

import blast_neuron.blast_neuron_comp as bn

import glob





def main():
    ###############################################################################
    preprocessing =0
    janelia =0
    taiwan=1
    if taiwan:
        data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/taiwan"
        original_dir = data_DIR + "/consensus_0330_anisosmooth"
        db_tags_csv_file = data_DIR + '/taiwan_smooth_features_with_tags.csv'
    if janelia:
        data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1"
        original_dir = data_DIR + "/consensus_0330_anisosmooth"
        db_tags_csv_file = data_DIR + '/j1_smooth_features_with_tags.csv'
    ###############################################################################

    print original_dir
    preprocessed_dir = data_DIR + "/preprocessed_consensus_smooth"
    if not os.path.exists(preprocessed_dir):
        os.system("mkdir -p  " + preprocessed_dir)

    if preprocessing==1:
        #preprocssing alignment
        count=0
        qsub_folder= "/data/mat/xiaoxiaol/work/qsub"
        os.system("rm "+qsub_folder+"/*.qsub")
        os.system("rm "+qsub_folder+"/*.o*")
        os.system("rm "+qsub_folder+"/jobs.txt")
        for input_swc_path in glob.glob(original_dir + "/*.eswc"):

            swc_fn = input_swc_path.split('/')[-1]

            preprocessed_swc_fn = preprocessed_dir+'/' + swc_fn
            if not os.path.exists(preprocessed_swc_fn):
               bn.pre_processing(input_swc_path, preprocessed_swc_fn,1,qsub_folder,count)
               count=count+1
        exit()  #run jobs on pstar

    print "done"
    preprocessed_ANO = preprocessed_dir + "/preprocessed.ano"
    bn.genLinkerFile(preprocessed_dir, preprocessed_ANO)

    ##batch computing  generate features
    feature_file = preprocessed_dir+'/features.nfb'
    bn.batch_compute(preprocessed_ANO, feature_file)

    ###  convert feature file into csv file
    nfb.generateALLFeatureCSV(feature_file, db_tags_csv_file)
    return


if __name__ == "__main__":
      main()


