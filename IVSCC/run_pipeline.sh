# This script stream line the morphology analysis process
# author: Xiaoxiao Liu  Sept, 2015

# 1. run sql qury :
# alignment paramters: pia_alignment_transforms.sql
# swc file database tags ( dendrite_type, layer info) : query_nr_va_ephys_qc.sql

# 2. download all swc files accordingly
#python ../utilities/downloadSWC_filter.py  -i ~/work/data/lims2/0923_pw_aligned/0923_pw_aligned.csv	 -o ~/work/data/lims2/0923_pw_aligned/original


# 3. apply pia white matter alignment
# the following script removes 464326095 (alignment paras are way too off) Pvalb-IRES-Cre_Ai14_IVSCC_-165172.06.02.01_475465899_m.swc
python applyPiaTransformToSWCs.py



remove_axon.py

chang axon_removed to preprocessed

#4. run feature calculation, make sure resample is disabled
python run_feature_calculation.py

features_with_db_tags
#generate three more features in excel  (change height, width,depth
add_features.py

# 5. clustering
# python merge_meta_info.py
# morph_cluster.sh
