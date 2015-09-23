# This script stream line the morphology analysis process
# author: Xiaoxiao Liu  Sept, 2015

# 1. run sql qury :
# alignment paramters: pia_alignment_transforms.sql
# swc file database tags ( dendrite_type, layer info) : query_nr_va_ephys_qc.sql

# 2. download all swc files accordingly
python ../utilities/downloadSWC_filter.py  -i ~/work/data/lims2/0923_pw_aligned/0923_pw_aligned.csv	 -o ~/work/data/lims2/0923_pw_aligned/original


# 3. apply pia white matter alignment
python applyPiaTransformToSWCs.py


#4. run feature calculation
python run_feature_calculation.py


# 5. clustering
python hierachical_cluster_vis.py