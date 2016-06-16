### after running the query to obtain the list of experiments into a csv file, run this script to download all swc files accordingly
#python  ../utilities/downloadSWC_filter.py  -i ~/work/data/lims2/june_25_filtered_all_condition.csv -o ~/work/data/lims2/nr_june_25_filter
python  ../utilities/downloadSWC_filter.py  -i ~/work/data/lims2/0729_filtered_ephys_qc.csv	 -o ~/work/data/lims2/0729_filtered_ephys_qc
python  ../utilities/downloadSWC_filter.py  -i ~/work/data/lims2/0903_filtered_ephys_qc.csv	 -o ~/work/data/lims2/0903_filtered_ephys_qc
python ../utilities/downloadSWC_filter.py  -i ~/work/data/lims2/0923_pw_aligned/0923_pw_aligned.csv	 -o ~/work/data/lims2/0923_pw_aligned/original