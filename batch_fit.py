#!/usr/bin/env python

import os
import argparse


def run_one_fit(specimen_id, down, up, result_dir, local_path = None):
	cmd ='python neuron_passive_fit.py ' + specimen_id + ' ' + down + ' ' + up + ' -o '+ result_dir
        if local_path:
           cmd += ' -i '+ local_path
        cmd += ' >' +result_dir+ '/'+specimen_id+'.log'
        os.system(cmd)
        return
 
def run_check(specimen_id, result_dir):
	cmd ='python average_cap_check.py '+specimen_id +' --noshow 1'+ ' -o '+ result_dir
        os.system(cmd)
        return

if __name__ == "__main__":
    default_input_para_file = '/local1/xiaoxiaol/work/data/lims2/usable_ephys_para.txt'
    default_result_dir = '/local1/xiaoxiaol/work/data/lims2/ephys_fit_result_modified'
    default_local_path = '/home/xiaoxiaol/work/data/lims2/modified_neurons/r0.5_x1.0_y1.0_z1.0_p0'

    parser = argparse.ArgumentParser(description='passive model fitting batch script ')

    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default = default_result_dir)
    parser.add_argument('-ii','-local_path', dest="local_path", default = None) # modified data
    args = parser.parse_args()
     
    input_para_file = args.input_para_file
    result_dir = args.result_dir
    local_path = args.local_path
    
    #read from the listn
    with open(input_para_file,'r') as f:
         count = 0 
         for line in f:
              count = count +1
              (specimen_id, down, up) = line.split()
              print "running: ",count, specimen_id,' ', down,' ', up
              run_check(specimen_id,result_dir)
              if local_path :
	          run_one_fit(specimen_id, down,up, result_dir, local_path)
              else:
	          run_one_fit(specimen_id, down,up, result_dir)

