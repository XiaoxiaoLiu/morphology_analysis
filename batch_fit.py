#!/usr/bin/env python

import os
import argparse


def run_one_fit(specimen_id, down, up, result_dir):
	cmd ='python neuron_passive_fit.py ' + specimen_id + ' ' + down + ' ' + up + ' -o '+ result_dir
        cmd += ' >' +result_dir+ '/'+specimen_id+'.log'
        os.system(cmd)
        return
 
def run_check(specimen_id, result_dir):
	cmd ='python average_cap_check.py '+specimen_id +' --noshow 1'+' -o '+ result_dir
        os.system(cmd)
        return

if __name__ == "__main__":
    default_input_para_file = '/local1/xiaoxiaol/work/data/lims2/usable_ephys_para.txt'
    default_result_dir = '/local1/xiaoxiaol/work/data/lims2/ephys_fit_result2'

    parser = argparse.ArgumentParser(description='passive model fitting batch script ')

    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default=default_result_dir)
    args = parser.parse_args()
     
    input_para_file = args.input_para_file
    result_dir = args.result_dir
    
    #read from the listn
    with open(input_para_file,'r') as f:
         for line in f:
              (specimen_id, down, up) = line.split()
              print "running: ",specimen_id,' ', down,' ', up
              run_check(specimen_id,result_dir)
	      run_one_fit(specimen_id, down,up, result_dir)

