#!/usr/bin/env python

import os
import argparse
import numpy as np



def collect_fit_results(specimen_id, result_dir):
        logfile= result_dir+ '/'+specimen_id+'.log'
        Ri=0
        Cm=0
        Ci=0
        err=0
        print logfile
        try:
             f = open(logfile,'r')
             lines=f.readlines()
             if lines:   
                if lines[-4].split()[0] == 'Ri':
	    	    t,Ri=lines[-4].split() 
		    t,Cm=lines[-3].split() 
		    t,Ci=lines[-2].split() 
		    t,s,err=lines[-1].split() 
        except :
            print "not found"

        return (float(Ri), float(Cm),float(Ci),float(err))

if __name__ == "__main__":
    default_result_dir = '/local1/xiaoxiaol/work/data/lims2/ephys_fit_result4' 
    default_input_para_file = '/local1/xiaoxiaol/work/data/lims2/usable_ephys_para.txt'
    parser = argparse.ArgumentParser(description='passive model fitting batch script ')

    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default=default_result_dir)
    args = parser.parse_args()

    input_para_file = args.input_para_file
    result_dir = args.result_dir
    
    #paras_all=['specimen id', 'Ri','Cm','Rm','Err' ]
    paras_all=[]
    with open(input_para_file,'r') as f:
         for line in f:
              specimen_id, down, up = line.split()
              ri, cm, ci,err = collect_fit_results(specimen_id, result_dir)
              paras_all.append([int(specimen_id),ri, cm,ci,err])
        
    paras_matrix = np.array(paras_all)     
    np.set_printoptions(precision=5) 
    np.set_printoptions(suppress=True)
    print (paras_matrix)
    np.savetxt(result_dir +"/passive_paras_results.csv",paras_matrix,delimiter=",")
    print "    mean   std"
    print "Ri:",np.mean(paras_matrix[:,1]),np.std(paras_matrix[:,0])
    print "Cm: ", np.mean(paras_matrix[:,2]),np.std(paras_matrix[:,1])
    print "Rm:",np.mean(paras_matrix[:,3]),np.std(paras_matrix[:,2])
    print "Err:",np.mean(paras_matrix[:,4]),np.std(paras_matrix[:,3])

