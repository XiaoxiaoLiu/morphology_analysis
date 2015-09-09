#!/usr/bin/env python

import os
import argparse
import numpy as np
from utilities.lims_orca_utils import *


VAA3D_EXE ='/local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d'
PRUNNING_EXE = '/home/xiaoxiaol/work/src/blastneuron/bin/prune_short_branch'

def run_radius_modifier(inputSWC_file, outputSWC_file,r_scale, x_scale, y_scale, z_scale, prune_level ):
        cmd =VAA3D_EXE + ' -x neuron_utilities/NeuronModifier -f modify -i ' \
               + inputSWC_file + ' '+ '-o '+ outputSWC_file + \
               ' -p '+ str(r_scale) +' '+ str(x_scale) + ' ' + str(y_scale)  + ' ' + str(z_scale) + ' '+ str(prune_level)
        print cmd
        os.system(cmd)
        return

def run_resort_swc(inputSWC_file, outputSWC_file):
        cmd =VAA3D_EXE + ' -x sort_neuron_swc -f sort_swc -i ' \
               + inputSWC_file + ' '+ '-o '+ outputSWC_file
        print cmd
        os.system(cmd)
        return
  
def run_prune_modifier(inputSWC_file, outputSWC_file, prune_thres):
        cmd =PRUNNING_EXE + ' -i ' \
               + inputSWC_file + ' '+ '-o '+ outputSWC_file + \
               ' -t '+ str(prune_thres)
        print cmd
        os.system(cmd)
        return

if __name__ == "__main__":
    default_input_para_file = '/local1/xiaoxiaol/work/data/lims2/usable_ephys_para.txt'
    default_result_dir = '/local1/xiaoxiaol/work/data/lims2/modified_neurons_p'

    parser = argparse.ArgumentParser(description='passive model fitting batch script ')
    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default = default_result_dir)
    args = parser.parse_args()
     
    input_para_file = args.input_para_file
    result_dir = args.result_dir
    
    r_scale = 1.0
    x_scale = 1.0
    y_scale = 1.0
    z_scale = 1.0 
    prune_level = 0

    if False:  # RADIUS
        for r_scale in np.arange(0.5,1.6,0.1):
            subfolder = 'r'+ str(r_scale) +'_x'+str(x_scale) + '_y'+str(z_scale)+ '_z'+str(z_scale) + '_p'+str(prune_level)
            outputfolder = result_dir +'/'+subfolder
            if not os.path.isdir(outputfolder):
                os.system('mkdir -p '+outputfolder)

            with open(input_para_file,'r') as f:
                for line in f:
                    (specimen_id, down, up) = line.split()
                    print "generate modified neuron with specimen id : ",specimen_id
                    swc_filename, swc_path = get_swc_from_lims(specimen_id)
                    inputSWC_file = swc_path
                    outputSWC_file = outputfolder +'/'+swc_filename
                    run_radius_modifier(inputSWC_file, outputSWC_file, r_scale, x_scale,y_scale, z_scale, prune_level )

    if True:  # BRANCH
        for prune_level in np.arange(0.05,0.35,0.05):
            subfolder = 'r'+ str(r_scale) +'_x'+str(x_scale) + '_y'+str(z_scale)+ '_z'+str(z_scale) + '_p'+str(prune_level)
            outputfolder = result_dir +'/'+subfolder
            if not  os.path.isdir(outputfolder):
                os.system('mkdir -p '+outputfolder)

            with open(input_para_file,'r') as f:
                for line in f:
                    (specimen_id, down, up) = line.split()
                    print "generate modified neuron with specimen id : ",specimen_id
                    print "prune:", prune_level
                    swc_filename, swc_path = get_swc_from_lims(specimen_id)
                    inputSWC_file = swc_path
                    outputSWC_file = outputfolder +'/'+swc_filename
                    run_prune_modifier(inputSWC_file, outputSWC_file, prune_level)
                    run_resort_swc(outputSWC_file, outputSWC_file)

