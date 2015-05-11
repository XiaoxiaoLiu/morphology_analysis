#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def  plotResults(data,result_path, XLABEL = "X", YLABEL = "Y",):

    for i in [1,2,3,4]:
        d = pd.DataFrame(data[:,i,:], columns = r_values)
        f = d.plot(kind='box',positions=[1, 3, 5, 6, 7, 9,11])
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)

        # Remove grid lines (dotted lines inside plot)
        f.grid(False)
        plt.title(param_names[i])

        fig_file_name = result_path +"/" + XLABEL + "_"+  YLABEL+ "_"+ param_names[i]+"_plot.png"
        plt.savefig(fig_file_name, bbox_inches="tight")
        print "saved figure:"+ fig_file_name

    plt.show()


def  plotResults2(data,result_path, XLABEL = "X", YLABEL = "Y",):

    for i in [1,2,3,4]:
        d = pd.DataFrame(data[:,i,:], columns = p_values)
        f = d.plot(kind='box')
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)

        # Remove grid lines (dotted lines inside plot)
        f.grid(False)
        plt.title(param_names[i])

        fig_file_name = result_path +"/" + XLABEL + "_"+  YLABEL+ "_"+ param_names[i]+"_plot.png"
        plt.savefig(fig_file_name, bbox_inches="tight")
        print "saved figure:"+ fig_file_name

    plt.show()



param_names =['PID','Ri','Cm','Rm','error']
sample_size = 55
r_values = [0.5,0.7,0.9,1.0,1.1,1.3,1.5]
p_values = [0.0,0.05,0.15,0.2,0.25,0.3] #np.arange(0,0.35,0.05)

# scale up and down radius of neurons to see how it affects passive parameters of the NEURON model fitting
if __name__ == "__main__":
    default_input_para_file = '/local1/xiaoxiaol/work/data/lims2/usable_ephys_para.txt'
    default_result_path = '/home/xiaoxiaol/work/data/lims2/modified_neurons_p'

    parser = argparse.ArgumentParser(description='passive model fitting batch script ')

    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default = default_result_path)
    args = parser.parse_args()

    input_para_file = args.input_para_file
    result_dir = args.result_dir

    if False:
        #results are stored as:  r0.5_x1.0_y1.0_z1.0_p0/passive_paras_results.csv
        m = np.zeros( (sample_size,len(param_names),len(r_values)) )
        i = 0
        for radius_scale in r_values:
            filename = result_dir+ '/ephys_fit_result_'+'r'+str(radius_scale)+'_x1.0_y1.0_z1.0_p0/passive_paras_results.csv'
            df = pd.read_csv(filename, sep=',',header=None)
            m[:,:,i] = df.values
            i = i+1


        # plot the raw results
        plotResults(m,default_result_path,'Radius Scale', 'Fitted Parameter Value')


        norm_m = np.zeros(m.shape)
        for i in range(len(r_values)):
            norm_m[:,:,i] = np.divide(m[:,:,i], m[:,:,3]) # r=1.0 is the original model fitting results
        plotResults(norm_m,default_result_path,'Radius Scale','Fitted Parameter Factor')

    if True:
        m = np.zeros( (sample_size,len(param_names),len(p_values)) )
        i = 0
        for prune_scale in p_values:
            filename = result_dir+ '/ephys_fit_result_'+'r1.0_x1.0_y1.0_z1.0_p'+str(prune_scale)+'/passive_paras_results.csv'
            df = pd.read_csv(filename, sep=',',header=None)
            m[:,:,i] = df.values
            i = i+1


        # plot the raw results
        print m.shape
        plotResults2(m,default_result_path,'Prune Threshold', 'Fitted Parameter Value')


        norm_m = np.zeros(m.shape)
        for i in range(len(p_values)):
            norm_m[:,:,i] = np.divide(m[:,:,i], m[:,:,0])  # p=0.0 is the original model fitting results
        plotResults2(norm_m,default_result_path,'Prune Threshold','Fitted Parameter Factor')

