#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


param_names =['PID','A1','A2','Ri','Cm','Rm','error']
sample_size = 55
r_values = [0.5,0.7,0.9,1.0,1.1,1.3,1.5]
p_values = np.arange(0,0.35,0.05)


def  plotResults(data,result_path, column_names,  XLABEL = "X",  YLABEL = "Y",):
    for i in [3,4,5,6]:
#3: Ri  49.2903334908
#4: Cm  2.3995967695
#5: Rm  2782.21957025
#6: Final error  0.000228200455599
        d = pd.DataFrame(data[:,i,:], columns = column_names)
        fig_file_name = result_path +"/" +param_names[i]+"_plot.png"
        my_box_plot(d, result_path, column_names, XLABEL, YLABEL, fig_file_name)


def  my_box_plot(data,result_path, column_names,  XLABEL = "X", YLABEL = "Y", output_fig_fn):
        d = pd.DataFrame(data, columns = column_names)
        f = sns.boxplot(data=d)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)

        fig_file_name = result_path +"/" + output_fig_fn
        plt.savefig(fig_file_name, bbox_inches="tight")
        print "saved figure:"+ fig_file_name



# scale up and down radius of neurons to see how it affects passive parameters of the NEURON model fitting
if __name__ == "__main__":
    default_input_para_file = '/local1/xiaoxiaol/work/data/lims2/usable_ephys_para.txt'

    parser = argparse.ArgumentParser(description='passive model fitting batch script ')

    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default = None)
    args = parser.parse_args()

    input_para_file = args.input_para_file
    result_dir = args.result_dir

    if True:
        result_dir = '/home/xiaoxiaol/work/data/lims2/modified_neurons'
        #results are stored as:  r0.5_x1.0_y1.0_z1.0_p0/passive_paras_results.csv
        m = np.zeros( (sample_size,len(param_names),len(r_values)) )
        i = 0
        for radius_scale in r_values:
            filename = result_dir+ '/ephys_fit_result_'+'r'+str(radius_scale)+'_x1.0_y1.0_z1.0_p0/passive_paras_results.csv'
            df = pd.read_csv(filename, sep=',',header=None)
            m[:,:,i] = df.values
            i = i+1


        # plot the raw results
        plotResults(m,result_dir, r_values,'Radius Scale', 'Fitted Parameter Value')


        norm_m = np.zeros(m.shape)
        for i in range(len(r_values)):
            norm_m[:,:,i] = np.divide(m[:,:,i], m[:,:,3]) # r=1.0 is the original model fitting results
        #plotResults(norm_m,result_dir,'Radius Scale','Fitted Parameter Factor')

        CmTotal = np.zeros((sample_size, len(r_values)))
        for i in range(len(r_values)):
            CmTotal[:,i] = m[:,4,i] * (m[:,1,i]+ m[:,2,i])  # Cm * (A1+A2) p=0.0 is the original model fitting results
        #plotResults3(CmTotal,result_dir,'Radius Scale','Total Capacitance')

        RmTotal = np.zeros((sample_size, len(r_values)))
        for i in range(len(p_values)):
            RmTotal[:,i] = m[:,5,i] / (m[:,1,i]+ m[:,2,i])  # Rm / (A1+A2) 
        #plotResults3(RmTotal,result_dir,'Prune Threshold','Total Membrane Resistance')

    if True:
        result_dir = '/home/xiaoxiaol/work/data/lims2/modified_neurons_p'
        m = np.zeros( (sample_size,len(param_names),len(p_values)) )
        i = 0
        for prune_scale in p_values:
            filename = result_dir+ '/ephys_fit_result_'+'r1.0_x1.0_y1.0_z1.0_p'+str(prune_scale)+'/passive_paras_results.csv'
            df = pd.read_csv(filename, sep=',',header=None)
            m[:,:,i] = df.values
            i = i+1


        # plot the raw results
        plotResults2(m,result_dir,'Prune Threshold', 'Fitted Parameter Value')


        norm_m = np.zeros(m.shape)
        for i in range(len(p_values)):
            norm_m[:,:,i] = np.divide(m[:,:,i], m[:,:,0])  # p=0.0 is the original model fitting results
        plotResults2(norm_m,result_dir,'Prune Threshold','Fitted Parameter Factor')

        CmTotal = np.zeros((sample_size, len(p_values)))
        for i in range(len(p_values)):
            CmTotal[:,i] = m[:,4,i] * (m[:,1,i]+ m[:,2,i])  # Cm * (A1+A2) 

        plotResults3(CmTotal,result_dir,'Prune Threshold','Total Capacitance')

        RmTotal = np.zeros((sample_size, len(p_values)))
        for i in range(len(p_values)):
            RmTotal[:,i] = m[:,5,i] / (m[:,1,i]+ m[:,2,i])  # Rm / (A1+A2) 
        plotResults3(CmTotal,result_dir,'Prune Threshold','Total Membrane Resistance')
