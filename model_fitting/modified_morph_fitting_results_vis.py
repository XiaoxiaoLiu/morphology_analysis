#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
import seaborn as sns
import platform


param_names =['PID','A1','A2','Ri','Cm','Rm','error']
sample_size = 55
r_values = [0.5,0.7,0.9,1.0,1.1,1.3,1.5]
p_values = [0.0,0.05,0.1,0.15,0.2,0.25,0.3]




if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

data_DIR = WORK_PATH + "/data/lims2/0923_pw_aligned"


def  my_box_plot(d,output_fig_file_name, XLABEL = "X", YLABEL = "Y"):

        #sns.boxplot(data=d)
        sns.set_context("poster")
        ax=sns.tsplot(data=d.values, time = d.columns,estimator=np.median)
        #plt.xticks(range(len(d.columns)), d.columns)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)

        plt.savefig(output_fig_file_name, bbox_inches="tight")
        plt.close()
        print "saved figure:"+ output_fig_file_name


def  plotResults(data,result_path, column_names,  XLABEL = "X",  YLABEL = "Y",):
    for i in [3,4,5,6]:
#3: Ri  49.2903334908
#4: Cm  2.3995967695
#5: Rm  2782.21957025
#6: Final error  0.000228200455599
        d = pd.DataFrame(data[:,i,:], columns = column_names)
        #d.index.name=param_names[i]
        ylabel = '_'.join(YLABEL.split(" "))+"_"
        output_fig_file_name = result_path +"/"+ ylabel+param_names[i]+"_plot.png"
        my_box_plot(d, output_fig_file_name, XLABEL, YLABEL + ": "+param_names[i])



def  visResults (m, result_dir, varying_para_values = r_values, xlabel= 'Radius Scale Factor', basenum = 0):

        # plot the raw results
        plotResults(m,result_dir, varying_para_values,xlabel, 'Fitted Parameter Value')

        # plot the normalized results
        norm_m = np.zeros(m.shape)
        for i in range(len(varying_para_values)):
            norm_m[:,:,i] = np.divide(m[:,:,i], m[:,:,basenum]) # r=1.0 is the original model fitting results
        plotResults(norm_m,result_dir,varying_para_values,xlabel,'Normalized Fitted Parameter')

        for i in range(len(varying_para_values)):
            norm_m[:,:,i] = np.divide(m[:,:,i]- m[:,:,basenum], m[:,:,basenum]) # r=1.0 is the original model fitting results
        plotResults(norm_m,result_dir,varying_para_values,xlabel,'Fitted Parameter Difference Percentage')



        CmTotal = np.zeros((sample_size, len(varying_para_values)))
        for i in range(len(varying_para_values)):
            CmTotal[:,i] = m[:,4,i] * (m[:,1,i]+ m[:,2,i])  # Cm * (A1+A2) p=0.0 is the original model fitting results
        df_CmTotal = pd.DataFrame(CmTotal, columns = varying_para_values)
        my_box_plot(df_CmTotal, result_dir+ "/total_capacitance.png", xlabel,'Total Capacitance')



        RmUnit = np.zeros((sample_size, len(varying_para_values)))
        for i in range(len(varying_para_values)):
            RmUnit[:,i] = m[:,5,i] / (m[:,1,i]+ m[:,2,i])  # Rm / (A1+A2)
        df_RmUnit = pd.DataFrame(RmUnit, columns = varying_para_values)
        my_box_plot(df_RmUnit, result_dir+ "/total_Rm.png", xlabel,'Unit Membrane Resistance')




# scale up and down radius of neurons to see how it affects passive parameters of the NEURON model fitting
if __name__ == "__main__":
    default_input_para_file = WORK_PATH +'/data/lims2/usable_ephys_para.txt'

    parser = argparse.ArgumentParser(description='passive model fitting batch script ')

    parser.add_argument('-i','-input_para_file',dest="input_para_file",default = default_input_para_file)
    parser.add_argument('-o','-result_dir', dest="result_dir",default = None)
    args = parser.parse_args()

    input_para_file = args.input_para_file
    result_dir = args.result_dir

    if 0:
        result_dir =  WORK_PATH +'/data/lims2/modified_neurons'
        #results are stored as:  r0.5_x1.0_y1.0_z1.0_p0/passive_paras_results.csv
        m = np.zeros( (sample_size,len(param_names),len(r_values)) )
        i = 0
        for radius_scale in r_values:
            filename = result_dir+ '/ephys_fit_result_'+'r'+str(radius_scale)+'_x1.0_y1.0_z1.0_p0/passive_paras_results.csv'
            df = pd.read_csv(filename, sep=',',header=None)
            m[:,:,i] = df.values
            i = i+1

        visResults (m, result_dir, r_values, 'Radius Scale Factor',3)




    if 0:
        result_dir =  WORK_PATH +'/data/lims2/modified_neurons_p'
        m = np.zeros( (sample_size,len(param_names),len(p_values)) )
        i = 0
        for prune_scale in p_values:
            filename = result_dir+ '/ephys_fit_result_'+'r1.0_x1.0_y1.0_z1.0_p'+str(prune_scale)+'/passive_paras_results.csv'
            print filename
            df = pd.read_csv(filename, sep=',',header=None)
            m[:,:,i] = df.values
            i = i+1
        visResults (m, result_dir, p_values, 'Prune Length',0)




    if 1:   # heatmap plot of all values in pairs
        df_pam= pd.read_csv(WORK_PATH +'/data/morph_alt_results/parameter_fit.csv',index_col=0)
        df_fea= pd.read_csv(WORK_PATH +'/data/morph_alt_results/feature_fit.csv',index_col=0)
        df_pam = df_pam.astype(float)
        df_fea = df_fea.astype(float)
        df_pam_diff = pd.DataFrame()
        for col in df_pam.columns[1:5]:
           df_pam_diff[col] =  (df_pam[col]-df_pam.Org0)/df_pam.Org0


        sns.heatmap(df_pam_diff)
        pl.yticks(rotation=90)
        pl.show()
        print df_pam_diff



        df_fea_diff = pd.DataFrame()
        for col in df_fea.columns[1:7]:
           df_fea_diff[col] =  (df_fea[col]-df_fea.Org0)/df_fea.Org0

        sns.heatmap(df_fea_diff)
        pl.yticks(rotation=90)
        print df_fea_diff
        pl.show()


