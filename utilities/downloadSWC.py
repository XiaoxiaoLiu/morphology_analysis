#!/usr/bin/env python

import lims_orca_utils as lu
import argparse
import pandas as pd
import os

def updateLocalDB():
   # grab new data from lims
   return

def copyBySpecimenId(df): 
       ids = df['specimen__id']
       for i in range(len(ids)):
            swc_filename, swc_path = lu.get_swc_from_lims(str(ids[i]))
            cmd = 'cp '+swc_path +' '+ save_data_dir+'/'
            print cmd
            os.system(cmd)

def copyBySpecimenName(df): 
       ids = df['specimen']
       for i in range(len(ids)):
            swc_filename, swc_path = lu.get_swc_from_lims_by_specimen_name(str(ids[i]))
            cmd = 'cp '+swc_path +' '+ save_data_dir+'/'
            print cmd
            os.system(cmd)

def copyByNrid(df): 
       ids = df['nrid']
       for i in range(len(ids)):
            swc_filename, swc_path = lu.get_swc_from_lims_by_nrid(ids[i])
            cmd = 'cp '+swc_path +' '+ save_data_dir+'/'
            print cmd
            os.system(cmd)

if __name__ == "__main__":
        #default_input_csv_file = '/local1/xiaoxiaol/work/data/lims2/pvalb.csv'
	#default_save_data_dir = '/local1/xiaoxiaol/work/data/lims2/pvalb'
        default_input_csv_file = '/local1/xiaoxiaol/work/data/lims2/custom_report-IVSCC_classification-April_2015.csv'
	default_save_data_dir = '/local1/xiaoxiaol/work/data/lims2/neuron_recon_latest'

	parser = argparse.ArgumentParser(description='download swc files from lims by specimen id in the csv file ')
	parser.add_argument('-i','-input_csv_file',dest="input_csv_file",default = default_input_csv_file)
	parser.add_argument('-o','-save_data_dir', dest="save_data_dir",default = default_save_data_dir)
	args = parser.parse_args()

	input_csv_file = args.input_csv_file
	save_data_dir = args.save_data_dir


        os.system('mkdir '+save_data_dir)
        
        df = pd.read_csv(input_csv_file)

        #copyBySpecimenId(df)
        copyBySpecimenName(df)

        #copyByNrid(df)
