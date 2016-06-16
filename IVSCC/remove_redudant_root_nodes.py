__author__ = 'xiaoxiaoliu'

import pandas as pd

def read_swc(infile):
    ###n,type,x,y,z,radius,parent
     swc_df = pd.read_csv(infile, sep=" ",skiprows = 3, header =None, names =['n','type','x','y','z','radius','parent'] )
     return swc_df



def write_swc(df_out,outfile):
     ###n,type,x,y,z,radius,parent
     df_out.to_csv(outfile, index=False, header= False, sep=" ")


def extract_swc_by_id(type_id, swc_file, output_swc_file = None):
      if (output_swc_file == None):
          output_swc_file = swc_file[0:-4] +'_'+ str(type_id) +'.swc'
      swc_df = read_swc(swc_file)
      extracted = swc_df[ swc_df['type'] == type_id ]  # keep the soma
      soma = swc_df[ swc_df['type'] == 1 ]

      frames = [soma, extracted]
      result = pd.concat(frames)

      write_swc(result,output_swc_file)



def remove_redundant_roots(swc_file, output_swc_file):
    swc_df = read_swc(swc_file)
    root_row = swc_df.iloc[0]
    print root_row
    out_swc_df = swc_df[swc_df['parent'] != -1]

    write_swc(out_swc_df, output_swc_file)



data_DIR = '/data/mat/xiaoxiaol/data/BBP/raw_converted_swcs/swc_files'
output_DIR = '/data/mat/xiaoxiaol/data/BBP/raw_converted_swcs/cleaned'
df_list = pd.read_csv('/data/mat/xiaoxiaol/data/BBP/raw_converted_swcs/swc_files/list.csv')
for i in range(df_list.shape[0]):
     swc_file = data_DIR +'/'+ df_list.iloc[i]['swc_file_name']
     output_swc_file  = output_DIR +'/'+ df_list.iloc[i]['swc_file_name']
     remove_redundant_roots(swc_file, output_swc_file)
     exit()
