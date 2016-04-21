__author__ = 'xiaoxiaol'


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
from scipy import stats



def read_swc(infile):
    ###n,type,x,y,z,radius,parent
     swc_df = pd.read_csv(infile, sep=" ", header =None, comment="#",names =['n','type','x','y','z','radius','parent'] )
     return swc_df


def write_swc(df_out,outfile):
     ###n,type,x,y,z,radius,parent
     df_out.to_csv(outfile, index=False, header= False, sep=" ")


def density2d_plot(swc_file):
    df_swc = read_swc(swc_file)
    sn.distplot(df_swc.z,vertical=True)
    return



def main():
    # read swc files
    data_DIR ="/data/mat/xiaoxiaol/data/lims2/ivscc_0411"
    SWC_DIR = data_DIR+'/SWC'
    df_features_with_tags = pd.read_csv(data_DIR+'/ivscc_0411_features_with_meta.csv')

    dfg = df_features_with_tags.groupby('cre_line')

    
    for cre in np.unique(df_features_with_tags['cre_line']):
        print cre
        anofile = data_DIR +'/'+cre+'.ano'
        f=open(anofile, 'w')
        df_cre = dfg.get_group(cre)
        z_coords=[]
        for swc_file in df_cre['swc_file_name']:
             swc_file = SWC_DIR + '/'+swc_file[:-4]+"_pia.swc"
             line = "SWCFILE="+swc_file+"\n"
             f.write(line)
             print swc_file
             df_swc = read_swc(swc_file)
             z_coords.extend(df_swc.z)
        print "total number of node:",len(z_coords)
        f.close()

        plt.figure()
        lm= sn.distplot(z_coords,vertical=True, hist=True, kde=False)
        lm.set(ylim=(-1000, 5000))

        plt.ylabel('depth')
        plt.xlabel('counts')
        plt.title(cre+"("+str(len(df_cre))+")")
        #plt.show()
        plt.savefig(data_DIR+'/counts_z_'+cre+'.png')
	plt.close()


if __name__ == "__main__":
        main()
