__author__ = 'xiaoxiaol'


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import numpy as np
from scipy import stats



def read_swc(infile):
    ###n,type,x,y,z,radius,parent
     swc_df = pd.read_csv(infile, sep=" ",comment="#", header =None, names =['n','type','x','y','z','radius','parent'] )
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
    data_DIR ="/data/mat/xiaoxiaol/data/BBP/BBP_data_uncurated_SWC"
    SWC_DIR = data_DIR+'/sorted'
    df_features_with_tags = pd.read_csv(data_DIR+'/bbp_inhibitory_cloud_xy_features.csv')
    metric = 'm-type'
    dfg = df_features_with_tags.groupby(metric)


    for cre in np.unique(df_features_with_tags[metric]):
        print cre
        anofile = data_DIR +'/'+cre+'.ano'
        f=open(anofile, 'w')
        df_cre = dfg.get_group(cre)
        coords_x=[]
        coords_y=[]
        for swc_file in df_cre['swc_file_name']:
             swc_file = SWC_DIR + '/'+swc_file
             line = "SWCFILE="+swc_file+"\n"
             f.write(line)
             print swc_file
             df_swc = read_swc(swc_file)

             #axon
             df_swc=df_swc[df_swc['type'] == 2]
             #flip y
             coords_y.extend(np.round(-df_swc.y))
             coords_x.extend(np.round(df_swc.x))
        print "total number of node:",len(coords_x), " and ", len(coords_y)
        f.close()

        plt.figure()
        #lm= sn.distplot(coords,vertical=True, hist=True, kde=True)
        g = sn.jointplot(np.array(coords_x), np.array(coords_y), kind="kde", size=7, space=0, xlim=(-1500, 1500), ylim=(-1500,1500))

        plt.ylabel('-y')
        plt.xlabel('x')
        plt.title(cre+"("+str(len(df_cre))+")")
        #plt.show()
        plt.savefig(data_DIR+'/density_2d_'+cre+'.png')
	plt.close()


if __name__ == "__main__":
        main()

