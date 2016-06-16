# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:40:51 2015

@author: xiaoxiaol
"""

#Standardized swc files (www.neuromorpho.org) -
#0 - undefined
#1 - soma
#2 - axon
#3 - (basal) dendrite
#4 - apical dendrite
#5+ - custom


import numpy
import matplotlib.pylab as pl
import pandas as pd
import seaborn as sns





data_DIR = "/data/mat/xiaoxiaol/data/IVSCC_testing/swc_pixel"
table = pd.DataFrame()




#df_ref = pd.read_csv('/data/mat/xiaoxiaol/data/IVSCC_testing/list.csv')
df_ref = pd.read_csv('/data/mat/xiaoxiaol/data/lims2/0729_filtered_ephys_qc.csv')

def  showMetric(metric_name, m_table, title):
    pl.figure()
    sns.boxplot(x='specimen_id',y=metric_name,data = m_table)
    pl.xticks(rotation=90)
    pl.title(title)
    pl.tight_layout()
    pl.savefig(data_DIR+'/'+metric_name+'.pdf')
    pl.show()


def  showBarPlotMetric(metric_name, m_table, title, output_fn =None):
    pl.figure()
    #sns.barplot(x='specimen_id', y=metric_name,hue='segment_type', data=m_table)
    sns.barplot(x='specimen_id', y=metric_name, data=m_table)

    pl.xticks(rotation=90)
    #pl.xticks()
    pl.title(title)
    pl.tight_layout()
    if ( output_fn is None):
        output_fn = data_DIR+'/'+metric_name+title+'.barplot.pdf'
    pl.savefig(data_DIR+'/'+output_fn)
    pl.show()




def  showMeanStdBarMetric(metric_name, m_table,ids, title):
    pl.figure()
    m_mean = numpy.zeros(ids.size)
    m_std = numpy.zeros(ids.size)
    for i in range(ids.size):
        numbers = m_table[metric_name][m_table['specimen_id']==ids[i]]
        m_mean[i]=  numbers.mean()
        m_std[i] = numbers.std()
    #sort by mean
    order = m_mean.argsort()
    pl.errorbar(range(ids.size),m_mean[order[::-1]],yerr=m_std[order[::-1]], fmt='o')
    pl.xticks(range(ids.size), ids[order[::-1]], rotation='vertical')
    pl.title(title)
    pl.savefig(data_DIR+'/'+metric_name+'.mean_std_bar.pdf')
    pl.show()


####################
with open(data_DIR+'/list.txt','r') as f:
    i = 0
    for line in f:
            csv_path =  data_DIR + "/"+ line.split("\n")[0]
            df = pd.read_csv(csv_path)
            df['csv_path'] = csv_path.split('/')[-1]
            swc_file = csv_path.split('/')[-1][0:-4]
            specimen_id = df_ref[  df_ref['orca_path'].split('/')[-1] == swc_file]['specimen_id']
            #if  numpy.count_nonzero(  df['segment_type'] == 2 )  < 2  :
            #        df = df[df['segment_type'] != 2]
            df['specimen_id'] =  specimen_id.values[0]
            frames= [table,df]
            table = pd.concat(frames,ignore_index=True)
            i = i+1;


# merge apical dendrite  (4) to dendrite (3)
table['segment_type'] = table['segment_type'].replace(4,3)


table = table[table['segment_type']> 1]
before = table.shape[0]

#clean up data
table = table[table['tubularity']<50]
table = table[table['cnr']<200]
table = table[table['snr']<200]
after = table.shape[0]
print "deleted "+ str(before-after) +'rows!'



############## append confocal studies
INCLUDE_CONFOCAL = 1
if INCLUDE_CONFOCAL:
    other_table = pd.DataFrame()
    with open(data_DIR+'/../compare_list.csv','r') as f:
       i = 0
       for line in f:
            csv_path =  data_DIR + "/../"+ line.split("\n")[0]
            df = pd.read_csv(csv_path)
            df['csv_path'] = csv_path.split('/')[-1]
            swc_file = csv_path.split('/')[-1][0:-4]

            #if  numpy.count_nonzero(  df['segment_type'] == 2 )  < 2  :
             #       df = df[df['segment_type'] != 2]

            df['specimen_id'] =  i
            frames= [other_table,df]
            other_table = pd.concat(frames,ignore_index=True)
            i = i+1;


############## append confocal studies

table_filtered = table
# remove soma, and axon
#table_filtered = table[table['segment_type'] == 2]  #axon
#table_filtered = table[table['segment_type'] == 3]  #dendrite


# rank swc by median

ids = numpy.unique(table_filtered['specimen_id'])
s = numpy.zeros(ids.size)
# size of s and ids should be the same

for iii in range(ids.size):
      cnrs = table_filtered['cnr'][table_filtered['specimen_id']==ids[iii]]
      #cnrs.sort()
      #chopoff = int(cnrs.size*0.1)  # to be more robust towards wrongly segmented neurites
      #s[iii]= cnrs[chopoff:-chopoff].mean()
      s[iii] = cnrs.mean()
rank = s.argsort().argsort()

id_checking_map = dict(zip(rank,ids))


#sort by cnr mean

my_map = dict(zip(ids,rank))
#v = table_filtered['specimen_id'].values
vv = numpy.zeros(table_filtered.shape[0])
for i in range(table_filtered.shape[0]):
       vv[i] = my_map[table_filtered['specimen_id'].values[i]]

table_filtered.insert(1,'cnr_rank',vv)

sorted_table = table_filtered.sort(['cnr_rank','segment_id'], ascending=[0,1])


frames= [other_table,sorted_table]
all_table = pd.concat(frames,ignore_index=True)

PLOT = 1
if PLOT:
    #showBarPlotMetric('snr',sorted_table,'AXON-Signal to Noise Ratio')
    #showBarPlotMetric('cnr',sorted_table,'AXON-Constrat to Noise Ratio')
    showBarPlotMetric('cnr',all_table,'Constrat to Noise Ratio','all_with_comapre_Barplot.pdf')

   #showBarPlotMetric('dynamic_range',sorted_table, 'AXON-Dynamic Range')
    #showBarPlotMetric('tubularity',sorted_table, 'AXON-Tubularity')

    #output
    #showMetric('snr',sorted_table,'AXON-Signal to Noise Ratio')
    ##showMetric('cnr',sorted_table,'AXON-Constrat to Noise Ratio')
    #showMetric('dynamic_range',sorted_table, 'AXON-Dynamic Range')
    #showMetric('tubularity',sorted_table, 'AXON-Tubularity')


    #showMeanStdBarMetric('snr',sorted_table,ids,'AXON-Signal to Noise Ratio')
    #showMeanStdBarMetric('cnr',sorted_table,ids,'AXON-Constrat to Noise Ratio')
    #showMeanStdBarMetric('dynamic_range',sorted_table,ids, 'AXON-Dynamic Range')
    #showMeanStdBarMetric('tubularity',sorted_table, ids,'AXON-Tubularity')










#pl.figure()
#ax1=sns.boxplot(x='specimen_id',y='dynamic_range', hue='segment_type',data = table_filtered)
#pl.title('Dynamic Range')

