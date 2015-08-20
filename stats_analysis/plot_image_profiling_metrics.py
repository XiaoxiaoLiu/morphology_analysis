# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 10:15:53 2015

@author: xiaoxiaoliu
"""


#Standardized swc files (www.neuromorpho.org) -
#0 - undefined
#1 - soma
#2 - axon
#3 - (basal) dendrite
#4 - apical dendrite
#5+ - custom


import matplotlib.pylab as pl
import pandas as pd
import seaborn as sns


#data_DIR = "/data/mat/xiaoxiaol/data/IVSCC_testing/swc_pixel"
data_DIR = "/data/mat/xiaoxiaol/data/IVSCC_testing/swc_enhanced_4thrun_preprocessed"
#df_ref = pd.read_csv('/data/mat/xiaoxiaol/data/lims2/0729_filtered_ephys_qc_edit.csv')


def showMetric(metric_name, m_table, title, output_fn = None):
    pl.figure()
   # siz =  m_table.shape[0]
    #pl.plot(range(siz), m_table[metric_name], color="r", lw=2)
    sns.factorplot(x='specimen_id',y=metric_name, data = m_table, kind="bar")
  #  labels = range(1,siz,30)
    pl.xticks(fontsize=3)
    pl.title(title)
    pl.xticks(rotation=90)
   # pl.xlabel('Image ID 1 ~ ' + str(siz))
    if ( output_fn is None):
        output_fn = data_DIR+'/'+metric_name+title+'.pdf'
    pl.savefig(output_fn)
    pl.show()


###################### plottting function #######################
def  showMeanSTDMetric(metric_name, m_table, title, output_fn = None):
    pl.figure()
    ax = pl.subplot(111)
    #pl.set_context("talk", font_scale=1, rc={"lines.linewidth": 0.2})
    siz = table.shape[0]
  #  ax=pl.errorbar(range(siz),m_table[metric_name],yerr=m_table[metric_name], fmt='o')

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    pl.fill_between(range(siz), m_table[metric_name+'_mean'] - m_table[metric_name+'_std'],
                m_table[metric_name+'_mean'] + m_table[metric_name+'_std'], color="#3F5D7D",alpha=0.8, lw=0.001)
    pl.plot(range(siz), m_table[metric_name+'_mean'], color="r", lw=2)

    labels = range(1,siz,30)

    #g = sns.FacetGrid(m_table)
    #g.map(pl.errorbar, "csv_path", metric_name, yerr= m_table[metric_name+'_std'] , marker="o")

    pl.xticks(labels,labels,fontsize=9)
    pl.title(title)
    pl.xlabel('Image ID 1 ~ ' + str(siz))
    if ( output_fn is None):
        output_fn = data_DIR+'/'+metric_name+title+'.mean_std.pdf'
    pl.savefig(output_fn)
    pl.show()


#####################  MAIN   ############################




df_list = pd.read_csv(data_DIR+'/list.csv')
table = pd.DataFrame()    
for i in range(df_list.shape[0]):
            swc_file =  df_list.swc_file[i]
            csv_path = data_DIR + '/'+swc_file + '.csv'
            df = pd.read_csv(csv_path)
            df = df[df['segment_type'] == -1]  # all-swcnode stats
            
          #  specimen_id = df_ref[  df_ref['swc_file'] == swc_file]['specimen_id']       
           # df['specimen_id'] =  specimen_id.values[0]
            df['specimen_id'] =  i
            df['swc_file'] =swc_file
            frames = [table,df]
            table = pd.concat(frames,ignore_index = True)


print "done reading " + str(i) + " csv files."



#sort by cnr mean
rank = table['cnr'].argsort()

sorted_table = table.sort(['cnr'], ascending=[0])



PLOT = 1
if PLOT:
    showMetric('cnr',sorted_table,'CNR',data_DIR+'/CNR.pdf')
    showMetric('dynamic_range',sorted_table,'Dynamic Range',data_DIR+ '/Dynamic_Range.pdf')
    showMeanSTDMetric('tubularity',sorted_table,'Tubularity',  data_DIR+'/Tubularity.pdf')




COMPARISON = 0
if COMPARISON:
  other_table = pd.DataFrame()
  all_table = pd.DataFrame()
  with open(data_DIR+'/../compare_list.csv','r') as f:
      i = 0
      for line in f:
            csv_path =  data_DIR + "/../"+ line.split("\n")[0]
            print csv_path
            df = pd.read_csv(csv_path)
            df = df[df['segment_type'] == -1]  # all-swcnode stats
             
            df['swc_file'] = csv_path[0:-4]
            df['specimen_id'] = i
            frames = [other_table ,df]
            other_table = pd.concat(frames,ignore_index = True)
            i = i + 1

  frames = [other_table, sorted_table]
  all_table = pd.concat(frames,ignore_index=True)

  showMetric('cnr',all_table,'CNR', data_DIR+'/compare_all_cnr.pdf')



