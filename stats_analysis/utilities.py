import os
import scipy
import numpy
import matplotlib.pylab as pl

#order_str= '2 22 28 50  7 21 42 19 44 36 43 46 34 39 17  1 16 15 30 26 12 13 25 33 38 41 48  8 49 35 27  4 11 31  9 45 37 10 32 14 29 23 40 47 18  3 20 24  5  6'

#order_str='2 22 28 50  7 21 42 19 44 36 43 46 34 39 17  1 16 15 30 26 12 13 25  3 33 38 41 48  8 49 35 27  4 11 31  9 45 37 10 32 14 29 23 20 18 24 40 47  5  6'


#order_str=' 22  2 28 21 50  7 33 42 19 30 26 12 13 25  3 43 46 34 39 17  1 36 38 32 14 44 35 49 48  4 27 11  8 18 40  9 23 37 47 10 41 20 16 15 24 29 31 45  5  6'



order_str=' 28 43 33 46 34 39 21 15 31 45 22 50 17  1  2  7  9 20 36 48 41 38 37 30 12 42 18 13 25  3  8 26 24 23 40 47 19 44 35 49 32 14 10 16 29 27 11  4  5  6'

order = order_str.split()

data_DIR = "/Users/xiaoxiaoliu/work/data/lims2/manual"

inputLinkerFile =  data_DIR+'/preprocessed/mylinker.ano'
outputLinkerFile =  data_DIR+'/preprocessed/hcluster_ordered_mylinker_16features_correlation.ano'



with open(inputLinkerFile,'r') as f:
   lines=f.readlines()

order = array(map(int, order)) -1
newLines = [ lines[i] for i in order[::-1] ]

with open(outputLinkerFile, 'w') as outf:
     outf.writelines(newLines)