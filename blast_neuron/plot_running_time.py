__author__ = 'xiaoxiaol'



import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd


time_csv="/data/mat/xiaoxiaol/data/reconstructions_2015_1207/auto_recons/running_time.csv"
df_time=pd.read_csv(time_csv)
#image,plugin,running_time
algorithm_time_csv ="/data/mat/xiaoxiaol/data/reconstructions_2015_1207/ ported_neuron_tracing_spreadsheet.csv"
df_check_table = pd.read_csv(algorithm_plugin_csv)