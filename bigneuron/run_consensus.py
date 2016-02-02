import pandas as pd
import os
import sys
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import pandas as pd
import numpy as np
import os
import blast_neuron.blast_neuron_comp as bn


data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver"

output_dir = "/data/mat/xiaoxiaol/data/big_neuron/silver/gold_163_all_somasorted"
os.system("mkdir "+output_dir)




