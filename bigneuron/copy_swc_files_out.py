__author__ = 'xiaoxiaol'
import os
### main
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"
output_dir = data_DIR



destnation_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort_consensus3_strict_swc_only"


subdirs = [x[0] for x in os.walk(data_DIR)]
for recon_dir in subdirs[1:]:
        folder_name = recon_dir.split('/')[-1]
        if 'processed' in  folder_name:
                  subfolderpath=recon_dir.split('/')[-2]
                  print subfolderpath
                 # os.system('mkdir -p '+destnation_DIR+'/'+subfolderpath+'/processed' )
                 # os.system('cp '+ recon_dir+'/*.strict.swc  ' + destnation_DIR+'/'+subfolderpath+'/processed/')
                 # os.system('mkdir -p '+destnation_DIR+'/'+subfolderpath )
                  os.system('cp '+ recon_dir+'/../consensus3.strict.swc  ' + destnation_DIR+'/'+subfolderpath+'/')
