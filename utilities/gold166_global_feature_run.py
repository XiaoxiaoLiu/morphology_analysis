import os
import scipy.stats
import glob
import numpy
import matplotlib.pylab as pl
import pandas as pd


# program path on this machine
#===================================================================
#BlastNeuron Plugin: rapid global retrieval and local alignent of 3D neuron morphologies. This plugin also includes pre-processing, inverse projection and batch feature computation. by Yinan Wan.
#
#Functions:
#	 pre_processing          prune, resample and align neuron swcs.
#	 global_retrieve         given a query neuron, retrieve morphologically similar neurons from a neuron feature database.
#	 local_alignment         point to point alignment of two morphologically similar neurosn.
#	 batch_compute           generate neuron feature database (.nfb) file first to speed up neuron retrieve
#	 invert_projection       find an optimal affine transform between neuron structures regardless of topology.
#
#Example: 
#	 vaa3d -x blastneuron -f pre_processing -p  "#i input.swc #o result.swc #s 2 #r 1 "
#	 vaa3d -x blastneuron -f batch_compute -p "#i mydatabase.ano #o features.nfb"
#	 vaa3d -x blastneuron -f global_retrieve -p "#d features.nfb #q query.swc #n 10 #m 1,2,4 #o retrieve.ano"
#	 vaa3d -x blastneuron -f local_alignment -p "#t target.swc #s subject.swc #o match.swc"
#	 vaa3d -x blastneuron -f inverse_projection -p "#i target.swc #o target_invp.swc"


#V3D="/Users/xiaoxiaoliu/work/v3d/v3d_external/bin/vaa3d64.app/Contents/MacOS/vaa3d64"
V3D="/local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d"

data_DIR = "/data/mat/xiaoxiaol/data/gold166/gold166"





#===================================================================
def pre_processing(inputswc_path, outputswc_path):
    output_dir = os.path.dirname( outputswc_path )
    if not os.path.exists (output_dir) :
        os.system("mkdir -p  "+ output_dir)
        print "create output dir: ",output_dir
    cmd = V3D +  " -x blastneuron -f pre_processing -p \"#i "+ inputswc_path +"  #o "+outputswc_path+  " #s 2 #r 1\" "
    print cmd
    os.system(cmd)
    return


def batch_compute(input_linker_fn, feature_fn):

    cmd =  V3D + " -x blastneuron  -f batch_compute -p \" #i "+input_linker_fn + " #o " + feature_fn + "\" "
    os.system(cmd)
    print cmd
    return
    
    
    
def global_retrieve(inputswc_fn, feature_fn, result_fn, retrieve_number, logfile):
   # cmd = V3D + " -x blastneuron  -f global_retrieve  -p " +  "#d " + feature_fn + " #q " +inputswc_fn + " #n "+ \
   #       str(retrieve_number) +" -o "+result_fn+" -m 1,3" + " >" + logfile
    print cmd
    os.system(cmd)
    return
    
    

def genLinkerFile(swcDir, linker_file):
    
    input_swc_paths = [os.path.join(dirpath, f)
        for dirpath, dirnames, files in os.walk(swcDir)
        for f in files if f.endswith('.swc')]

    with open(linker_file, 'w') as f:
         for input_swc_path in input_swc_paths:
             print input_swc_path
             line = "SWCFILE="+input_swc_path+'\n'
             f.write(line)
    f.close()
    return

    
    
    
#--------------------------------------------------------
def removeLinkerFilePath(inputLinkFile, outputLinkFile):
    with open(outputLinkFile, 'w') as out_f:
        with open (inputLinkFile,'r') as in_f:
            for inline in in_f:
                outline = 'SWCFILE=' + inline.split('/')[-1]
                out_f.write(outline)
        in_f.close()
    out_f.close()
    return


def genLinkerFileFromList(listCSVFile, linkFile):
    df = pd.read_csv(listCSVFile, sep=',',header=0)
    fns = df.orca_path
    with open(linkFile, 'w') as f:
      for i in range(len(fns)):
          line = "SWCFILE="+fns[i]+'\n'
          f.write(line)
    f.close()
    return



#==================================================================================================
def main():
    if  not os.path.exists(data_DIR+'/preprocessed'):
        os.mkdir(data_DIR+'/preprocessed')

    preprocessed_dir = data_DIR +"/../gold166_preprocessed"

    for input_swc_path in glob.glob(data_DIR+"/*/*/*.swc"):
       print input_swc_path
       if( os.path.getsize(input_swc_path) > 1000  ):
            swc_fn = "/".join (input_swc_path.split("/")[-3:])
            preprocessed_swc_fn = preprocessed_dir+ '/'+swc_fn
    #        pre_processing(input_swc_path, preprocessed_swc_fn)

    preprocessed_ANO = preprocessed_dir+"/preprocessed.ano"
    #genLinkerFile( preprocessed_dir, preprocessed_ANO)

    ##batch computing
    feature_file =  preprocessed_dir+ "/features.nfb"
    batch_compute( preprocessed_ANO,feature_file)



if __name__ == "__main__":
      main()
