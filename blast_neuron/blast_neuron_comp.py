#!/usr/bin/env python


import os

# ===================================================================
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






import subprocess, threading

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            print 'Thread started'
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            print 'Thread finished'

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print 'Terminating process'
            self.process.terminate()
            thread.join()
        print self.process.returncode


#V3D="/Users/xiaoxiaoliu/work/v3d/v3d_external/bin/vaa3d64.app/Contents/MacOS/vaa3d64"
V3D= "/local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d"


#===================================================================
def sort_swc(inputswc_path, outputswc_path):
    output_dir = os.path.dirname(outputswc_path)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir
    cmd = V3D + " -x sort_neuron_swc -f sort_swc -i "+ inputswc_path + " -o " + outputswc_path + " -p 0 "
    print cmd
    #os.system(cmd)

    command = Command(cmd)
    command.run(timeout=60*5)

    return


def resample(inputswc_path, outputswc_path):
    output_dir = os.path.dirname(outputswc_path)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir
    cmd = V3D + " -x resample_swc -f resample_swc -i " + inputswc_path + " -o " + outputswc_path + " -p 1 "
    print cmd
    #os.system(cmd)
    command = Command(cmd)
    command.run(timeout=60*5)
    return


def pre_processing(inputswc_path, outputswc_path):  # prune align x-axis and resample
    output_dir = os.path.dirname(outputswc_path)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir
    cmd = V3D + " -x blastneuron -f pre_processing -p \"#i " + inputswc_path + "  #o " + outputswc_path + " #s 2 #r 1\" "
    print cmd
    os.system(cmd)
    return


def batch_compute(input_linker_fn, feature_fn):
    cmd = V3D + " -x blastneuron  -f batch_compute -p \" #i " + input_linker_fn + " #o " + feature_fn + "\" "
    os.system(cmd)
    print cmd
    return


def global_retrieve(inputswc_fn, feature_fn, result_fn, retrieve_number, logfile):
    # cmd = V3D + " -x blastneuron  -f global_retrieve  -p " +  "#d " + feature_fn + " #q " +inputswc_fn + " #n "+ \
    #       str(retrieve_number) +" -o "+result_fn+" -m 1,3" + " >" + logfile
    print cmd
    os.system(cmd)
    return


def neuron_dist(inputswc_path1, inputswc_path2, logfile='./test.log'):
    #Distance between neuron 1 /home/xiaoxiaol/work/data/test_frog2-2.swc and neuron 2 /home/xiaoxiaol/work/data/test_frog2-2.swc is:
    #entire-structure-average = 8.20009e-07
    #differen-structure-average = 0
    #percent of different-structure = 0

    # log file format
    # file1 file2   8.20009e-07  0 0

    cmd = V3D + " -x neuron_distance -f neuron_distance -i " + inputswc_path1 + " " + inputswc_path2 + " -o " + logfile
    os.system(cmd)
    print cmd

    # read log file
    f = open(logfile, 'r')
    line = f.readline()

    ave = float(line.split()[2])  # read the average difference number
    diff = float(line.split()[3])  # read the average difference number
    perc = float(line.split()[4])  # read the average difference number
    return {'ave': ave, 'diff': diff, 'perc': perc}


###########################################################


def genLinkerFile(swcDir, linker_file):
    input_swc_paths = [os.path.join(dirpath, f)
                       for dirpath, dirnames, files in os.walk(swcDir)
                       for f in files if f.endswith('.swc')]

    with open(linker_file, 'w') as f:
        for input_swc_path in input_swc_paths:
            print input_swc_path
            line = "SWCFILE=" + input_swc_path + '\n'
            f.write(line)
    f.close()
    return


#--------------------------------------------------------
def removeLinkerFilePath(inputLinkFile, outputLinkFile):
    with open(outputLinkFile, 'w') as out_f:
        with open(inputLinkFile, 'r') as in_f:
            for inline in in_f:
                outline = 'SWCFILE=' + inline.split('/')[-1]
                out_f.write(outline)
        in_f.close()
    out_f.close()
    return


def genLinkerFileFromList(listCSVFile, linkFile):
    df = pd.read_csv(listCSVFile, sep=',', header=0)
    fns = df.orca_path
    with open(linkFile, 'w') as f:
        for i in range(len(fns)):
            line = "SWCFILE=" + fns[i] + '\n'
            f.write(line)
    f.close()
    return





