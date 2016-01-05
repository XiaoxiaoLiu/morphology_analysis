#!/usr/bin/env python


import os
import random
import subprocess, threading
import platform


V3D=""
if (platform.system() == "Linux"):
    V3D= "/local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d"
else:
    V3D="/Users/xiaoxiaoliu/work/v3d/v3d_external/bin/vaa3d64.app/Contents/MacOS/vaa3d64"

QMasterV3D = "/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters/vaa3d"


class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            #print 'Thread started'
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            #print 'Thread finished'

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print 'Terminating process'
            print self.cmd
            self.process.terminate()
            thread.join()
        #print self.process.returncode




def gen_qsub_script(cmd, job_name, script_fn):
# example.qsub
# ## There are several queues available.  Please check with !ITsupport to verify which queue you should use
# #PBS -q mindscope
# # Declare that your job will use no more than some amount of memory _at_peak_
# #PBS -l vmem=16g
# # Allow up to 2 hours of walltime.  Default is 12 hours
# #PBS -l walltime=00:30:00
# # Request just one core on the host
# #PBS -l ncpus=1
# # Give your job a descriptive name. This is visible in qstat and other job reports.  Also serves as the default basename for log files
# #PBS -N /data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted/checked6_human_culturedcell_Cambridge_in_vitro_confocal_GFP/Series009/Series009.v3dpbd_x216_y266_z36_app2.swc
# # Should torque automatically re-run the job on error?
# #PBS -r n
# # Merge STDOUT into STDERR
# #PBS -j oe
# # location for stderr/stdout log files _after_ job completion
# #
# #
# export DISPLAY=:13918
# Xvfb :13918 -auth /dev/null &
# export LD_LIBRARY_PATH=/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters
# /data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters/vaa3d -x neuron_distance -f neuron_distance -i /data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted/checked6_human_culturedcell_Cambridge_in_vitro_confocal_GFP/Series009/Series009.v3dpbd_x216_y266_z36_app2.swc  /data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted/checked6_human_culturedcell_Cambridge_in_vitro_confocal_GFP/Series009/Series009.v3dpbd_x216_y266_z36_app2.swc -o /data/mat/xiaoxiaol/work/testqmaster.log
# kill %1
    output_dir = os.path.dirname(script_fn)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    FILE = open(script_fn, 'w')
    FILE.write("#PBS -q mindscope\n")
    FILE.write("#PBS -l vmem=16g\n")
    FILE.write("#PBS -l walltime=00:30:00\n")
    FILE.write("#PBS -l ncpus=1\n")
    FILE.write("#PBS -N %s\n" % job_name)
    FILE.write("#PBS -r n\n")
    FILE.write("#PBS -j oe\n")
    random_number = random.randint(5000,10000)
    FILE.write("export DISPLAY=:%s\n" % random_number)
    FILE.write("Xvfb :%s -auth /dev/null & \n" % random_number)
    FILE.write("export LD_LIBRARY_PATH=/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters\n")
    FILE.write("%s\n" % cmd)
    FILE.write("kill %1\n")

    FILE.close()




#===================================================================
def sort_swc(inputswc_path, outputswc_path, GEN_QSUB = 0, qsub_script_dir= "."):
    output_dir = os.path.dirname(outputswc_path)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    arguments = "  -x sort_neuron_swc -f sort_swc -i "+ inputswc_path + " -o " + outputswc_path + " -p 0 "

    if "spanningtree" in inputswc_path:# spanningtree algirhtm generate swc files without root type!!! the first one is used as the root
        arguments = arguments + " 1"

    if GEN_QSUB :
        cmd = QMasterV3D + arguments
        print cmd
        script_fn = qsub_script_dir +'/'+inputswc_path.split('/')[-1]+'.qsub'
        jobname = qsub_script_dir+inputswc_path.split('/')[-1]
        gen_qsub_script(cmd, jobname, script_fn)
    else:
        cmd = V3D + arguments
        print cmd
        command = Command(cmd)
        command.run(timeout=60*10)

    return

def consensus(input_ano_path, output_eswc_path, method=1, GEN_QSUB = 0, qsub_script_dir= "."):
    output_dir = os.path.dirname(output_eswc_path)
    logfile = output_eswc_path+'.log'
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    arguments = " -x consensus_swc -f consensus_swc -i " + input_ano_path + " -o " + output_eswc_path + " -p "+ str(method)+" >"+logfile

    if GEN_QSUB :
        cmd = QMasterV3D + arguments
        #print cmd
        script_fn = qsub_script_dir +'/'+input_ano_path.split('/')[-1]+'.qsub'
        jobname = qsub_script_dir+input_ano_path.split('/')[-1]
        gen_qsub_script(cmd, jobname, script_fn)
    else:
        cmd = V3D + arguments
        print cmd
        command = Command(cmd)
        command.run(timeout=60*5)
    return


def median_swc(input_ano_path, output_swc_path, GEN_QSUB = 0, qsub_script_dir= "."):
    output_dir = os.path.dirname(output_swc_path)
    logfile = output_swc_path+'.log'
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    arguments = " -x consensus_swc -f median_swc -i " + input_ano_path + " -o " + output_swc_path + " >"+logfile

    if GEN_QSUB :
        cmd = QMasterV3D + arguments
        #print cmd
        script_fn = qsub_script_dir +'/'+input_ano_path.split('/')[-1]+'.qsub'
        jobname = qsub_script_dir+input_ano_path.split('/')[-1]
        gen_qsub_script(cmd, jobname, script_fn)
    else:
        cmd = V3D + arguments
        #print cmd
        command = Command(cmd)
        command.run(timeout=60*3)
    return

def vote_map(input_ano_path, output_img_path, GEN_QSUB = 0, qsub_script_dir= "."):
    output_dir = os.path.dirname(output_img_path)
    logfile = output_img_path+'.log'
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    arguments = " -x consensus_swc -f vote_map -i " + input_ano_path + " -o " + output_img_path + " >"+logfile

    if GEN_QSUB :
        cmd = QMasterV3D + arguments
        #print cmd
        script_fn = qsub_script_dir +'/'+input_ano_path.split('/')[-1]+'.qsub'
        jobname = qsub_script_dir+input_ano_path.split('/')[-1]
        gen_qsub_script(cmd, jobname, script_fn)
    else:
        cmd = V3D + arguments
        print cmd
        command = Command(cmd)
        command.run(timeout=60*10)
    return




def resample(inputswc_path, outputswc_path,step_len = 1, GEN_QSUB = 0, qsub_script_dir= "."):
    output_dir = os.path.dirname(outputswc_path)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir
    arguments = " -x resample_swc -f resample_swc -i " + inputswc_path + " -o " + outputswc_path + " -p "+ str(step_len) + " >tmp.log"

    if GEN_QSUB :
        cmd = QMasterV3D + arguments
        #print cmd
        script_fn = qsub_script_dir +'/'+inputswc_path.split('/')[-1]+'.qsub'
        jobname = qsub_script_dir+inputswc_path.split('/')[-1]
        gen_qsub_script(cmd, jobname, script_fn)
    else:
        cmd = V3D + arguments
        #print cmd
        command = Command(cmd)
        command.run(timeout=60*10)
    return

def run_neuron_dist(inputswc_path1, inputswc_path2, logfile='./test.log',GEN_QSUB = 0, qsub_script_dir= "."):
    #Distance between neuron 1 /home/xiaoxiaol/work/data/test_frog2-2.swc and neuron 2 /home/xiaoxiaol/work/data/test_frog2-2.swc is:
    #entire-structure-average = 8.20009e-07
    #differen-structure-average = 0
    #percent of different-structure = 0

    # log file format
    # file1 file2   8.20009e-07  0 0

    arguments = " -x neuron_distance -f neuron_distance -i " + inputswc_path1 + " " + inputswc_path2 + " -o " + logfile + " >tmp.log"

    if GEN_QSUB :
        cmd = QMasterV3D + arguments
        print cmd
        script_fn = qsub_script_dir +'/'+str(random.randint(1000000,9999999))+'.qsub'
        jobname = qsub_script_dir+inputswc_path1.split('/')[-1]
        gen_qsub_script(cmd, jobname, script_fn)
    else:
        cmd = V3D + arguments
        #print cmd
        command = Command(cmd)
        command.run(timeout=60*5)
    return





def read_neuron_dist_log(logfile):
    # read log file
    f = open(logfile, 'r')
    line = f.readline()

    ave = float(line.split()[2])  # entire-structure-average
    diff = float(line.split()[3])  #differen-structure-average
    perc = float(line.split()[4])  #percent of different-structure
    return {'ave': ave, 'diff': diff, 'perc': perc}

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

    dist = read_neuron_dist_log(logfile)
    return dist



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


def pre_processing(inputswc_path, outputswc_path):  # prune align x-axis and resample
    output_dir = os.path.dirname(outputswc_path)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    #s =2  resampling size   #r = 0  --skip rotation
    cmd = V3D + " -x blastneuron -f pre_processing -p \"#i " + inputswc_path + "  #o " + outputswc_path + " #s 2 #r 0\" "
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




###########################################################


def genLinkerFile(swcDir, linker_file):
    input_swc_paths = [os.path.join(dirpath, f)
                       for dirpath, dirnames, files in os.walk(swcDir)
                       for f in files if f.endswith('.swc')]
    if len(input_swc_paths)> 0:
        with open(linker_file, 'w') as f:
            for input_swc_path in input_swc_paths:
                #print input_swc_path
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



WORK_PATH = "/local1/xiaoxiaol/work"
MRMR= WORK_PATH+"/src/mrmr_c_src/mrmr"
def selectFeatures_MRMR(df_all, feature_names,  threshold=0, number_of_features=10, selection_method='MID', data_DIR="."):
    #write out feature array into a csv file, then execute MRMR


    csvfile = data_DIR+"/zscore_for_mrmr.csv"
    np.savetxt(csvfile, featureArray, delimiter=",")
    # call MRMR
    cmd = MRMR +  " -i "+ csvfile + " -t "+ threshold + " -n " + number_of_features
    print cmd
    os.system(cmd)
    return

