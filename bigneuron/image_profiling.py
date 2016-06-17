__author__ = 'coriannaj'

###

import os
import subprocess, threading
import pandas as pd

V3D =  "/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters/start_vaa3d.sh"

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

subfolder="0401_gold163_all_soma_sort"
data_DIR="/local1/home/coriannaj/Desktop/"+subfolder

def RUN_Vaa3d_Job(arguments):

    # run in local python env
    cmd = V3D + arguments
    print cmd
    command = Command(cmd)
    command.run(timeout=60*10)
    return

def profiling(input_img, input_swc, output_file, dilation_ratio = 3, flip = 0, invert = 0, cutoff_ratio=0.05):

    arguments = " -x profiling -f profile_swc -i "+input_img+" "+input_swc+" -o "+output_file+" -p "+str(dilation_ratio)+" "+str(flip)+"  "+str(invert)+" "+str(cutoff_ratio)

    RUN_Vaa3d_Job(arguments)
    return

def runGold(out_file):

    #create df table with headers
    cols = ['Image Number', 'CNR', 'Dynamic Range', 'Mean FG', 'Mean BG', 'Mean Tubularity']
    overall_profile = pd.DataFrame(columns=cols)

    print data_DIR

    folder_num = 0

    # go through all image directories in gold163
    for dirName, subdirList, fileList in os.walk(data_DIR, topdown=False):
        if dirName != data_DIR:
            if not os.path.basename(dirName).isdigit():
                #if not an image folder do nothing
                continue

            subfolder_path = os.path.join(data_DIR, dirName)
            folder_num += 1

            swc = None
            img = None

            profile_path = os.path.join(subfolder_path, 'profile.csv')

            for f in fileList:

                #checking all files in directory for img and trace
                if f.endswith('.swc'):
                    swc = os.path.join(subfolder_path, f)
                if f.endswith(('.v3dpbd','.v3draw')):
                    img = os.path.join(subfolder_path, f)

            #create profile if all files necessary were found
            if swc != None and img != None:
                profiling(img, swc, profile_path)

                #read in CSV output file
                try:
                    profile_df = pd.read_csv(profile_path)
                except IOError:
                    print "Folder %d did not create a profile.csv" %folder_num
                else:
                    #if file exists get data
                    stats = [int(os.path.basename(dirName)), profile_df.at[0,'cnr'], profile_df.at[0,'dynamic_range'], profile_df.at[0,'fg_mean'], profile_df.at[0,'bg_mean'], profile_df.at[0,'tubularity_mean']]
                    print stats
                    overall_profile.loc[len(overall_profile)]=stats

        else: #topdown = FALSE so will happen after all subdirectories explored
            print "exporting csv"
            print "visited %d folders" %folder_num
            print overall_profile
            overall_profile.sort_values(by='Image Number', inplace=True)
            overall_profile.to_csv(out_file, index=False)
    return

runGold("/local1/home/coriannaj/Desktop/0401_gold163_all_soma_sort/img_profiling.csv")

#to add additional folders to data compilation
def runAddGold(data_dir, out_file, images):
    out_path = os.path.join(data_DIR, out_file)

    data_df = pd.read_csv(out_path)

    for img_d in images:
        data_path = os.path.join(data_dir, img_d)
        files = os.listdir(data_path)
        for f in files:
            profile_path = os.path.join(data_path, 'profile.csv')

            if f.endswith('.swc'):
                swc = os.path.join(data_path, f)
            if f.endswith(('.v3dpbd', '.v3draw')):
                img = os.path.join(data_path, f)

        if swc != None and img != None:
            profiling(img, swc, profile_path)

            try:
                profile_df = pd.read_csv(profile_path)
            except IOError:
                print "Folder %s did not create a profile.csv" %img_d
            else:
                # if file exists get data
                stats = [158, img_d, profile_df.at[0, 'cnr'], profile_df.at[0, 'dynamic_range'],
                         profile_df.at[0, 'fg_mean'], profile_df.at[0, 'bg_mean'], profile_df.at[0, 'tubularity_mean']]
                print stats
                data_df.loc[len(data_df)] = stats

    data_df.to_csv(out_file)

#server_data = '/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort/'

#runAddGold(server_data, "/local1/home/coriannaj/Desktop/0401_gold163_all_soma_sort/img_profiling.csv", ["292", "293"])
