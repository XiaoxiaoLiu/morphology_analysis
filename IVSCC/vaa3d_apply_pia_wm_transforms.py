# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 18:13:25 2015

@author: xiaoxiaol
"""

import pandas as pd
import os


def generateTransformFilesfromCSV(data_DIR, csv_file):
    df = pd.read_csv(csv_file)
    num_samples, num_cols = df.shape
    # transform = np.array([ [ float(a3d.find('tvr-00').text), float(a3d.find('tvr-01').text), float(a3d.find('tvr-02').text), float(a3d.find('tvr-09').text) ],
    # [ float(a3d.find('tvr-03').text), float(a3d.find('tvr-04').text), float(a3d.find('tvr-05').text), float(a3d.find('tvr-10').text) ],
    #                           [ float(a3d.find('tvr-06').text), float(a3d.find('tvr-07').text), float(a3d.find('tvr-08').text), float(a3d.find('tvr-11').text) ],
    for i in range(num_samples):
        orca_path = df['orca_path'][i]
        fn = orca_path.split('/')[-1]
        text_file = open(data_DIR + "/%s.txt" % fn.split('.swc')[0], "w")

        SCALE = 1000  # from  mm to microns

        text_file.write("%f %f %f %f \n" % (
        SCALE * df['tvr_00'][i], SCALE * df['tvr_01'][i], SCALE * df['tvr_02'][i], SCALE * df['tvr_09'][i]))

        text_file.write("%f %f %f %f \n" % (
        SCALE * df['tvr_03'][i], SCALE * df['tvr_04'][i], SCALE * df['tvr_05'][i], SCALE * df['tvr_10'][i]))

        text_file.write("%f %f %f %f \n" % (
        SCALE * df['tvr_06'][i], SCALE * df['tvr_07'][i], SCALE * df['tvr_08'][i], SCALE * df['tvr_11'][i]))

        text_file.close()
    return


def applyTransformBySpecimenName(swc_path, transform_path, output_path):
    cmd = ' /local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d -x affine_transform_swc -f apply_transform  -i ' + swc_path + ' -o ' + output_path + ' -p ' + transform_path
    print cmd
    os.system(cmd)



def main(argv):

    #####################################################################################################
    data_dir = '/local1/xiaoxiaol/work/data/lims2/0923_pw_aligned'
    csv_file = data_dir + "/0923_pw_aligned.csv"  ## where the transform parameters are, obtained from lims
    origin_data_DIR = data_dir + "/original"
    ####################################################################################################



    #########
    # output dirs
    transform_DIR = data_dir + "/transforms"  ## where to store the transform.txt
    output_DIR = data_dir + "/pw_aligned"


    #########
    # generate transform txt files
    if not os.path.exists(transform_DIR):
        os.makedirs(transform_DIR)

    generateTransformFilesfromCSV(transform_DIR, csv_file)

    ##########
    # apply transforms to swc files
    if not os.path.exists(output_DIR):
        os.makedirs(output_DIR)

    df = pd.read_csv(csv_file)
    df = df[df.specimen_id != 464326095]  # this specimenid's alignments are wrong
    data_table = df.values
    num_samples, num_cols = data_table.shape


    for i in range(num_samples):
        orca_path = data_table[i][num_cols - 1]
        fn = orca_path.split('/')[-1]
        transform_fn = transform_DIR + "/%s.txt" % fn.split('.swc')[0]
        output_fn = output_DIR + "/" + fn
        applyTransformBySpecimenName(origin_data_DIR + '/' + fn, transform_fn, output_fn)

if __name__ == "__main__":
    main()