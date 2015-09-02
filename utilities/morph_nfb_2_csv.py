import numpy as np
import pandas as pd


GL_FEATURE_TAGS = np.array(
    ['num_nodes', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips', 'overall_width',
     'overall_height', 'overall_depth', 'average_diameter', 'total_length', 'total_surface', 'total_volume',
     'max_euclidean_distance', 'max_path_distance', 'max_branch_order', 'average_contraction', 'average fragmentation',
     'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])
GMI_FEATURE_TAGS = np.array(
    ['moment1', 'moment2', 'moment3', 'moment4', 'moment5', 'moment6', 'moment7', 'moment8', 'moment9', 'moment10',
     'moment11', 'moment12', 'moment13', 'avgR'])



# ===================================================================
def readDBFeatures(FEATURE_FILE):
    # TODO: detect nan values
    glf_featureList = []  # each row is a feature vector
    gmi_featureList = []
    swc_file_nameList = []
    with open(FEATURE_FILE, 'r') as  f:
        for fn_line in f:  # ignore the SWCFILE=* line

            swc_file = fn_line[8:].strip()
            swc_file_nameList.append(swc_file)

            line_globalFeature = (f.next()).strip()
            glf = map(float, line_globalFeature.split('\t'))
            glf_featureList.append(glf)

            line_GMI = (f.next()).strip()
            gmi = map(float, line_GMI.split('\t'))
            gmi_featureList.append(gmi)

    return swc_file_nameList, np.array(glf_featureList), np.array(gmi_featureList)


def generateALLFeatureCSV(feature_file, feature_csv_file):
    swc_file_nameList, glFeatures, gmiFeatures = readDBFeatures(feature_file)

    allFeatures = np.append(glFeatures, gmiFeatures, 1)
    allColumns = np.append(GL_FEATURE_TAGS, GMI_FEATURE_TAGS, 0)

    df = pd.DataFrame(allFeatures, columns=allColumns)

    df['swc_file'] = pd.Series(swc_file_nameList, index=df.index)

    algorithmList = []
    imageList = []
    for swc_file in swc_file_nameList:
        fn = swc_file.split('/')[-1]
        algorithm = fn.split('v3dpbd_')[-1]
        algorithm = algorithm.split('.')[0]
        if "app1" in algorithm:   # for patterns like *x245_y234_z234_app1.swc
              algorithm = "app1"
        if "app2" in algorithm:
              algorithm = "app2"
        if  "spanningtree" in algorithm: # fastmarching_spanningtree is too long
              algorithm = "spanningtree"
        image = fn.split('.v3dpbd')[0]
        image = image.split('sorted_')[-1]  # for sorted_* swc_files
        algorithmList.append(algorithm)
        imageList.append(image)

    df['algorithm'] = pd.Series(algorithmList, index=df.index)
    df['image'] = pd.Series(imageList, index=df.index)

    allColumns = np.append(np.array(['image', 'algorithm', 'swc_file']), allColumns, 0)

    df = df[allColumns]

    df.to_csv(feature_csv_file, index=False)

    print 'output all feature csv file to :', feature_csv_file
    return


def generateLinkerFileFromCSV(result_dir, csvfile, column_name):
    df = pd.read_csv(csvfile)
    types = df[column_name]
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file']
        with open(result_dir + '/' + atype + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                line = 'SWCFILE=' + afile + '\n'
                outf.write(line)
        outf.close()
    return

def normalizeFeatures(features):
    meanFeatures = np.mean(features, 0)
    stdFeatures = np.std(features, 0)
    if np.count_nonzero(stdFeatures) < len(stdFeatures):
        print "zero detected"
        print stdFeatures
    normalized = (features - meanFeatures) / stdFeatures
    return normalized


def concatCSVs(csv1, csv2, outcsv):
    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)

    #out_df = pd.merge(df1, df2)
    out_df = pd.concat([df1, df2], axis=1)
    out_df.to_csv(outcsv, index=False)
    return


def saveScoreMatrix(featureArray, scoreMatrix_file, REMOVE_OUTLIER=1):
    scoreMatrix = []
    normalized = normalizeFeatures(featureArray)

    # remove outliers!!!
    normalized[normalized < -3] = -3
    normalized[normalized > 3] = 3

    for i in range(len(normalized)):
        queryFeature = normalized[i]  # each row is a feature vecto
        #scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)
        scores = np.sum(np.abs(normalized - queryFeature) ** 2, 1)
        scoreMatrix.append(scores)

    df = pd.DataFrame(scoreMatrix)
    df.to_csv(scoreMatrix_file)

