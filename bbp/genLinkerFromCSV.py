__author__ = 'xiaoxiaol'


def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True, fpath='.'):
    df = pd.read_csv(csvfile)
    if (column_name == None):
        swc_files = df['swc_file']
        with open(result_dir + '/all.ano', 'w') as outf:
            for afile in swc_files:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
        return

    types = df[column_name]
    print np.unique(types)
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file_name']
        with open(result_dir + '/' + atype + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + fpath +"/"+filename + '\n'
                outf.write(line)
            outf.close()
    return