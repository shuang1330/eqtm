import shlex, subprocess


def callBedtools(findOverlap_path,
                 project_rootdir,
                 feature_dirpath,
                 feature_list,
                 feature_filetype,
                 input_bedtoolFormat_filepath,
                 input_filename):
    print("Processing input file %s" % input_bedtoolFormat_filepath, flush=True)
    for feature in feature_list:
        feature_filename = feature + ".imputed.gappedPeak.bed.gPk.gz"
        subprocess.check_call([findOverlap_path,
                               '%s' % project_rootdir,
                               '%s' % feature_dirpath,
                               '%s' % feature_filename,
                               '%s' % feature,
                               '%s' % input_bedtoolFormat_filepath,
                               '%s' % input_filename])
        print('Built the overlap matrix for feature: %s' % feature)
