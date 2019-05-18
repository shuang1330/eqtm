import shlex, subprocess
import os


def callBedtools(findOverlap_path,
                 project_rootdir,
                 feature_dirpath,
                 feature_list,
                 feature_filetype,
                 input_bedtoolFormat_filepath,
                 input_filename):
    print("First check how many existing intersection files.")
    saved_feature_dirpath = os.path.join(project_rootdir, 'data', 'temp',
                       'intersect_features', input_filename)
    print(saved_feature_dirpath)
    if os.path.exists(saved_feature_dirpath):
        processed_features = set(['.'.join(name.split(".")[:-2]) for name in os.listdir(saved_feature_dirpath)])
        if len(processed_features) > 0:
            processing_feature_list = [item for item in feature_list if
                                       item not in processed_features]
            print("Continue creating intersection files, %d files left to be processed."%len(processing_feature_list))
    else:
        print("None is processed, start to process %d features."%len(feature_list))
        processing_feature_list = feature_list
    print("Processing input file %s" % input_bedtoolFormat_filepath, flush=True)
    for feature in processing_feature_list:
        feature_filename = feature + ".imputed.gappedPeak.bed.gPk.gz"
        subprocess.check_call([findOverlap_path,
                               '%s' % project_rootdir,
                               '%s' % feature_dirpath,
                               '%s' % feature_filename,
                               '%s' % feature,
                               '%s' % input_bedtoolFormat_filepath,
                               '%s' % input_filename])
        print('Feature {} processed.'.format(feature))
