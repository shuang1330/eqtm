import __init__path
import os
from lib.read.read_for_cnn import read_individual_image

if __name__=='__main__':
    data_dirpath = '/home/shuang/projects/development_eqtm/data/cpgSites/seperate_cpgFiles'
    check_dir = {}
    for filename in os.listdir(data_dirpath):
        filepath = os.path.join(data_dirpath, filename)
        image = read_individual_image(filepath)
        shape = '\t'.join([str(item) for item in image.shape])
        if shape not in check_dir:
            check_dir[shape] = 1
        else:
            check_dir[shape] += 1
    print(check_dir)
