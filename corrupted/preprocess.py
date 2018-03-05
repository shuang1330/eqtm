import csv
import os
import gzip

folder = os.path.join('ori_files','dna_files')
for file_name in os.listdir(folder):
    file_dir = os.path.join(folder,file_name)
    print(file_dir)
    with gzip.open(file_dir,'r') as f:
        for line in f.readlines():
            print(line.split('\t'))
        # reader = csv.reader(f.read(),delimiter='\t')
        # for row in reader:
            # print(row)
