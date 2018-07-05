import __init__path
from lib.model.finetune_model import examineModel
import os
import argparse

if __name__ == '__main__':
    # locally
    # project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon project root dir
    project_rootdir = '/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm'

    parser = argparse.ArgumentParser()
    parser.add_argument("--embedding_filename", dest="embedding_filename",
                        default=".")
    args = parser.parse_args()
    print("Processing embeddfing file: ", args.embedding_filename)
    # raise NotImplementedError

    eqtm_path = os.path.join(project_rootdir, 'data',
                             'features', 'embeddings',
                              args.embedding_filename+'.txt')
    model_name = 'ranfor'
    exclude = ['SNPName', 'ProbeName']
    keep = ['OverallZScore']
    train_path = eqtm_path
    test_path = train_path
    examineModel(model_name, train_path, test_path,
                 keep=keep, exclude=exclude, display=True)
