import __init__path
import os
from lib.model.finetune_model import examineModel
from lib.read.read_eqtm import cpgFile_DefaultExcludeAndKeep

if __name__=='__main__':
    data_dir = '/home/shuang/projects/development_eqtm/data/eqtmZscores'
    withOverlapTssMethy = os.path.join(data_dir,
    'withExpressionTSSMethyCpgOverlap',
    '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap.txt')
    withOverlapTssMethyExpressionGene = os.path.join(data_dir,
    'withExpressionTSSMethyCpgOverlapGene',
    '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')
    model_name='ranfor'
    exclude,keep = cpgFile_DefaultExcludeAndKeep()
    for filepath in [withOverlapTssMethy,withOverlapTssMethyExpressionGene]:
        print(filepath.split('/')[-1])
        train_path = filepath
        test_path = train_path
        examineModel(model_name,train_path,test_path,iteration=1,
                     keep=keep,exclude=exclude,display=True)
