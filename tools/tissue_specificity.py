import __init__path
import os
from sklearn import metrics
from sklearn.metrics import confusion_matrix
import pickle
from lib.read.read_data import load_data


def calculate_output(classifier, data, labels):
    pred = classifier.predict(data)
    tn, fp, fn, tp = confusion_matrix(labels, pred).ravel()
    sensitivity = tp / (fn + tp)
    specificity = tn / (fp + tn)
    prods = classifier.predict_proba(data)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(labels, prods)
    score = metrics.auc(fpr, tpr)  # auc score
    print(round(sensitivity, 2), round(specificity, 2), round(score, 2))


def examine(model_path, data_path, keep, exclude):
    model = pickle.load(open(model_path, "rb"))
    data = load_data(data_path, keep=keep, exclude=exclude, test_size=0)
    features = [col for col in data.train.values.columns if col not in keep]
    calculate_output(model, data.train.values[features], data.train.labels)
    return None


if __name__ == '__main__':
    project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    eqtm_dirpath = os.path.join(project_rootdir,
                                "data",
                                "eqtmZscores",
                                "withExpressionTSSMethyCpgOverlapGenePromoter")

    blood_et_filepath = os.path.join(eqtm_dirpath,
                                     "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap_PromoterOverlap.txt")
    non_blood_et_filepath = os.path.join(eqtm_dirpath,
                                         "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap_PromoterOverlap.txt")
    blood_gt_filepath = os.path.join(eqtm_dirpath,
                                     "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap_PromoterOverlap.txt")
    non_blood_gt_filepath = os.path.join(eqtm_dirpath,
                                         "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap_PromoterOverlap.txt")

    exclude = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr', 'SNPName_ProbeName',
               'TssSite', 'chr']
    keep = ['OverallZScore']
    model_dirpath = os.path.join(project_rootdir, "model")
    blood_et_model = os.path.join(model_dirpath, "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap_PromoterOverlap.txt.pkl")
    non_blood_et_model = os.path.join(model_dirpath,
                                      "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap_PromoterOverlap.txt.pkl")
    blood_gt_model = os.path.join(model_dirpath,
                                  "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap_PromoterOverlap.txt.pkl")
    non_blood_gt_model = os.path.join(model_dirpath,
                                      "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap_PromoterOverlap.txt.pkl")

    print("Et blood model.")
    print("Test on et non-blood.")
    examine(blood_et_model, non_blood_et_filepath, keep, exclude)
    print("Et non blood model.")
    print("Test on et blood.")
    examine(non_blood_et_model, blood_et_filepath, keep, exclude)
    print("Gt blood model.")
    print("Test on gt non-blood.")
    examine(blood_gt_model, non_blood_gt_filepath, keep, exclude)
    print("Gt non-blood model.")
    print("Test on gt blood.")
    examine(non_blood_gt_model, blood_gt_filepath, keep, exclude)
    print("Gt blood model.")
    print("Test on et non-blood")
    examine(blood_gt_model, non_blood_et_filepath, keep, exclude)
    print("Gt non blood model,")
    print("Test on et blood")
    examine(non_blood_gt_model, blood_et_filepath, keep, exclude)
