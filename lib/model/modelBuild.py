from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn import metrics


class convModelObject(object):

    def __init__(self, name, pipeline, param_grid, savefilename=None):
        self.modelname = name
        self.pipeline = pipeline
        self.param_grid = param_grid
        self.classifier = None
        self.savefilename = savefilename

    @staticmethod
    def calculate_output(classifier, test_data):
        """
        get confusion matrix and auc score for test dataset

        :param classifier:
        :param test_data:
        :return:
        """
        pred = classifier.predict(test_data.values)
        tn, fp, fn, tp = confusion_matrix(test_data.labels,pred).ravel()
        sensitivity = tp/(fn+tp)
        specificity = tn/(fp+tn)
        prods = classifier.predict_proba(test_data.values)[:,1]
        fpr, tpr, _ = metrics.roc_curve(test_data.labels,prods)
        score = metrics.auc(fpr,tpr) #auc score
        return round(sensitivity,2), round(specificity,2), round(score,2)

    def train_and_display_output(self, train_data, valid_data,
                                 returnRes=False, display=True):
        """
        (optional) save the model in path given in savefilename
        :param train_data:
        :param valid_data:
        :param returnRes:
        :param display:
        :return:
        """

        # train the model and finetune the hyperparameters
        print('Start training...')
        # print("test GridSearchCV")
        self.classifier = GridSearchCV(estimator=self.pipeline,
                                       param_grid=self.param_grid,
                                       n_jobs=8)
        self.classifier.fit(train_data.values, train_data.labels)
        sensitivity, specificity, score = self.calculate_output(self.classifier, valid_data)
        if self.savefilename is not None:
            pickle.dump(self.classifier, open(self.savefilename, 'wb'))
            print('Saved model to path:', self.savefilename)
        if display:
            print('Model Description:\n', self.classifier.best_estimator_)
            print('>>> best model results: sensitivity: {:.{prec}}\tspecificity: {:.{prec}f}\tauc:{}'.\
                  format(sensitivity, specificity, score, prec=3))
        if returnRes:
            return sensitivity, specificity, score
