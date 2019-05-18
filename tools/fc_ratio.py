import __init__path
import os
import tensorflow as tf
from lib.read.read_roadmap_features import roadmap
import numpy as np
from lib.read.read_data import load_data
from lib.read.read_allCpgFile import cpgFile
from lib.read.read_for_cnn import dense_to_one_hot
import keras
from keras.layers import Input, Dense
from keras.models import Model
from keras import backend as K
import pandas as pd
from sklearn.ensemble import RandomForestClassifier


class iter_cpgs(object):
    def __init__(self, cpg_names, feature_vectors, labels):
        self.cpg_names = cpg_names
        self.feature_vectors = feature_vectors
        self.labels = labels
        self._num_examples = len(cpg_names)

        self._epochs_completed = 0
        self._index_in_epoch = 0

        self._train_names = []
        self._train_values = []
        self._train_labels = []

    def next_batch(self, batch_size, shuffle=True):
        start = self._index_in_epoch
        # shuffle for the first epoch
        if self._epochs_completed == 0 and start == 0 and shuffle:
            perm0 = np.arange(self._num_examples)
            np.random.shuffle(perm0)
            self._train_names = self.cpg_names[perm0]
            self._train_labels = self.labels[perm0]
            self._train_values = self.feature_vectors[perm0]
        if start + batch_size > self._num_examples:
            # finish epoch
            self._epochs_completed += 1
            # get the rest examples in this epoch
            rest_num_examples = self._num_examples - start
            names_rest_part = self.cpg_names[start:self._num_examples]
            labels_rest_part = self.labels[start:self._num_examples]
            values_rest_part = self.feature_vectors[start:self._num_examples]
            # shuffle the dataset
            if shuffle:
                perm = np.arange(self._num_examples)
                np.random.shuffle(perm)
                self._train_names = self.cpg_names[perm]
                self._train_labels = self.labels[perm]
                self._train_values = self.feature_vectors[perm]
            # start the new epoch
            start = 0
            self._index_in_epoch = batch_size - rest_num_examples
            end = self._index_in_epoch
            names_new_part = self.cpg_names[start:end]
            labels_new_part = self.labels[start:end]
            values_new_part = self.feature_vectors[start:end]
            return_names = np.concatenate((names_rest_part, names_new_part), axis=0)
            return_labels = np.concatenate((labels_rest_part, labels_new_part), axis=0)
            return_values = np.concatenate((values_rest_part, values_new_part), axis=0)
            return return_names, return_values.reshape([-1, 32, 1]), return_labels
        else:
            self._index_in_epoch += batch_size
            end = self._index_in_epoch
            self._train_names = self.cpg_names[start:end]
            self._train_values = self.feature_vectors[start:end]
            self._train_labels = self.labels[start:end]
            return self._train_names, self._train_values.reshape([-1, 32, 1]), self._train_labels


def auc(y_true, y_pred):
    ptas = tf.stack([binary_PTA(y_true,y_pred,k) for k in np.linspace(0, 1, 1000)],axis=0)
    pfas = tf.stack([binary_PFA(y_true,y_pred,k) for k in np.linspace(0, 1, 1000)],axis=0)
    pfas = tf.concat([tf.ones((1,)) ,pfas],axis=0)
    binSizes = -(pfas[1:]-pfas[:-1])
    s = ptas*binSizes
    return K.sum(s, axis=0)


def binary_PFA(y_true, y_pred, threshold=K.variable(value=0.5)):
    y_pred = K.cast(y_pred >= threshold, 'float32')
    # N = total number of negative labels
    N = K.sum(1 - y_true)
    # FP = total number of false alerts, alerts from the negative class labels
    FP = K.sum(y_pred - y_pred * y_true)
    return FP/N


# P_TA prob true alerts for binary classifier
def binary_PTA(y_true, y_pred, threshold=K.variable(value=0.5)):
    y_pred = K.cast(y_pred >= threshold, 'float32')
    # P = total number of positive labels
    P = K.sum(y_true)
    # TP = total number of correct alerts, alerts from the positive class labels
    TP = K.sum(y_pred * y_true)
    return TP/P


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm"
    test_eqtm_ratio_filepath = os.path.join(project_rootdir, "data",
                                            "eqtmZscores",
                                            "withExpressionTSSMethyCpgOverlapGene",
                                            "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    feature_list_filepath = os.path.join(project_rootdir,
                                         "data",
                                         "temp",
                                         "meta",
                                         "feature_list.txt")
    histone_list_filepath = os.path.join(project_rootdir,
                                         "data",
                                         "temp",
                                         "meta",
                                         "histone_list.txt")

    roadmap_features = roadmap()
    feature_list = roadmap_features.readFeatureList(feature_list_filepath)
    histone_list = roadmap_features.readHistoneList(histone_list_filepath)
    # test_eqtm_ratio = pd.read_csv(test_eqtm_ratio_filepath, index_col=0, sep="\t")

    test_eqtm_ratio = cpgFile(allCpg_filepath=test_eqtm_ratio_filepath)
    exclude, keep = cpgFile.cpgFile_DefaultExcludeAndKeep()
    data = load_data(test_eqtm_ratio.allCpg_filepath, keep=keep, exclude=exclude)

    train_cpgnames = data.train.values.index
    train_feature_vectors = data.train.values[histone_list]
    additional_features = ["expressionMean", "expressionVar",
                           "TssDistance", "methyMean", "methyVar"]
    train_feature_additional = data.train.values[additional_features]
    train_labels = dense_to_one_hot(data.train.labels, 2)
    iter_train_data = iter_cpgs(cpg_names=np.array(train_cpgnames.values),
                                feature_vectors=np.array(train_feature_vectors.values),
                                labels=np.array(train_labels))
    print("\n\n\n", train_cpgnames.shape, train_feature_vectors.shape, train_labels.shape)
    test_cpgnames = np.array(data.test.values.index.values)
    test_features = np.array(data.test.values[histone_list].values)
    test_labels = np.array(data.test.labels.values)

    # build the model
    # Headline input: meant to receive sequences of 100 integers, between 1 and 10000.
    # Note that we can name any layer by passing it a "name" argument.

    main_input = Input(shape=(37,), dtype='float32', name='main_input')

    # x = Dense(64, activation='relu')(main_input)
    x = Dense(64, activation='relu')(main_input)
    x = Dense(32, activation='relu')(x)
    x = Dense(8, activation='relu')(x)
    main_output = Dense(2, activation='sigmoid', name='main_output')(x)
    model = Model(inputs=[main_input], outputs=[main_output])
    print(model.summary())
    print(train_feature_vectors.head())
    print(train_feature_additional.head())
    print(train_labels[:5])
    train_allfeatures = pd.concat([train_feature_additional, train_feature_vectors], axis=1)
    print(train_allfeatures.shape)
    rms = keras.optimizers.RMSprop(lr=0.01, rho=0.9, epsilon=None, decay=0.0)
    model.compile(optimizer=rms, loss='binary_crossentropy', metrics=['acc'])
    model.fit([train_allfeatures], [train_labels],
              epochs=200, batch_size=128)

    # randomforest model
    ranfor_model = RandomForestClassifier()
    ranfor_model.fit(train_allfeatures, train_labels)

    # TODO: test dataset evaluation

    raise NotImplementedError
    print(x.shape)

    auxiliary_input = Input(shape=(5,), name='aux_input')
    x = keras.layers.concatenate([x, auxiliary_input])

    # We stack a deep densely-connected network on top
    x = Dense(4, activation='relu')(x)
    # x = Dense(4, activation='relu')(x)

    # And finally we add the main logistic regression layer
    main_output = Dense(2, activation='sigmoid', name='main_output')(x)
    model = Model(inputs=[main_input, auxiliary_input], outputs=[main_output])
    print(model.summary())
    print(train_feature_vectors.head())
    print(train_feature_additional.head())
    print(train_labels[:5])
    model.compile(optimizer='rmsprop', loss='binary_crossentropy', metrics=['acc'])
    model.fit([train_feature_vectors, train_feature_additional], [train_labels],
              epochs=50, batch_size=128)
