import __init__path
from lib.model.net import conv_1d
import tensorflow as tf
import os
import pandas as pd
from lib.read.read_roadmap_features import roadmap
import numpy as np
from lib.read.read_data import load_data
from lib.read.read_allCpgFile import cpgFile
from lib.read.read_for_cnn import dense_to_one_hot


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
    train_labels = dense_to_one_hot(data.train.labels, 2)
    iter_train_data = iter_cpgs(cpg_names=np.array(train_cpgnames.values),
                                feature_vectors=np.array(train_feature_vectors.values),
                                labels=np.array(train_labels))
    print("\n\n\n", train_cpgnames.shape, train_feature_vectors.shape, train_labels.shape)
    test_cpgnames = np.array(data.test.values.index.values)
    test_features = np.array(data.test.values[histone_list].values)
    test_labels = np.array(data.test.labels.values)
    # print(train_feature_vectors.shape, train_labels.shape)
    # iter_names, iter_features, iter_labels = iter_train_data.next_batch(100)
    # print(iter_names.shape, iter_features.shape, iter_labels.shape)
    # print("labels are: ", iter_labels)

    # build the model
    x_image = tf.placeholder(tf.float32, shape=[None, 32, 1], name="inputs")
    y_ = tf.placeholder(tf.float32, shape=[None, 2], name="output")
    lr = tf.placeholder(tf.float32, shape=[], name="learning_rate")
    y_conv = conv_1d(x_image)
    cross_entropy = tf.reduce_mean(
        tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y_conv))
    train_step = tf.train.AdamOptimizer(lr).minimize(cross_entropy)
    correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    batch_size = 200

    # start training
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    for iter in range(5000):
        batch_names, batch_features, batch_labels = \
            iter_train_data.next_batch(batch_size)
        sess.run(train_step, feed_dict={x_image: batch_features.reshape([-1, 32, 1]),
                                        y_: batch_labels,
                                        lr: 10e-4})
        if iter % 100 == 0:
            accuracy_ = accuracy.eval(session=sess,
                                      feed_dict={x_image: batch_features.reshape([-1, 32, 1]),
                                                 y_: batch_labels})
            print("Iter {}, accuracy {}".format(iter, accuracy_))
