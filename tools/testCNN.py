import __init__path
import os
import tensorflow as tf
import pandas as pd
import numpy as np
from lib.read.read_finalCpgFiles import cpgFile_DefaultExcludeAndKeep
from lib.read.read_data import load_data

import tensorflow as tf
import tensorflow.contrib.layers as lays
import tensorflow.contrib.slim as slim
from tensorflow.contrib.slim import losses
from tensorflow.contrib.slim import arg_scope

def dense_to_one_hot(labels_dense, num_classes):
    return np.eye(num_classes)[np.array(labels_dense).reshape(-1)]
def read_individual_image(img_path):
    f = open(img_path,'r').readlines()
    res = np.array([[int(element.strip()) for element in row.strip().split('\t')] for row in f])
    return res


class img(object):
    def __init__(self,imgname,img_dirpath):
        self._imgname = imgname
        self._dirpath = img_dirpath

    def imgname(self):
        return self._imgname

    def imgpath(self):
        return self._imgpath

    def img_data(self):
        self._imgpath = os.path.join(self._dirpath,self._imgname+'.txt')
        return read_individual_image(self._imgpath)


def find_correspondingImgs(item):
        cpg_img = img(item[0],cpgSites_dirpath).img_data()
        gene_img = img(item[1],geneSites_dirpath).img_data()
        return np.stack((cpg_img,gene_img),axis=-1)


class iterlist(object):
    def __init__(self,namearray,all_labels):
        # namearray here needs to be 'cpgname_probename', dataframe[['SNPName','ProbeName']]
        self._namearray = namearray
        self._all_labels = dense_to_one_hot(all_labels, 2)
        self._epochs_completed = 0
        self._index_in_epoch = 0
        self._num_examples = namearray.shape[0]

    def all_labels(self):
        return self._all_labels

    def all_images(self):
        return np.array([find_correspondingImgs(item) for item in self._namearray])

    def next_batch(self,batch_size,shuffle=True):
        start = self._index_in_epoch
        # shuffle for the first epoch
        if self._epochs_completed == 0 and start == 0 and shuffle:
            perm0 = np.arange(self._num_examples)
            np.random.shuffle(perm0)
            self._names = self._namearray[perm0]
            self._labels = self._all_labels[perm0]
        if start + batch_size > self._num_examples:
            # finish epoch
            self._epochs_completed += 1
            # get the rest examples in this epoch
            rest_num_examples = self._num_examples - start
            names_rest_part = self._namearray[start:self._num_examples]
            labels_rest_part = self._all_labels[start:self._num_examples]
            # shuffle the dataset
            if shuffle:
                perm = np.arange(self._num_examples)
                np.random.shuffle(perm)
                self._names = self._namearray[perm]
                self._labels = self._all_labels[perm]
            # start the new epoch
            start = 0
            self._index_in_epoch = batch_size - rest_num_examples
            end = self._index_in_epoch
            names_new_part = self._namearray[start:end]
            labels_new_part = self._all_labels[start:end]
            return_names = np.concatenate((names_rest_part,names_new_part),
                                          axis=0)
            return_labels = np.concatenate((labels_rest_part,labels_new_part),
                                           axis=0)
            return [np.array([find_correspondingImgs(item) for item in return_names]),
                    return_labels]
        else:
            self._index_in_epoch += batch_size
            end = self._index_in_epoch
            self._data = np.array([find_correspondingImgs(item) for item in self._names[start:end]])
            self._labels = self._all_labels[start:end]
            return [self._data,self._labels]


def read_namelist(namelist_filepath):
    with open(namelist_filepath,'r') as f:
        allnames = f.readlines()
        return np.array([img(col.strip(),cpgSites_dirpath) for col in allnames])


def test_cnn(x_image):
    def weight_variable(shape):
        initial = tf.truncated_normal(shape, stddev=0.1)
        return tf.Variable(initial)

    def bias_variable(shape):
        initial = tf.constant(0.1, shape=shape)
        return tf.Variable(initial)

    def conv2d(x, W):
        return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

    def max_pool_2x2(x):
        return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                            strides=[1, 2, 2, 1], padding='SAME')

    # first layer
    W_conv1 = weight_variable([5, 5, 2, 32])
    b_conv1 = bias_variable([32])
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool_2x2(h_conv1)
    # second layer
    W_conv2 = weight_variable([5, 5, 32, 32])
    b_conv2 = bias_variable([32])
    h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
    h_pool2 = max_pool_2x2(h_conv2)

    # fc layer
    W_fc1 = weight_variable([8192,512])
    b_fc1 = bias_variable([512])
    h_pool2_flat = tf.reshape(h_pool2, [-1, 8192])
    h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

    # drop out layer
    keep_prob = tf.placeholder(tf.float32)
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # readout layer
    W_fc2 = weight_variable([512, 2])
    b_fc2 = bias_variable([2])
    y_conv = tf.matmul(h_fc1_drop, W_fc2) + b_fc2

    return y_conv,keep_prob


if __name__=='__main__':
    project_rootdir = '/home/shuang/projects/development_eqtm'
    data_dirpath = os.path.join(project_rootdir, 'data')
    cpgSites_dirpath = os.path.join(data_dirpath, 'eqtmZscores', 'testCNN', 'cpgSites')
    geneSites_dirpath = os.path.join(data_dirpath, 'eqtmZscores', 'testCNN', 'geneSites')
    eqtm_path = os.path.join(data_dirpath, 'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')
    namelist_dirpath = os.path.join(data_dirpath, 'eqtmZscores', 'testCNN', 'namelist')


    # cpg_namelist_filepath = os.path.join(namelist_dirpath,'cpgnames.txt')
    # cpg_namelist = read_namelist(cpg_namelist_filepath)
    # gene_namelist_filepath = os.path.join(namelist_dirpath,'genenames.txt')
    # gene_namelist = read_namelist(gene_namelist_filepath)
    exclude = ['SNPChr','PValue','SNPChrPos','ProbeChr',
               'ProbeCenterChrPos','CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)','FoldChange', 'FDR','checkChr',
               'TssSite', 'chr','SNPName_ProbeName']
    keep = ['SNPName', 'ProbeName', 'OverallZScore']
    data = load_data(eqtm_path, keep=keep, exclude=exclude)
    test = iterlist(data.train.values[['SNPName', 'ProbeName']].values,data.train.labels)
    for ep in range(5):
        print(test.next_batch(50)[0].shape,test.next_batch(50)[1].shape)
    # build the model
    x = tf.placeholder(tf.float32, shape=[None, 32, 127, 2])
    y_ = tf.placeholder(tf.float32, shape=[None, 2])
    y_conv,keep_prob = test_cnn(x)
    cross_entropy = tf.reduce_mean(
        tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y_conv))
    train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)
    correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

    # start training
    train_data = iterlist(data.train.values[['SNPName', 'ProbeName']].values,
                          data.train.labels)
    test_data = iterlist(data.test.values[['SNPName', 'ProbeName']].values,
                         data.test.labels)

    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        for i in range(20000):
            batch = train_data.next_batch(10)
            if i % 100 == 0:
                train_accuracy = accuracy.eval(feed_dict={x: batch[0],
                                                          y_: batch[1],
                                                          keep_prob: 1.0})
                print('step %d, training accuracy %g' % (i, train_accuracy))
            train_step.run(feed_dict={x: batch[0],
                                      y_: batch[1],
                                      keep_prob: 0.5})
        print(accuracy.eval(feed_dict={x: test_data.all_images(),
                                       y_: test_data.all_labels(),
                                       keep_prob: 1.0}))
