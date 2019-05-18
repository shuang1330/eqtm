'''
shallow network
'''

import pandas as pd
import os

import tensorflow as tf
from math import ceil
import numpy as np
from lib.read_data import dataset,Datasets
from lib.net import autoencoder,feedforward_net
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from copy import copy
from sklearn import preprocessing


def dense_to_one_hot(labels_dense, num_classes):
  '''
  Convert class labels from scalars to one-hot vectors.
  '''
  return np.eye(num_classes)[np.array(labels_dense).reshape(-1)]

def read_data_set(data_table,test_size=0.25,normalization=True):
    '''
    convert a pandas dataframe data table into Datasets(dataset,dataset)
    '''
    train, test = train_test_split(data_table,test_size=0.25)
    train_x = train[[col for col in train.columns
                     if col not in ['zscore','direction','cpgName']]]
    features = train_x.columns
    if normalization:
        minMaxScaler = preprocessing.MinMaxScaler()
        train_x = minMaxScaler.fit_transform(train_x)
        test_x = minMaxScaler.fit_transform(test[[col for col in
                                                  train.columns
                      if col not in ['zscore','direction','cpgName']]])
    else:
        train_x = np.array(train_x)
        test_x = np.array(test[[col for col in train.columns
                      if col not in ['zscore','direction','cpgName']]])
    train_y = np.array(train['direction'],dtype=np.int8)
    test_y = np.array(test['direction'],dtype=np.int8)

    return Datasets(train=dataset(train_x,train_y),
                    test=dataset(test_x,test_y))


if __name__=='__main__':

    def load_bonderWestraData(path):
        data = pd.read_csv(path,sep=',')
        # print(data.head())
        def binarize(row):
            if row > 0:
                return 1
            else:
                return 0
        data['direction'] = data['zscore'].apply(binarize)
        dataset = read_data_set(data)
        return dataset


    bonder_path = 'data/bonder_withzscore.csv'
    westra_allFeat_path = 'data/westra_all_with_zscore.csv'
    westra_bonderFeat_path = 'data/westra_bonderfeat_with_zscore.csv'

    bonder = load_bonderWestraData(bonder_path)
    westra_allFeat = load_bonderWestraData(westra_allFeat_path)
    westra_bonderFeat = load_bonderWestraData(westra_bonderFeat_path)

    print('Bonder dataset loaded.',bonder.train.values.shape)
    print('Westra with all features dataset loaded.',westra_allFeat.train.values.shape)
    print('Westra with bonder features dataset loaded.',
           westra_bonderFeat.train.values.shape)

    # change eqtm into bonder/westra_allFeat/westra_bonderFeat
    eqtm = copy(bonder)


    train_fn = eqtm.train
    test_fn = eqtm.test

    # constant
    batch_size = 100
    epoch_num = 100
    lr = 0.01

    batch_per_ep = ceil(train_fn.num_examples/batch_size)


    # ================== feedforward shallow network ===========================
    # model
    inputs = tf.placeholder(tf.float32,(None, train_fn.num_features))
    labels = tf.placeholder(tf.float32,(None, 2))
    fn_outputs,fn_probs = feedforward_net(inputs)
    # loss and training options
    loss_fn = tf.reduce_mean(
                    tf.nn.weighted_cross_entropy_with_logits(
                    targets=labels,
                    logits=fn_outputs,
                    pos_weight=4))
    train_op = tf.train.AdamOptimizer(learning_rate=lr).minimize(loss_fn)

    # initializer
    init = tf.global_variables_initializer()

    # start training
    with tf.Session() as sess:
        sess.run(init)
        for ep in range(epoch_num):
            for batch_no in range(batch_per_ep):
                batch_data, batch_label = train_fn.next_batch(batch_size)
                batch_label_onehot = dense_to_one_hot(batch_label,2)
                _, probs,error = sess.run([train_op,
                                           fn_probs,loss_fn],
                                           feed_dict={inputs:batch_data,
                                           labels:batch_label_onehot})
                print('Epoch: {0}\tIteration:{1}\tError: {2}\t'.format(
                ep, batch_no, error
                ))

        # test
        batch_label_onehot = dense_to_one_hot(test_fn.labels,2)
        _,probs,_ = sess.run([fn_outputs,
                              fn_probs,
                              loss_fn],
                              feed_dict={inputs:test_fn.values,
                                         labels:batch_label_onehot})

        pred = probs[:,1]
        pred[pred>=0.5] = 1
        pred[pred<0.5] = 0

        tn, fp, fn, tp = confusion_matrix(test_fn.labels,pred).ravel()
        sensitivity = tp/(fn+tp)
        specificity = tn/(fp+tn)
        print(tn,fp,fn,tp)
        print(sensitivity,specificity)
