import __init__path
import tensorflow as tf
import os
from lib.model.net import test_cnn
from lib.read.read_data import read_eqtm_names
from lib.read.read_for_cnn import iterlist
import pandas as pd
import argparse


if __name__=='__main__':
    # locally
    # project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon project root dir
    project_rootdir = '/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm'

    parser = argparse.ArgumentParser()
    parser.add_argument("--ckpt_filename", dest="ckpt_filename", type=str)
    parser.add_argument("--trainIter_batch_lr", dest="train_info", type=str)
    args = parser.parse_args()
    print("Processing ckpt file: ", args.ckpt_filename)
    # raise NotImplementedError

    # ckpt saved path
    ckpt_path = os.path.join(project_rootdir,
                             'data',
                             'model',
                             'test_cnn_model',
                             args.ckpt_filename+'.ckpt')
    # input data path
    cpgSites_dirpath = os.path.join(project_rootdir,
                                    'data',
                                    'cpgSites',
                                    'seperate_cpgFiles')
    eqtm_path = os.path.join(project_rootdir, 'data',
                             'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')
    eqtm_name = '_'.join(['2017-12-09-eQTLsFDR-gt0_withExpressionTssMethy', args.train_info])
    # save path
    save_embeddings_dirpath = os.path.join(project_rootdir,
                                           'data',
                                           'features',
                                           'embeddings')

    # build the model
    x = tf.placeholder(tf.float32, shape=[None, 127, 32, 1])
    y_ = tf.placeholder(tf.float32, shape=[None, 2])
    y_conv, keep_prob, embedding = test_cnn(x)
    cross_entropy = tf.reduce_mean(
        tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y_conv))
    train_step = tf.train.AdamOptimizer(0.001).minimize(cross_entropy)
    correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

    # load back model
    saver = tf.train.Saver()

    # start the session
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    saver.restore(sess, ckpt_path)
    print("Loaded ckpt files.", flush=True)

    # load input data
    def read_embeddings(sess, eqtm_path, features_to_save):
        [train_values, train_labels, test_values, test_labels] = read_eqtm_names(eqtm_path)
        train_data = iterlist(train_values[['SNPName', 'ProbeName']].values, train_labels, cpgSites_dirpath)
        test_data = iterlist(test_values[['SNPName', 'ProbeName']].values, test_labels, cpgSites_dirpath)
        # get embeddings for data
        train_embeddings = embedding.eval(session=sess,
                                          feed_dict={x: train_data.all_images(),
                                                     y_: train_data.all_labels(),
                                                     keep_prob: 1.0})
        test_embeddings = embedding.eval(session=sess,
                                         feed_dict={x: test_data.all_images(),
                                                     y_: test_data.all_labels(),
                                                     keep_prob: 1.0})
        print("Calculated embeddings for inputs.", flush=True)
        print(train_embeddings.shape)
        return train_values[features_to_save], \
               test_values[features_to_save],\
               train_embeddings, test_embeddings
    features_to_save = ['SNPName', 'ProbeName', 'OverallZScore', 'expressionMean',
                        'expressionVar', 'TssDistance', 'methyMean', 'methyVar']
    train_eqtm, test_eqtm, train_embeddings, test_embeddings = \
        read_embeddings(sess, eqtm_path, features_to_save)

    # save embeddings to feature path
    def build_embeddings_dataframe(eqtm, embeddings):
        dataframe = pd.DataFrame(data=embeddings,
                                 index=eqtm.index,
                                 columns=['feature_%d' % ind for ind in range(embeddings.shape[1])])
        return pd.concat([dataframe, eqtm], axis=1)
    train_dataframe = build_embeddings_dataframe(train_eqtm, train_embeddings)
    print(train_dataframe.shape)
    print(train_dataframe.head())
    test_dataframe = build_embeddings_dataframe(test_eqtm, test_embeddings)

    save_dataframe = train_dataframe.append(test_dataframe)
    save_dataframe.to_csv(os.path.join(save_embeddings_dirpath, eqtm_name+'_embeddings.txt'))
    print("Saved results to path: ",
          os.path.join(save_embeddings_dirpath, eqtm_name+'_embeddings.txt'),
          flush=True)
