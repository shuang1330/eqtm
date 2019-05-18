import __init__path
import os

from lib.model.net import large_kernel_cnn
from lib.read.read_data import read_eqtm_names
from lib.read.read_for_cnn import iterlist

import tensorflow as tf
import argparse


def train_model(cpgSites_dirpath,
                train_values, train_labels,
                test_values, test_labels,
                model_save_dirpath,
                batch_size=100,
                lr=1e-4,
                train_iter=20000,
                save_step=5000):
    # build the model
    x = tf.placeholder(tf.float32, shape=[None, 127, 32, 1])
    y_ = tf.placeholder(tf.float32, shape=[None, 2])
    y_conv, keep_prob, embedding = large_kernel_cnn(x)
    cross_entropy = tf.reduce_mean(
        tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y_conv))
    train_step = tf.train.AdamOptimizer(lr).minimize(cross_entropy)
    correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

    # start training
    train_data = iterlist(train_values, train_labels, cpgSites_dirpath)
    test_data = iterlist(test_values, test_labels, cpgSites_dirpath)
    batch_size = batch_size

    # record the results and save the models
    record_filepath = os.path.join(project_rootdir, 'commands', 'records',
                                   'test_cnn.txt')
    record_file = open(record_filepath, 'w')
    saver = tf.train.Saver()

    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        for i in range(train_iter+1):
            batch = train_data.next_batch(batch_size)
            train_step.run(feed_dict={x: batch[0],
                                      y_: batch[1],
                                      keep_prob: 0.5})
            if i % 100 == 0:
                train_accuracy = accuracy.eval(feed_dict={x: batch[0],
                                                          y_: batch[1],
                                                          keep_prob: 1.0})

                print('step %d, training accuracy %g' % (i, train_accuracy),
                      flush=True)
                record_file.write('step %d, training accuracy %g' % (i, train_accuracy))
                record_file.flush()
            if i > 0 and i % save_step == 0:
                save_model_path = os.path.join(model_save_dirpath,
                                               'test_CNN_trainIter_%d_lr%f_batchsize%d.ckpt' % (i, lr, batch_size))
                saver.save(sess, save_model_path)
                print('Model train in iteration %d saved in path %s' % (
                    i, save_model_path
                ),
                      flush=True)


        hidden_results = embedding.eval(feed_dict={x: test_data.all_images(),
                                                   y_: test_data.all_labels(),
                                                   keep_prob: 1.0})
        print(accuracy.eval(feed_dict={x: test_data.all_images(),
                                       y_: test_data.all_labels(),
                                       keep_prob: 1.0}),
              flush=True)
        record_file.close()
        return hidden_results


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project_rootdir",
        dest="project_rootdir",
        default="/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm"
    )
    parser.add_argument(
        "--eqtm_filename",
        dest="eqtm_filename",
        type=str,
        default="",
        help="eqtm file path"
    )
    parser.add_argument(
        "--batch_size",
        dest="batch_size",
        type=int,
        default=100,
        help="batch size for training"
    )
    parser.add_argument(
        "--lr",
        dest="lr",
        type=float,
        default=0.0001,
        help="learning rate for training"
    )
    parser.add_argument(
        "--save_prefix",
        dest="save_prefix",
        type=str,
        default="test_cnn_model",
        help="prefix for saving models"
    )
    parser.add_argument(
        "--save_step",
        dest="save_step",
        type=int,
        default=50000,
        help="prefix for saving models"
    )
    parser.add_argument(
        "--train_iter",
        dest="train_iter",
        type=int,
        default=20000,
        help="max training iterations"
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse()
    print("Execute with arguments: ", args, flush=True)
    project_rootdir = args.project_rootdir
    eqtm_path = os.path.join(project_rootdir, 'data', 'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             args.eqtm_filename + '.txt')
    # '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')

    # data_dirpath = os.path.join(project_rootdir,
    #                             'data', 'cpgSites',
    #                             'seperate_cpgFiles')
    cpgSites_dirpath = os.path.join(project_rootdir,
                                    'data', 'cpgSites',
                                    'seperate_cpgFiles')
    geneSites_dirpath = os.path.join(project_rootdir,
                                     'data', 'features',
                                     'geneOverlap',
                                     'seperate_geneFiles')
    model_save_dirpath = os.path.join(project_rootdir,
                                      'data',
                                      'model',
                                      args.save_prefix)
    # 'test_cnn')

    # load training and test data
    [train_values, train_labels, test_values, test_labels] = \
        read_eqtm_names(eqtm_path)

    _ = train_model(cpgSites_dirpath,
                    train_values[['SNPName', 'ProbeName']].values,
                    train_labels,
                    test_values[['SNPName', 'ProbeName']].values,
                    test_labels,
                    model_save_dirpath,
                    args.batch_size,
                    lr=args.lr,
                    train_iter=args.train_iter,
                    save_step=args.save_step)

    # for batch_size in [100, 200]:
    #     for lr in [0.01, 0.001, 0.1]:
    #         print('Using learning rate %f, batch_size %d' % (lr, batch_size), flush=True)
    #         _ = train_model(cpgSites_dirpath,
    #                                  train_values[['SNPName', 'ProbeName']].values,
    #                                  train_labels,
    #                                  test_values[['SNPName', 'ProbeName']].values,
    #                                  test_labels,
    #                                  model_save_dirpath,
    #                                  batch_size,
    #                                  lr=lr)
