
# TODO: maybe I should couple it with the dataset using
# instead of defining the node number?


import tensorflow as tf
import tensorflow.contrib.layers as lays
import tensorflow.contrib.slim as slim
from tensorflow.contrib.slim import losses
from tensorflow.contrib.slim import arg_scope

def autoencoder(x):
    with arg_scope([slim.fully_connected],
    weights_initializer=tf.truncated_normal_initializer(stddev=0.01),
    weights_regularizer=slim.l2_regularizer(0.5)):
        net1 = slim.fully_connected(x,50)
        # net1 = slim.fully_connected(net1,50)
        # net2 = slim.fully_connected(net1,100)
        net = slim.fully_connected(net1,214)
    return net, net1


def feedforward_net(x):
    with arg_scope([slim.fully_connected],
    weights_initializer=tf.truncated_normal_initializer(stddev=0.01),
    weights_regularizer=slim.l2_regularizer(0.005)):
        net = slim.fully_connected(x,10)
        # net = slim.fully_connected(net,4)
        net = slim.fully_connected(net,2)
        cls_prob = tf.nn.softmax(net, name="cls_prob")
    return net, cls_prob

def test_cnn(x):
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
    W_conv1 = weight_variable([5, 5, 1, 32])
    b_conv1 = bias_variable([32])
    x_image = tf.reshape(x, [-1, 28, 28, 2])
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool_2x2(h_conv1)
    # second layer
    W_conv2 = weight_variable([5, 5, 32, 64])
    b_conv2 = bias_variable([64])
    h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
    h_pool2 = max_pool_2x2(h_conv2)

    # fc layer
    W_fc1 = weight_variable([7 * 7 * 64, 1024])
    b_fc1 = bias_variable([1024])
    h_pool2_flat = tf.reshape(h_pool2, [-1, 7*7*64])
    h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

    # drop out layer
    keep_prob = tf.placeholder(tf.float32)
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # readout layer
    W_fc2 = weight_variable([1024, 10])
    b_fc2 = bias_variable([10])
    y_conv = tf.matmul(h_fc1_drop, W_fc2) + b_fc2

    return y_conv
