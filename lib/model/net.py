# TODO: maybe I should couple it with the dataset?
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


def test_cnn(x_image):
    def weight_variable(shape, name=None):
        initial = tf.truncated_normal(shape, stddev=0.1)
        return tf.Variable(initial, name=name)

    def bias_variable(shape, name=None):
        initial = tf.constant(0.1, shape=shape)
        return tf.Variable(initial, name=name)

    def conv2d(x, W):
        return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

    def max_pool_2x2(x):
        return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                              strides=[1, 2, 2, 1],
                              padding='SAME')

    # first layer
    W_conv1 = weight_variable([5, 5, 1, 32], name='w_conv1')
    b_conv1 = bias_variable([32],name='b_conv1')
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool_2x2(h_conv1)
    # second layer
    W_conv2 = weight_variable([5, 5, 32, 64], name='w_conv2')
    b_conv2 = bias_variable([64], name='b_conv2')
    h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
    h_pool2 = max_pool_2x2(h_conv2)

    # second layer
    W_conv3 = weight_variable([5, 5, 64, 32], name='w_conv2')
    b_conv3 = bias_variable([32], name='b_conv2')
    h_conv3 = tf.nn.relu(conv2d(h_pool2, W_conv3) + b_conv3)
    h_pool3 = max_pool_2x2(h_conv3)

    # fc layer
    W_fc1 = weight_variable([2048, 1024], name='w_fc1') #8192
    b_fc1 = bias_variable([1024], name='b_fc1')
    h_pool2_flat = tf.reshape(h_pool3, [-1, 2048])
    h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

    # drop out layer
    keep_prob = tf.placeholder(tf.float32)
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # fc layer2
    W_fc2 = weight_variable([1024, 32], name='w_fc2')
    b_fc2 = bias_variable([32], name='b_fc2')
    h_fc2 = tf.matmul(h_fc1_drop, W_fc2) + b_fc2

    # drop out layer
    h_fc2_drop = tf.nn.dropout(h_fc2, keep_prob)

    # readout layer
    W_fc3 = weight_variable([32, 2], name='w_fc3')
    b_fc3 = bias_variable([2], name='b_fc3')
    y_conv = tf.matmul(h_fc2_drop, W_fc3)+b_fc3

    return y_conv, keep_prob, h_fc2


def large_kernel_cnn(x_image):
    def weight_variable(shape, name=None):
        initial = tf.truncated_normal(shape, stddev=0.1)
        return tf.Variable(initial, name=name)

    def bias_variable(shape, name=None):
        initial = tf.constant(0.1, shape=shape)
        return tf.Variable(initial, name=name)

    def conv2d(x, W):
        return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

    def max_pool_2x2(x):
        return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                              strides=[1, 2, 2, 1],
                              padding='SAME')

    # first layer
    W_conv1 = weight_variable([5, 32, 1, 32], name='w_conv1')
    b_conv1 = bias_variable([32],name='b_conv1')
    h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
    h_pool1 = max_pool_2x2(h_conv1)
    # second layer
    W_conv2 = weight_variable([5, 8, 32, 64], name='w_conv2')
    b_conv2 = bias_variable([64], name='b_conv2')
    h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2) + b_conv2)
    h_pool2 = max_pool_2x2(h_conv2)

    # second layer
    W_conv3 = weight_variable([2, 4, 64, 32], name='w_conv2')
    b_conv3 = bias_variable([32], name='b_conv2')
    h_conv3 = tf.nn.relu(conv2d(h_pool2, W_conv3) + b_conv3)
    h_pool3 = max_pool_2x2(h_conv3)

    # fc layer
    W_fc1 = weight_variable([2048, 1024], name='w_fc1') #8192
    b_fc1 = bias_variable([1024], name='b_fc1')
    h_pool2_flat = tf.reshape(h_pool3, [-1, 2048])
    h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

    # drop out layer
    keep_prob = tf.placeholder(tf.float32)
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    # fc layer2
    W_fc2 = weight_variable([1024, 32], name='w_fc2')
    b_fc2 = bias_variable([32], name='b_fc2')
    h_fc2 = tf.matmul(h_fc1_drop, W_fc2) + b_fc2

    # drop out layer
    h_fc2_drop = tf.nn.dropout(h_fc2, keep_prob)

    # readout layer
    W_fc3 = weight_variable([32, 2], name='w_fc3')
    b_fc3 = bias_variable([2], name='b_fc3')
    y_conv = tf.matmul(h_fc2_drop, W_fc3)+b_fc3

    return y_conv, keep_prob, h_fc2


def conv_1d(x_vector):
    """

    :param x_vector: 32-length vector
    :return:
    """
    # input = tf.reshape(x_vector, [None, 32, 1])
    net = tf.layers.conv1d(x_vector, 8, 3, name="conv1d_1")
    print(net.shape)
    net = tf.layers.max_pooling1d(net, 2, 2)
    print(net.shape)
    net = tf.layers.conv1d(net, 8, 3, name="conv1d_2")
    print(net.shape)
    net = tf.layers.max_pooling1d(net, 4, 2)
    print(net.shape)
    net = tf.reshape(net, [-1, 40])
    net = tf.layers.dense(net, 16)
    net = tf.layers.dense(net, 2)
    print(net.shape)
    return net