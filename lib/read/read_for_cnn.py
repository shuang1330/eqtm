import numpy as np
import os


def read_individual_image(img_path):
    # print(img_path)
    f = open(img_path, 'r').readlines()
    res = np.array([[int(element.strip()) for element in row.strip().split('\t')] for row in f])
    return res.reshape([res.shape[0], res.shape[1], 1])


class img(object):
    def __init__(self, imgname, img_dirpath):
        self._imgname = imgname
        self._dirpath = img_dirpath

    def imgname(self):
        return self._imgname

    def imgpath(self):
        return self._imgpath

    def img_data(self):
        self._imgpath = os.path.join(self._dirpath, self._imgname+'.txt')
        return read_individual_image(self._imgpath)


def find_correspondingImgs(item, cpgSites_dirpath):
    cpg_img = img(item[0], cpgSites_dirpath).img_data()
    # gene_img = img(item[1], geneSites_dirpath).img_data()
    # return np.stack((cpg_img, gene_img), axis=-1)
    return cpg_img


def dense_to_one_hot(labels_dense, num_classes):
    return np.eye(num_classes)[np.array(labels_dense).reshape(-1)]


class iterlist(object):
    def __init__(self, namearray, all_labels, img_dirpath):
        # namearray here is dataframe[['SNPName','ProbeName']]
        self._namearray = namearray
        self._all_labels = dense_to_one_hot(all_labels, 2)
        self._epochs_completed = 0
        self._index_in_epoch = 0
        self._num_examples = namearray.shape[0]
        self._img_dirpath = img_dirpath
        self._names = None
        self._data = None
        self._labels = None

    def all_labels(self):
        return self._all_labels

    def all_images(self):
        return np.array([find_correspondingImgs(
            item,
            self._img_dirpath
        ) for item in self._namearray])

    def next_batch(self, batch_size, shuffle=True):
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
            return_names = np.concatenate((names_rest_part, names_new_part),
                                          axis=0)
            return_labels = np.concatenate((labels_rest_part, labels_new_part),
                                           axis=0)
            return [np.array([find_correspondingImgs(
                item,
                self._img_dirpath
            ) for item in return_names]),
                    return_labels]
        else:
            self._index_in_epoch += batch_size
            end = self._index_in_epoch
            self._data = np.array([find_correspondingImgs(
                item,
                self._img_dirpath
            ) for item in self._names[start:end]])
            self._labels = self._all_labels[start:end]
            return [self._data, self._labels]


