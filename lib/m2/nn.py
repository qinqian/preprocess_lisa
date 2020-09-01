"""
Module for Traditional Neural Networks
"""

from __future__ import absolute_import, print_function, division
from theano import pp
from theano.misc.pkl_utils import dump
from theano.misc.pkl_utils import load
import six.moves.cPickle as pickle
import os
from sklearn.cross_validation import KFold
import timeit
import numpy as np
import sys
import theano.tensor as T
import theano
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from theano.sandbox.rng_mrg import MRG_RandomStreams
from theano.tensor.shared_randomstreams import RandomStreams

def get_pickling_errors(obj,seen=None):
    """
    try to see which property of the object cannot be pickled
    :param obj: object
    :param seen:
    :return:
    """
    if seen == None:
        seen = []
    try:
        state = obj.__getstate__()
    except AttributeError:
        return
    if state == None:
        return
    if isinstance(state,tuple):
        if not isinstance(state[0],dict):
            state=state[1]
        else:
            state=state[0].update(state[1])
    result = {}
    for i in state:
        try:
            pickle.dumps(state[i],protocol=2)
        except pickle.PicklingError:
            if not state[i] in seen:
                seen.append(state[i])
                result[i]=get_pickling_errors(state[i],seen)
    return result


class LogisticRegression(object):
    def __init__(self, rng, input, n_in, n_out, param=None, activation = T.nnet.softmax):
        # Randomly generate weight matrix W
        W_values = np.asarray(
            rng.uniform(
                low=-np.sqrt(6. / (n_in + n_out)),
                high=np.sqrt(6. / (n_in + n_out)),
                size=(n_in, n_out)
            ),
            dtype=theano.config.floatX
        )
        # W_values = n_in**(-0.5)*np.asarray(
        #     rng.randn(
        #         n_in, n_out
        #     ),
        #     dtype=theano.config.floatX
        # )
        # W_values = 10**((np.log10(1e-2) - np.log10(1e-5))*np.asarray(
        #     rng.uniform(
        #         low=1e-5,
        #         high=1e-2,
        #         size=(n_in, n_out)
        #     ),
        #     dtype=theano.config.floatX
        # )+np.log10(1e-5))

        if activation == T.nnet.sigmoid:
            W_values *= 4
        self.activation = activation
        self.input = input

        # Randomly initialize bias vector b
        b_values = np.zeros((n_out,), dtype=theano.config.floatX)

        self.W = theano.shared(value=W_values, name='W', borrow=True)
        self.b = theano.shared(value=b_values, name='b', borrow=True)

        # Define feed-forward function
        lin_output = T.dot(self.input, self.W) + self.b
        self.p_y_given_x = (
                lin_output if activation is None
                else self.activation(lin_output)
        )

        # Specify parameters of layer
        self.params = [self.W, self.b]

        # prediction of y labels
        self.y_pred = T.argmax(self.p_y_given_x, axis=1)
        # theano.printing.pprint(self.y_pred)
        # theano.printing.pydotprint(self.y_pred, outfile="my_nn.png", var_with_name_simple=True)

    def negative_log_likelihood(self, y):
        # y corresponds to a vector that gives for each data sample
        # the correct label
        return -T.mean(T.log(self.p_y_given_x)[T.arange(y.shape[0]), y])
        # return T.nnet.binary_crossentropy(self.p_y_given_x, y).mean()

    def predict(self, y):
        return self.p_y_given_x[T.arange(y.shape[0]), 1] ## use the 1st column's probability

    def dropout_prediction(self):
        self.p_y_given_x = self.activation(T.dot(self.input, self.W) + self.b)
        self.y_pred = T.argmax(self.p_y_given_x, axis=1)

    def errors(self, y):
        if y.ndim != self.y_pred.ndim:
            raise TypeError(
                ' y should have the same shape as y_pred '
            )
        if y.dtype.startswith('int'):
            return T.mean(T.neq(self.y_pred, y))
        else:
            raise NotImplementedError()


class HiddenLayer(object):
    def __init__(self, rng, input, n_in, n_out, W=None, b=None, activation=T.tanh):

        self.input = input

        # Randomly generate weight matrix W
        if W is None:
            W_values = np.asarray(
                rng.uniform(
                    low=-np.sqrt(6. / (n_in + n_out)),
                    high=np.sqrt(6. / (n_in + n_out)),
                    size=(n_in, n_out)
                ),
                dtype=theano.config.floatX
            )
    
            if activation == T.nnet.sigmoid:
                W_values *= 4
    
            self.W = theano.shared(value=W_values, name='W', borrow=True)        

        if b is None:
            # Randomly initialize bias vector b
            b_values = np.zeros((n_out,), dtype=theano.config.floatX)
            self.b = theano.shared(value=b_values, name='b', borrow=True)
        
        # Define feed-forward function
        lin_output = T.dot(input, self.W) + self.b
        self.output = (
            lin_output if activation is None
            else activation(lin_output)
        )

        # Specify parameters of layer
        self.params = [self.W, self.b]
        
    def feedforward(self, input):
        # Function to apply layer feedforward to new input data
        return T.tanh(T.dot(input, self.W) + self.b)


class FeatureMapHiddenLayer(object):
    """ Hidden-Layer class with customized feature mapping
     for a certain number of groups each with several features(samples),
     it map features from the same group to certain number of
     hidden units. The motivation is to add nonlinear property
     to the same feature source, e.g. 10 samples from H3K27ac,
     the class present certain number of feature maps to next
     layer
    """
    def __init__(self, rng, input, n_group, n_in, n_out, seed, dropout, pre_training=None, activation=T.tanh):
        """ Initialize the parameters for the customized hidden layer

        :type rng: numpy.random.RandomState
        :param rng: a random number generator used to initialize weights and biases

        :param inputs: a T.tensor, whose rows is examples of a mini-batch,
                                   whos cols is feature_num = n_group * n_in(features of each group)

        :param n_group: the number of feature groups, such as 4 histone mark types
        :param n_in: a list of ints, for the number of features in one mark group
        :param n_out: the number of output for a feature group, that is a feature map with n_out hidden units,
                      such as map 10 features(samples) of h3k4me3 to 2 hidden units with nonlinear transformation

        :type pre_training: numpy array
        :param pre_training: use logistic regression weights to replace one column of W (from input to hidden layer)

        """
        self.input = input
        msrng = MRG_RandomStreams(seed=seed)

        # dropout
        self.dropout = T.switch(msrng.binomial(size=self.input.shape, p=dropout), self.input, 0)

        self.activation = activation
        self.Ws = []
        self.bs = []
        self.feature_maps = []
        self.n_group = n_group
        self.n_in = n_in
        self.n_out = n_out

        feature_end = 0
        print(n_group)
        print(self.n_in)
        print(self.n_out)
        for i in np.arange(n_group):
            W_values = np.asarray(
                    rng.uniform(
                            low=-np.sqrt(6. / (self.n_in[i] + self.n_out)),
                            high=np.sqrt(6. / (self.n_in[i] + self.n_out)),
                            size=(self.n_in[i], self.n_out)
                    ),
                    dtype=theano.config.floatX
            )
            # W_values = self.n_in[i]**(-0.5)*np.asarray(
            #     rng.randn(
            #         self.n_in[i], self.n_out
            #     ),
            #     dtype=theano.config.floatX
            # )
            # r = rng.uniform(
            #         low=1e-7,
            #         high=1e-3,
            #         size=(n_in[i], n_out)
            #     )
            # W_values = 10**((np.log10(1e-3) - np.log10(1e-7))*np.asarray(r, dtype=theano.config.floatX)+np.log10(1e-7))
            
            if self.activation == T.nnet.sigmoid:
                W_values *= 4

            if i == 0:
                feature_start = 0
                feature_end = self.n_in[i]
            else:
                feature_start = feature_end
                feature_end = feature_start + self.n_in[i]

            # obsoleted, shutdonw insertion of logistic regression coefficients
            # if pre_training is not None:
            #     # inject into one position of W weight matrix
            #     inject_W_index = rng.choice(n_out, 1, False)[0]
            #     W_values[:,inject_W_index] = pre_training[:, feature_start:feature_end]

            W = theano.shared(value=W_values, name='W', borrow=True)
            self.Ws.append(W)

            # Randomly initialize bias vector b
            b_values = np.zeros((n_out,), dtype=theano.config.floatX)
            b = theano.shared(value=b_values, name='b', borrow=True)
            self.bs.append(b)

            lin_output = T.dot(self.dropout[:, T.arange(feature_start, feature_end)], W) + b

            self.feature_maps.append(
                lin_output if self.activation is None
                else self.activation(lin_output)
            )

        self.output = T.concatenate(self.feature_maps, axis=1)
        self.params = self.Ws + self.bs


    def feedforward(self):
        """
        Function to apply layer feedforward to new input data
        :return:
        """
        feature_maps = []
        feature_end = 0

        for i in np.arange(self.n_group):
            # half the input layer in testing stage
            if i == 0:
                feature_start = 0
                feature_end = self.n_in[i]
            else:
                feature_start = feature_end
                feature_end = feature_start + self.n_in[i]

            lin_output = T.dot(0.5*self.input[:, T.arange(feature_start, feature_end)], self.Ws[i]) + self.bs[i]

            feature_maps.append(
                lin_output if self.activation is None
                else self.activation(lin_output)
            )

        return T.concatenate(feature_maps, axis=1)


class MyMLP(object):
    """Customized Multi-Layer Neural network class
    hidden_layer number of layers neural network
    """
    def __init__(self, rng, input, n_group, n_in, n_hidden, dropout_in, dropout_h, activation, n_out, seed, pre_training=None):
        """Initialize the parameters for the multilayer perceptron

        :type rng: numpy.random.RandomState
        :param rng: a random number generator used to initialize weights

        :type input: theano.tensor.TensorType
        :param input: symbolic variable that describes the input of the
        architecture (one minibatch)

        :type n_group: int
        :param n_group: n groups of features, such as different histone mark number

        :type n_in: a list of int
        :param n_in: number of input units for each of the group, it has the same len as n_group

        :type n_hidden: int
        :param n_hidden: number of hidden units

        :type n_out: int
        :param n_out: number of output units, the dimension of the space in which the labels lie

        :type pre_training: numpy array
        :param pre_training: use logistic regression weights to replace one column of W (from input to hidden layer)
        """
        self.input = input
        srng = MRG_RandomStreams(seed=seed)        
        if n_hidden >= 1:
            self.hiddenLayer = FeatureMapHiddenLayer(
                    rng=rng,
                    input=self.input,
                    seed=seed,
                    n_group=n_group,
                    n_in=n_in,
                    n_out=n_hidden,
                    dropout=dropout_in, 
                    pre_training=pre_training,
                    activation=activation,
            )

            self.dropout_output = T.switch(srng.binomial(size=self.hiddenLayer.output.shape, p=dropout_h),
                                           self.hiddenLayer.output, 0)
        else:
            self.dropout_output = T.switch(srng.binomial(size=self.input.shape, p=dropout_in),
                                           self.input, 0)


        if n_hidden >= 1:
            # The logistic regression layer gets as input the hidden units
            # of the hidden layer
            self.logRegressionLayer = LogisticRegression(
                    rng=rng,
                    input=self.dropout_output,
                    n_in=n_hidden * n_group,
                    n_out=n_out,
                    activation=T.nnet.softmax
            )
            self.params = self.hiddenLayer.params + self.logRegressionLayer.params
            self.L1 = T.sum(abs(T.concatenate(self.hiddenLayer.Ws, axis=1)), axis=None) + abs(self.logRegressionLayer.W).sum()
            self.L2_sqr = T.sqr(T.concatenate(self.hiddenLayer.Ws, axis=1)).sum() + (T.sqr(self.logRegressionLayer.W)).sum()
        else:
            # The logistic regression layer gets as input the hidden units
            # of the hidden layer
            self.logRegressionLayer = LogisticRegression(
                    rng=rng,
                    input=self.dropout_output,
                    n_in=sum(n_in),
                    n_out=n_out,
                    activation=T.nnet.softmax
            )
            self.params = self.logRegressionLayer.params
            self.L1 = abs(self.logRegressionLayer.W).sum()
            self.L2_sqr = (T.sqr(self.logRegressionLayer.W)).sum()


        # only regularize the hidden layer
        # self.L1 = T.sum(abs(T.concatenate(self.hiddenLayer.Ws, axis=1)), axis=None)
        # self.L2_sqr = T.sqr(T.concatenate(self.hiddenLayer.Ws, axis=1)).sum()
        # negative log likelihood of the MLP is given by the negative
        # log likelihood of the output of the model, computed in the
        # logistic regression layer
        self.negative_log_likelihood = (
            self.logRegressionLayer.negative_log_likelihood
        )

        if n_hidden >= 1:
            dropout_prediction = 0.5 * self.hiddenLayer.feedforward()
        else:
            dropout_prediction = 0.5 * self.input
        self.logRegressionLayer.input = dropout_prediction
        self.logRegressionLayer.dropout_prediction()
        # same holds for the function computing the number of errors
        self.errors = self.logRegressionLayer.errors
        # holds for the function predicting the probability score
        self.predict = self.logRegressionLayer.predict
    # def dropout(self):
    #     """ drop out input layer and hidden layer by random ratio
    #     """
    #     self.logRegressionLayer.input = 0.5 * self.hiddenLayer.feedforward()
    #     self.logRegressionLayer.dropout()
    #     self.errors = self.logRegressionLayer.errors
    #     self.predict = self.logRegressionLayer.predict
    def __getstate__(self):
        """ remove memebers which cannot be pickled
        """
        state = dict(self.__dict__)
        # remove unwanted property
        del state['predict']
        del state['errors']
        del state['negative_log_likelihood']
        return state

    def __setstate__(self, state):
        """ assign back the instancemethod which cannot be pickled
        """
        # somehow this would cause the issue
        # donot asssign input or cannot find input values for the graph....!!!
        self.__dict__.update(state)
        self.predict = self.logRegressionLayer.predict

class MLP(object):
    """Multi-Layer Perceptron Class

    A multilayer perceptron is a feedforward artificial neural network model
    that has one layer or more of hidden units and nonlinear activations.
    Intermediate layers usually have as activation function tanh or the
    sigmoid function (defined here by a ``HiddenLayer`` class)  while the
    top layer is a softmax layer (defined here by a ``LogisticRegression``
    class).
    """
    def __init__(self, rng, input, n_in, n_hidden, n_out):
        """Initialize the parameters for the multilayer perceptron

        :type rng: numpy.random.RandomState
        :param rng: a random number generator used to initialize weights

        :type input: theano.tensor.TensorType
        :param input: symbolic variable that describes the input of the
        architecture (one minibatch)

        :type n_in: int
        :param n_in: number of input units, the dimension of the space in
        which the datapoints lie

        :type n_hidden: int
        :param n_hidden: number of hidden units

        :type n_out: int
        :param n_out: number of output units, the dimension of the space in
        which the labels lie

        """

        # Since we are dealing with a one hidden layer MLP, this will translate
        # into a HiddenLayer with a tanh activation function connected to the
        # LogisticRegression layer; the activation function can be replaced by
        # sigmoid or any other nonlinear function
        self.hiddenLayer = HiddenLayer(
                rng=rng,
                input=input,
                n_in=n_in,
                n_out=n_hidden,
                activation=T.tanh
        )

        # The logistic regression layer gets as input the hidden units
        # of the hidden layer
        self.logRegressionLayer = LogisticRegression(
                rng=rng,
                input=self.hiddenLayer.output,
                n_in=n_hidden,
                n_out=n_out,
                activation=T.nnet.sigmoid
        )

        self.L1 = (
            abs(self.hiddenLayer.W).sum()
            + abs(self.logRegressionLayer.W).sum()
        )
        self.L2_sqr = (
            (self.hiddenLayer.W**2).sum()
            +(self.logRegressionLayer.W**2).sum()
        )

        # negative log likelihood of the MLP is given by the negative
        # log likelihood of the output of the model, computed in the
        # logistic regression layer
        self.negative_log_likelihood = (
            self.logRegressionLayer.negative_log_likelihood
        )
        # same holds for the function computing the number of errors
        self.errors = self.logRegressionLayer.errors
        # holds for the function predicting the probability score
        self.predict = self.logRegressionLayer.predict
        self.params = self.hiddenLayer.params + self.logRegressionLayer.params

    def __getstate__(self):
        """ remove memebers which cannot be pickled
        """
        state = dict(self.__dict__)
        # remove unwanted property
        del state['predict']
        del state['errors']
        del state['negative_log_likelihood']
        return state

    def __setstate__(self, state):
        """ assign back the instancemethod which cannot be pickled
        """
        # somehow this would cause the issue
        # donot asssign input or cannot find input values for the graph....!!!
        self.__dict__.update(state)
        self.predict = self.logRegressionLayer.predict


def RMSprop(cost, params, lr=0.001, rho=0.9, epsilon=1e-6):
    grads = T.grad(cost=cost, wrt=params)
    updates = []
    for p, g in zip(params, grads):
        acc = theano.shared(p.get_value(borrow=True) * 0.)
        acc_new = rho * acc + (1 - rho) * g ** 2
        gradient_scaling = T.sqrt(acc_new + epsilon)
        g = g / gradient_scaling
        updates.append((acc, acc_new))
        updates.append((p, p - lr * g))
    return updates

def Nesterov(cost, params, lr=0.001, mu=0.9):
    updates = []
    for p, g in zip(params, T.grad(cost=cost, wrt=params)):
        d = theano.shared(p.get_value() * 0., name="momentum", borrow=True)
        updates.append((d, mu*d - lr *g))
        updates.append((p, p + mu*mu*d - (1+mu) * lr * g))
    return updates


class SMOTE:
    """ https://www.jair.org/media/953/live-953-2037-jair.pdf
    synthetic minority class samples by using insertion values of the k nearest samples
    """
    def __init__(self, nnarray, N, x, y, seed):
        """
        T: minority class sample index
        N: the number of synthetic samples
        K: K nearest samples
        """
        self.nnarray = nnarray
        self.N = N
        self.x = x
        self.y = y
        self.seed = seed

    def populate(self):
        """
        return train model
        """
        srng = RandomStreams(seed=self.seed)
        rng = MRG_RandomStreams(seed=self.seed)

        added_train_index = srng.choice((self.N, ), T.arange(self.x.shape[0]), True)
        # copy state for synthetic x and y
        synthetic_x = self.x[added_train_index]
        state_after_train_x_index = added_train_index.rng.get_value(borrow=True).get_state()
        srng2 = added_train_index.rng.get_value(borrow=True)
        srng2.set_state(state_after_train_x_index)
        added_train_index.rng.set_value(srng2, borrow=True)
        synthetic_y = self.y[added_train_index]

        srng3 = added_train_index.rng.get_value(borrow=True)
        srng3.set_state(state_after_train_x_index)
        added_train_index.rng.set_value(srng3, borrow=True)

        # gaps size equal to feature number
        gaps = rng.uniform(low=0, high=1, size=(self.N, synthetic_x.shape[1]), dtype=T.config.floatX)
        # choose one randomly from the k nearest samples, not yet
        nn = srng.choice((self.N, ), self.nnarray.shape[0], True)
        synthetic_x += gaps * (self.x[self.nnarray[nn,added_train_index]]-synthetic_x)
        return synthetic_x, synthetic_y


def distance():
    """
    squared_euclidean_distances
    """
    X = T.matrix('X')
    Y = T.matrix('Y')
    output = T.sqrt(abs((X ** 2).sum(1).reshape((X.shape[0], 1)) + (Y ** 2).sum(1).reshape((1, Y.shape[0])) - 2 * X.dot(Y.T)))
    return theano.function([X, Y], output, allow_input_downcast=True)


def expand_positive(seed, nnarray, train_pos_x, train_pos_y, batch_size, index, train_set_x, train_set_y, x, y):
    """ different methods for expanding the positive examples for a mini-batch
    default method is SMOTE for over-sampling

    train_batch_x: a mini-batch of the original training x (train_set_x)
    train_batch_y: corresponding mini-batch of the original training y (train_set_y)

    return: tensor symbolic variable of the input train x and y
    """
    # oversampling with repeated positive samples
    # would cause serious overfitting
    train_batch_x = train_set_x[index*batch_size:(index+1)*batch_size]
    train_batch_y = train_set_y[index*batch_size:(index+1)*batch_size]

    def expand(train_batch_x, train_batch_y):
        s = SMOTE(nnarray, abs(batch_size*4//5-train_batch_y.sum()), train_pos_x, train_pos_y, seed)
        add_train_x, add_train_y = s.populate()
        return T.concatenate([train_batch_x, add_train_x], axis=0), T.concatenate([train_batch_y, add_train_y], axis=0)

    from theano.ifelse import ifelse
    train_x, train_y = ifelse(T.lt(1.0*train_batch_y.sum()//batch_size, 0.5),
                              expand(train_batch_x, train_batch_y),
                              (train_batch_x, train_batch_y))
    return [(x, train_x), (y, train_y)]


def shared_data_fold(seed, rng, X_test, Y_test, X_train, Y_train, Kfold):
    """
    :rtype theano shared variable
    :param seed: global random seed
    :param x: predictor, such as H3K4me3 regulatory potential pandas DataFrame
    :param y: response, such as E2 stimuli differential gene list pandas Series
    :return: shared variable for train_set_x, valid_set_x, test_set_x, train_set_y, valid_set_y, test_set_y
    """
    # K-fold as l1 logistic regression cross validation
    ra = rng.permutation(X_train.shape[0])
    X_train = X_train[ra]
    Y_train = Y_train[ra]

    kf = KFold(X_train.shape[0], n_folds=Kfold, random_state=seed)
    for train_index, valid_index in kf:
        X_train2,  X_valid = X_train[train_index], X_train[valid_index]
        Y_train2,  Y_valid = Y_train[train_index], Y_train[valid_index]

        train_set_x = theano.shared(value=X_train2,
                                    name="train_set_x", borrow=True)
        valid_set_x = theano.shared(value=X_valid,
                                    name="valid_set_x", borrow=True)
        test_set_x = theano.shared(value=X_test,
                                   name="test_set_x", borrow=True)

        train_set_y = theano.shared(value=Y_train2, name="train_set_y", borrow=True)
        valid_set_y = theano.shared(value=Y_valid, name="valid_set_y", borrow=True)
        test_set_y = theano.shared(value=Y_test, name="test_set_y", borrow=True)
        yield train_set_x, valid_set_x, test_set_x, train_set_y, valid_set_y, test_set_y


def shared_data(seed, X_test, Y_test, X_train, Y_train):
    rng = np.random.RandomState(seed)

    ra = rng.permutation(X_train.shape[0])

    X_train = X_train[ra]
    Y_train = Y_train[ra]

    train_set_x = theano.shared(value=X_train,
                                name="train_set_x", borrow=True)
    test_set_x = theano.shared(value=X_test,
                               name="test_set_x", borrow=True)
    train_set_y = theano.shared(value=Y_train, name="train_set_y", borrow=True)
    test_set_y = theano.shared(value=Y_test, name="test_set_y", borrow=True)
    return train_set_x, test_set_x, train_set_y, test_set_y

def mtanh(x):
    """lecun's symmetric sigmoid"""
    return 1.7159*T.tanh((2.0/3.0)*x)

def generate(seed):
    np.random.seed(seed)
    uniform = lambda low, high: np.asarray(np.random.uniform(low=low, high=high, size=1), dtype='float32')

    # neural network
    hidden_activation = T.tanh # np.random.choice((T.tanh, T.nnet.softplus), 1)[0]
    n_hidden = np.random.choice([0, 1, 2, 4, 8, 16], 1)[0]
    batch_size = np.random.choice([20, 32, 64, 128], 1)[0]
    n_epochs = 500 # np.random.choice([500, 1000, 2000], 1)[0]
    # gradient descent
    rho = ((0.99-0.90)*np.sqrt(uniform(0.90, 0.99))+0.90)[0]
    learning_rate = (10**((np.log10(0.05)-np.log10(0.0005))*uniform(0.0005, 0.05)+np.log10(0.0005)))[0]
    epsilon = ((1e-6-1e-8)*np.sqrt(uniform(1e-8, 1e-6))+1e-8)[0]

    # regulatorization
    # L2_reg = (10**((np.log10(1e-3)-np.log10(1e-15))*uniform(1e-15, 1e-3)+np.log10(1e-15)))[0]
    # L1_reg = (10**((np.log10(1e-3)-np.log10(1e-15))*uniform(1e-15, 1e-3)+np.log10(1e-15)))[0]
    L2_reg = (10**((np.log10(1e-3)-np.log10(1e-5))*uniform(1e-5, 1e-3)+np.log10(1e-5)))[0]
    L1_reg = (10**((np.log10(1e-3)-np.log10(1e-5))*uniform(1e-5, 1e-3)+np.log10(1e-5)))[0]

    dropout_in = np.random.choice((0.5, 0.7, 0.8, 0.95), size=1)[0]
    dropout_h = 0.5 # np.random.choice((0.5, 0.7, 0.8, 0.95, 1.00), size=1)[0]

    ## early stopping
    patience = 1000 # np.random.choice([5000, 10000, 15000], 1)[0]
    patience_increase = 2
    improvement_threshold = 0.995  # a relative improvement of this much is considered significant
    return locals()


def run_mlp( seed,
             X_test, Y_test, X_train, Y_train, 
             features_num, repeat_times,
             diff_label,
             sparse_connect=False,
             pre_training=None,
             ## this is the factor type number, available with sparse_connect=True
             sparse_group = None,
             ):
    """
     use deep bind mode of traing
     different from theano tutorial
     first use cross valiation to choose the best paramters set
     according to maximum average validation AUCs
     then train on all traning genes with the best paramters set 3 times
     test with the best trained models in terms of the best trained AUCs
     to see if it is better than theano tutorial
     since generalization of theano tutorial way is
     largely determined by validation set choice

     Mini batch stochastic gradient descent optimization for a multilayer perceptron
    :type data: tuples of pandas DataFrame
    :param t: train_set_x, valid_set_x, test_set_x, train_set_y, valid_set_y, test_set_y

    :type pre_training: numpy array
    :param pre_training: use logistic L1 regression weight to replace one of the weights mapping input layer to hidden layer

    :type sparse_connect: bool
    :param sparse_connect: sparse connection or not

    :type sparse_group: int
    :param sparse_group: group of factor number

    :type t: int
    :param t: learning repeat times

    :type learning_rate: float
    :param learning_rate: learning rate used for the stochastic gradient

    :type L1_reg: float
    :param L1_reg: L1-norm's weight when added to the cost (see
    regularization)

    :type L2_reg: float
    :param L2_reg: L2-norm's weight when added to the cost (see
    regularization)

    :type n_epochs: int
    :param n_epochs: maximal number of epochs to run the optimizer

    :type patience: int
    :param patience: early stopping tweak parameter
    """

    maxsamples = 0
    if features_num:
        maxsamples = features_num[0] / repeat_times

    if hasattr(X_test, 'values'):
        X_train, X_test = np.asarray(X_train.values, dtype=theano.config.floatX), \
                          np.asarray(X_test.values, dtype=theano.config.floatX)
        Y_train, Y_test = np.asarray(Y_train.values, dtype='int32'), \
                          np.asarray(Y_test.values, dtype='int32')
    else:
        X_train, X_test = np.asarray(X_train.values, dtype=theano.config.floatX), \
                          np.asarray(X_test, dtype=theano.config.floatX)
        Y_train, Y_test = np.asarray(Y_train.values, dtype='int32'), \
                          np.asarray(Y_test, dtype='int32')


    print(X_test.shape, Y_test.shape, X_train.shape, Y_train.shape)
    # prepare shared data for neural network
    model_index = 0
    models = []; models_performance = []
    models_calibration = []

    parameters_set = 2
    for ps in xrange(parameters_set):
        seed = seed + ps
        params = generate(seed)
        del params['uniform']
        print(params)
        rng = np.random.RandomState(seed)
        Kfold = 3
        validation_aucs = np.zeros((Kfold, ))
        validation_prs = np.zeros((Kfold, ))
        # work through each of the K folds
        for fold, data in enumerate(list(shared_data_fold(seed, rng, X_test, Y_test, X_train, Y_train, Kfold=Kfold))): 
            train_set_x, valid_set_x, test_set_x, train_set_y, valid_set_y, test_set_y = data
            # to balance within the mini-batch
            # each mini-batch sampling from train_pos_index
            # to get equal number of pos and neg label
            # preprocessing
            train_pos_index = np.asarray(np.where(train_set_y.get_value(borrow=True)==1)[0], dtype='int32')
            train_pos_x = train_set_x.get_value(borrow=True)[train_pos_index]
            train_pos_y = train_set_y.get_value(borrow=True)[train_pos_index]
            dist = distance()
            dist = dist(train_pos_x, train_pos_x)

            # 3 nearest neighbor for SMOTE
            nnarray = np.argsort(dist, axis=1)[:,1:3].transpose()
            nnarray = theano.shared(value = nnarray, name = 'k_nearest', borrow=True)
            train_pos_x = theano.shared(value=train_pos_x, name = "train_pos_x", borrow=True)
            train_pos_y = theano.shared(value=train_pos_y, name = "train_pos_y", borrow=True)

            # compute number of minibatches for training, validation and testing
            n_train_batches = train_set_x.get_value(borrow=True).shape[0] // params['batch_size']
            # n_valid_batches = valid_set_x.get_value(borrow=True).shape[0] // params['batch_size']
            # n_test_batches = test_set_x.get_value(borrow=True).shape[0] // params['batch_size']

            print('... building the model')

            index = T.lscalar()  # index to a [mini]batch
            # allocate symbolic variables for the data
            x = T.matrix('x')
            y = T.ivector('y')  # the labels are presented as 1D vector of [int] labels

            if not sparse_connect:
                # construct the MLP class of full connection NN
                classifier = MLP(
                        rng=rng,
                        input=x,
                        n_in=train_set_x.get_value(borrow=True).shape[1],
                        n_hidden=params['n_hidden'],
                        n_out=2
                )

            else:
                # construct the MyMLP class of sparse connection NN
                classifier = MyMLP(
                        rng=rng,
                        input=x,
                        seed=seed,
                        n_group=sparse_group,
                        n_in=features_num,
                        n_hidden=params['n_hidden'],
                        dropout_in=params['dropout_in'],
                        dropout_h=params['dropout_h'],
                        activation=params['hidden_activation'],
                        n_out=2,
                        pre_training=pre_training,
                )

            # the cost we minimize during training is the negative log likelihood of
            # the model plus the regularization terms (L1 and L2); cost is expressed
            # here symbolically
            cost = (
                classifier.negative_log_likelihood(y)
                + params['L1_reg'] * classifier.L1
                + params['L2_reg'] * classifier.L2_sqr
            )

            # specify how to update the parameters of the model as a list of
            # (variable, update expression) pairs
            # gradient descent
            # gparams = [T.grad(cost, param) for param in classifier.params]
            # updates = [
            #     (param, param - learning_rate * gparam)
            #     for param, gparam in zip(classifier.params, gparams)
            #     ]

            # momentum gradient descent
            #momentum = 0.9
            #updates = []
            #for p, g in zip(classifier.params, gparams):
            #   mparam_i = theano.shared(np.zeros(p.get_value().shape, dtype=theano.config.floatX), name="momentum", borrow=True)
            #   v = momentum * mparam_i - learning_rate * g
            #   updates.append((mparam_i, v))
            #   updates.append((p, p + v))

            # Nesterov momentum
            # updates = Nesterov(cost, classifier.params, lr=params['learning_rate'], mu=params['rho'])
            updates = RMSprop(cost, classifier.params, lr=params['learning_rate'], rho=params['rho'], epsilon=params['epsilon'])
            # compiling a Theano function `train_model` that returns the cost, but
            # in the same time updates the parameter of the model based on the rules
            # defined in `updates`
            # training with ordered training examples
            # we need to balance positive and negative class proportion within a mini-batch
            # sampling the same number as batch size of positive example from training gene
            # to inject into the mini-batch

            # srng = RandomStreams(seed=seed)
            # added_train_index = srng.choice((batch_size*9//10, ), train_pos_index, False)
            # state_after_train_x_index = added_train_index.rng.get_value(borrow=True).get_state()
            # added_train_x = train_set_x[added_train_index]

            # rng = added_train_index.rng.get_value(borrow=True)
            # rng.set_state(state_after_train_x_index)
            # added_train_index.rng.set_value(rng, borrow=True)
            # added_train_y = train_set_y[added_train_index]

            batch_size = params['batch_size']

            # train_model = theano.function(
            #         inputs=[index],
            #         outputs=cost,
            #         updates=updates,
            #         givens={
            #             # x: T.concatenate([train_set_x[index * batch_size: (index + 1) * batch_size], added_train_x], axis=0),
            #             # y: T.concatenate([train_set_y[index * batch_size: (index + 1) * batch_size], added_train_y], axis=0)
            #             x: train_set_x[index * batch_size:(index + 1) * batch_size],
            #             y: train_set_y[index * batch_size:(index + 1) * batch_size]
            #         }
            # )


            # train_imb = theano.function(
            #     inputs=[index],
            #     outputs=[x.shape, y.sum(), y.shape],
            #     givens={
            #         x: train_set_x[index * batch_size: (index + 1) * batch_size],
            #         y: train_set_y[index * batch_size: (index + 1) * batch_size],
            #     })

            train_model = theano.function(
                    inputs=[index],
                    outputs=[],
                    updates=updates,
                    givens=expand_positive(seed, nnarray, train_pos_x, train_pos_y, batch_size, index,
                                           train_set_x,
                                           train_set_y, x, y),
                    mode=theano.Mode(linker='vm')
            )

            # compiling a Theano function that computes the mistakes that are made
            # by the model on a minibatch
            # test_model = theano.function(
            #         inputs=[index],
            #         outputs=classifier.errors(y),
            #         givens={
            #             x: test_set_x[index * batch_size:(index + 1) * batch_size],
            #             y: test_set_y[index * batch_size:(index + 1) * batch_size]
            #         }
            # )
            validate_model = theano.function(
                    inputs=[index],
                    outputs=classifier.errors(y),
                    givens={
                        x: valid_set_x[index * batch_size:(index + 1) * batch_size],
                        y: valid_set_y[index * batch_size:(index + 1) * batch_size]
                    }
            )
            valid_pred = theano.function(
                    inputs=[],
                    outputs=classifier.predict(y),
                    givens={
                        x: valid_set_x,
                        y: valid_set_y
                    },
            )
            #test_hidden = theano.function(
            #        inputs=[],
            #        outputs=[classifier.logRegressionLayer.input, y],
            #        givens={
            #            x: test_set_x,
            #            y: test_set_y,
            #        },
            #)
            test_pred = theano.function(
                    inputs=[],
                    outputs=classifier.predict(y),
                    givens={
                        x: test_set_x,
                        y: test_set_y
                    },
            )

            ###############
            # TRAIN MODEL #
            ###############
            print('... training')
            # go through this many
            # minibatches before checking the network
            # on the validation set; in this case we
            # check every epoch
            n_epochs = params['n_epochs']
            patience = params['patience']
            patience_increase = params['patience_increase']
            validation_frequency = min(n_train_batches, patience // patience_increase)
            best_validation_loss = np.inf
            best_validation_auc = 0.
            best_validation_pr = 0.
            best_iter = 0
            # test_score = 0.
            # test_auc = 0.

            start_time = timeit.default_timer()

            epoch = 0
            done_looping = False

            # val_aucs = []
            # Ws = []
            # iters = []
            while (epoch < n_epochs) and (not done_looping):
                epoch += 1
                for minibatch_index in range(n_train_batches):
                    # print("#positive - negative=\t",train_imb(minibatch_index))
                    minibatch_avg_cost = train_model(minibatch_index)
                    # print(minibatch_avg_cost[1])
                    # iteration number
                    iter = (epoch - 1) * n_train_batches + minibatch_index
                    # if iter % 50 == 0:
                    #     iters.append(iter)
                    #     # Ws.append(classifier.params[-2].get_value(borrow=True).T)
                    #     val_aucs.append(roc_auc_score(valid_set_y.get_value(borrow=True), valid_pred()))
                    #     test_aucs.append(roc_auc_score(test_set_y.get_value(borrow=True), test_pred()))

                    if (iter + 1) % validation_frequency == 0:
                        # compute zero-one loss on validation set
                        # validation_losses = [validate_model(i) for i
                        #                      in range(n_valid_batches)]
                        # this_validation_loss = np.mean(validation_losses)

                        # this_validation_auc = roc_auc_score(valid_set_y.get_value(borrow=True), valid_pred())
                        this_validation_pr = average_precision_score(valid_set_y.get_value(borrow=True), valid_pred())

                        ## add validation AUC to evaluate generalization
                        # print(
                        #         'epoch %i, minibatch %i/%i, validation error %f %%' %

                        #             epoch,
                        #             minibatch_index + 1,
                        #             n_train_batches,
                        #             this_validation_loss * 100.
                        #         )
                        # )
                        #
                        # print(
                        #         'epoch %i, minibatch %i/%i, validation error %f %%' %
                        #         (
                        #             epoch,
                        #             minibatch_index + 1,
                        #             n_train_batches,
                        #             this_validation_loss * 100.
                        #         )
                        # )

                        # if we got the best validation score until now
                        # add auc score for validation
                        # if this_validation_loss < best_validation_loss:
                        if this_validation_pr > best_validation_pr:
                        # if this_validation_auc > best_validation_auc:
                            #improve patience if loss improvement is good enough
                            if (
                                        # this_validation_loss < best_validation_loss *
                                        # params['improvement_threshold']
                                        # this_validation_auc > best_validation_auc *
                                        # (2-params['improvement_threshold'])
                                        this_validation_pr > best_validation_pr *
                                        (2-params['improvement_threshold'])
                            ):
                                patience = max(patience, iter * patience_increase)

                            # best_validation_loss = this_validation_loss

                            # this_validation_auc = roc_auc_score(valid_set_y.get_value(borrow=True), valid_pred())
                            # best_validation_auc = this_validation_auc

                            this_validation_pr = average_precision_score(valid_set_y.get_value(borrow=True), valid_pred())
                            best_validation_pr = this_validation_pr
                            # print(best_validation_pr)

                            # validation_aucs[fold] = best_validation_auc

                            validation_prs[fold] = best_validation_pr

                            # print('test ROC:', roc_auc_score(test_set_y.get_value(borrow=True), test_pred()))
                            # print('test PR:', average_precision_score(test_set_y.get_value(borrow=True), test_pred()))

                            best_iter = iter

                            if len(models) == model_index:
                                models.append([classifier, best_iter])
                            else:
                                models[model_index] = [classifier, best_iter]

                            # use all test set genes instead
                            # try:
                            #     test_auc = roc_auc_score(test_set_y.get_value(borrow=True), test_pred())
                            # except ValueError:
                            #     test_auc = np.nan

                            #print(('     epoch %i, minibatch %i/%i, test error of '
                            #       'best model %f %%') %
                            #      (epoch, minibatch_index + 1, n_train_batches,
                            #       test_score * 100.))
                            #print('best iter %i, best test genes auc %f' % (best_iter, test_auc))

                            #print(get_pickling_errors(classifier))
                            # with open('iter_%s__time_%s__best_model_training__%s.pkl' % (best_iter, t, diff_label), 'wb') as fout:
                            #    pickle.dump(classifier.logRegressionLayer.W.get_value(borrow=True), fout)
                            # plot the hidden unit against differential gene list
                            # last_hidden = test_hidden()

                    if patience <= iter:
                        done_looping = True
                        break

            end_time = timeit.default_timer()

            # print(('Optimization complete. Best validation score of %f %% '
            #        'obtained at iteration %i, with test performance %f %%, AUC score %f') %
            #       (best_validation_loss * 100., best_iter + 1, test_score * 100., test_auc))
            # from plot import pheatmap
            # pheatmap(last_hidden, diff_label, fold, t, test_auc)

            # learning curves
            # fig = plt.figure()
            # plt.plot(np.arange(len(val_aucs)), val_aucs, 'ro', label="validation")
            # plt.plot(np.arange(len(val_aucs)), test_aucs, 'b+', label="test")
            # plt.title("Test and validation AUC score trough n_epoch")
            # plt.legend(loc=4, borderaxespad=0., fontsize='small')
            # fig.savefig('%s_fold__time_%s__best_model_training__%s__maxs_%i.pdf' % (fold, t, diff_label, maxsamples))
            # visualize weights
            # from plot import weights_heatmap
            # weights_heatmap(Ws, iters, "time_%i_%s__maxs_%i" %(t, diff_label, maxsamples))

            print(('The code for file ' +
                   os.path.split(__file__)[1] +
                   ' ran for %.2fm' % ((end_time - start_time) / 60.)), file=sys.stderr)
        model_index += 1
        # models_performance.append(validation_aucs.mean())
        models_performance.append(validation_prs.mean())
        models_calibration.append(params)

    # store all the models and parameters
    with open('all_models__%s__%s.zip' % (diff_label, maxsamples), 'wb') as fout:
        dump((models, models_performance, models_calibration), fout) # Highest protocol

    best_model_index = np.argmax(np.array(models_performance))
    best_model = models[best_model_index]
    best_param = models_calibration[best_model_index]
    best_performance = models_performance[best_model_index]
    print(best_param)
    print(best_performance)

    params = best_param
    batch_size = params['batch_size']

    train_set_x, test_set_x, train_set_y, test_set_y = shared_data(seed, X_test, Y_test, X_train, Y_train)
    train_pos_index = np.asarray(np.where(train_set_y.get_value(borrow=True)==1)[0], dtype='int32')
    train_pos_x = train_set_x.get_value(borrow=True)[train_pos_index]
    train_pos_y = train_set_y.get_value(borrow=True)[train_pos_index]
    dist = distance()
    dist = dist(train_pos_x, train_pos_x)
    # 3 nearest neighbor for SMOTE
    nnarray = np.argsort(dist, axis=1)[:,1:3].transpose()
    nnarray = theano.shared(value = nnarray, name = 'k_nearest', borrow=True)
    train_pos_x = theano.shared(value=train_pos_x, name = "train_pos_x", borrow=True)
    train_pos_y = theano.shared(value=train_pos_y, name = "train_pos_y", borrow=True)
    # compute number of minibatches for training, validation and testing
    n_train_batches = train_set_x.get_value(borrow=True).shape[0] // params['batch_size']
    n_test_batches = test_set_x.get_value(borrow=True).shape[0] // params['batch_size']
    print('... building the model')
    index = T.lscalar()  # index to a [mini]batch
    # allocate symbolic variables for the data
    x = T.matrix('x')
    y = T.ivector('y')  # the labels are presented as 1D vector of [int] labels

    # train_aucs = np.zeros(repeat_times)
    train_prs = np.zeros(repeat_times)
    final_models = []
    for t in range(repeat_times):
        if not sparse_connect:
            # construct the MLP class of full connection NN
            classifier = MLP(
                    rng=rng,
                    input=x,
                    n_in=train_set_x.get_value(borrow=True).shape[1],
                    n_hidden=params['n_hidden'],
                    n_out=2
            )

        else:
            # construct the MyMLP class of sparse connection NN
            classifier = MyMLP(
                    rng=rng,
                    input=x,
                    seed=seed,
                    n_group=sparse_group,
                    n_in=features_num,
                    n_hidden=params['n_hidden'],
                    dropout_in=params['dropout_in'],
                    dropout_h=params['dropout_h'],
                    activation=params['hidden_activation'],
                    n_out=2,
                    pre_training=pre_training,
            )

        cost = (
            classifier.negative_log_likelihood(y)
            + params['L1_reg'] * classifier.L1
            + params['L2_reg'] * classifier.L2_sqr
        )

        batch_size = params['batch_size']
        # updates = Nesterov(cost, classifier.params, lr=params['learning_rate'], mu=params['rho'])
        updates = RMSprop(cost, classifier.params, lr=params['learning_rate'], rho=params['rho'], epsilon=params['epsilon'])
        
        # train_model = theano.function(
        #         inputs=[index],
        #         outputs=cost,
        #         updates=updates,
        #         givens={
        #             # x: T.concatenate([train_set_x[index * batch_size: (index + 1) * batch_size], added_train_x], axis=0),
        #             # y: T.concatenate([train_set_y[index * batch_size: (index + 1) * batch_size], added_train_y], axis=0)
        #             x: train_set_x[index * batch_size:(index + 1) * batch_size],
        #             y: train_set_y[index * batch_size:(index + 1) * batch_size]
        #         }
        # )

        train_model = theano.function(
                inputs=[index],
                outputs=[],
                updates=updates,
                givens=expand_positive(seed, nnarray, train_pos_x, train_pos_y, batch_size, index,
                                       train_set_x,
                                       train_set_y, x, y),
                mode=theano.Mode(linker='vm')
        )

        train_loss = theano.function(
            inputs=[index],
            outputs=classifier.errors(y),
            givens={
                x: train_set_x[index * batch_size: (index + 1) * batch_size],
                y: train_set_y[index * batch_size: (index + 1) * batch_size],
            },
        )
        train_pred = theano.function(
            inputs=[],
            outputs=classifier.predict(y),
            givens={
                x: train_set_x,
                y: train_set_y
            },
        )
        ###############
        # TRAIN MODEL #
        ###############
        print('final training %s' % t)
        n_epochs = params['n_epochs']
        validation_frequency = min(n_train_batches, patience // patience_increase)
        best_iter = 0
        start_time = timeit.default_timer()
        epoch = 0
        # final_model = best_model[0]
        min_losses = np.inf
        best_train_auc = 0.
        best_train_pr = 0.
        done = False
        while (epoch < n_epochs) and (not done):
            epoch += 1
            for minibatch in range(n_train_batches):
                train_model(minibatch)
                iter = (epoch - 1) * n_train_batches + minibatch
                if (iter+1) % validation_frequency == 0:
                    # train_losses = np.mean([train_loss(index) for index
                    #                         in range(n_train_batches)])

                    # if train_losses < min_losses:
                    #     min_losses = train_losses
                    train_pr = average_precision_score(train_set_y.get_value(borrow=True), train_pred())
                    # if train_auc > best_train_auc:
                        # train_auc = roc_auc_score(train_set_y.get_value(borrow=True), train_pred())
                    if train_pr > best_train_pr:
                        # best_train_auc = train_auc
                        best_train_pr = train_pr
                        # train_aucs[t] = best_train_auc
                        train_prs[t] = best_train_pr
                        if len(final_models) < t+1:
                            final_models.append(classifier)
                        else:
                            final_models[t] = classifier
                        # test_pred = theano.function(
                        #     inputs=[x, y],
                        #     outputs=final_model.predict(y))
                        # print('middle test auc: ', roc_auc_score(test_set_y.get_value(borrow=True), test_pred(test_set_x.get_value(borrow=True), test_set_y.get_value(borrow=True))))
                # use best iter from cross validation as a variant of early stopping
                # since samples increase 1.25 fold, increase early stopping best_iter to 1.25 fold
                # if (iter >= 1.25*best_model[1]):
                if (iter > best_model[1]):
                    done = True
                    break

    # best_index = np.argmax(train_aucs)
    best_index = np.argmax(train_prs)
    final_model = final_models[best_index]
    # print('final best train auc:', train_aucs[best_index])
    print('final best train pr:', train_prs[best_index])

    # test_pred = theano.function(
    #         inputs=[],
    #         outputs=final_model.predict(y),
    #         givens={
    #             x: test_set_x,
    #             y: test_set_y
    #         })
    test_pred = theano.function(
            inputs=[x, y],
            outputs=final_model.predict(y))
    print('final test auc: ', roc_auc_score(test_set_y.get_value(borrow=True), test_pred(test_set_x.get_value(borrow=True), test_set_y.get_value(borrow=True))))
    print('final test pr: ', average_precision_score(test_set_y.get_value(borrow=True), test_pred(test_set_x.get_value(borrow=True), test_set_y.get_value(borrow=True))))

    # detect what instance method cannot be pickled
    # remove in the __getstate__
    # recover in the __setstate__
    print(get_pickling_errors(final_model))

    with open('final_model__%s__%s.zip' % (diff_label, maxsamples), 'wb') as fout:
        dump((final_model, params), fout) # Highest protocol
    with open('final_model_func__%s__%s.pkl' % (diff_label, maxsamples), 'wb') as fout:
        pickle.dump(test_pred, fout) # Highest protocol

    return(predict_and_eval(test_set_x, test_set_y, diff_label, maxsamples))

def predict_and_eval(test_set_x, test_set_y, diff_label, ms):

    with open("final_model__%s__%s.zip" % (diff_label, ms), 'rb') as fin:
        classifier, params = load(fin)

    # issue: theano.gof.fg.MissingInputError: A variable that is an input to the graph was neither provided as an input to the function nor given a value. A chain of variables leading from this input to an output is [x, DimShuffle{1,0}.0, AdvancedSubtensor1.0, DimShuffle{1,0}.0, Elemwise{mul,no_inplace}.0, dot.0, Elemwise{add,no_inplace}.0, Elemwise{tanh,no_inplace}.0, Join.0, Elemwise{mul,no_inplace}.0, dot.0, Elemwise{add,no_inplace}.0, Softmax.0, AdvancedSubtensor.0]. This chain may not be unique
    # analysis:    # this is a funny thing
    # the same tensor variable with the same name
    # actually is two instance of tensor
    # which cannot be shared for a function
    # see code example below
    # one x created in the function is actually different
    # from the outer x which caused the above issue
    '''
    def test():
        x = theano.tensor.dscalar('x') # x1
        y = x**3
        f = theano.function([x], y) # f use the x1 
        return y, f(3)
    z, a = test()
    print(a)
    x = theano.tensor.dscalar('x') # f do not use the x2
    f2 = theano.function([x], z)
    print(f2(3))
    '''

    y = T.ivector('y')
    # pp(classifier.predict(y))
    test_pred = theano.function(
        inputs=[y],
        outputs=classifier.predict(y),
        givens={
            classifier.input: test_set_x, 
        }
    )
    pp(test_pred.maker.fgraph.outputs[0])
    print('recover', roc_auc_score(test_set_y.get_value(borrow=True), test_pred(test_set_y.get_value())))

    with open("final_model_func__%s__%s.pkl" % (diff_label, ms), 'rb') as fin:
        classifier = pickle.load(fin)
    auc = roc_auc_score(test_set_y.get_value(borrow=True), classifier(test_set_x.get_value(borrow=True), test_set_y.get_value(borrow=True)))

    print(auc)
    return auc
