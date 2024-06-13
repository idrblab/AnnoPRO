# 要求的库与参数
import os
from tensorflow.keras.utils import Sequence, plot_model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import models, layers
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, CSVLogger, TensorBoard
# import matplotlib.pyplot as plt
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import (
    Input, Dense, Embedding, Conv2D, Flatten, Concatenate, TimeDistributed,
    MaxPool2D, Dropout, RepeatVector, Layer, Reshape, SimpleRNN, LSTM, BatchNormalization, GRU, Reshape,
    GlobalAveragePooling2D, GlobalMaxPooling2D, multiply, Permute, Add, Activation, Lambda, Permute, Multiply
)
import tensorflow as tf
import numpy as np
import pandas as pd
import math
import pickle
from sklearn.metrics import roc_curve, auc, matthews_corrcoef, precision_score, recall_score, roc_auc_score
from collections import deque, Counter
from tqdm import tqdm
from tensorflow.keras import backend as K

# AAINDEX = {'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5, 'Q': 6, 'E': 7, 'G': 8, 'H': 9, 'I': 10, 'L': 11,
#  'K': 12, 'M': 13, 'F': 14, 'P': 15, 'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20}
# MAXLEN = 2000


# def to_onehot(seq, start=0):
#     onehot = np.zeros((MAXLEN, 21), dtype=np.int32)
#     l = min(MAXLEN, len(seq))
#     for i in range(start, start + l):
#         onehot[i, AAINDEX.get(seq[i - start], 0)] = 1
#     onehot[0:start, 0] = 1
#     onehot[start + l:, 0] = 1
#     return onehot

# gpus = tf.config.experimental.list_physical_devices('GPU')

# tf.config.experimental.set_memory_growth(gpus[0], True)
os.environ["CUDA_VISIBLE_DEVICES"] = "2"

print(tf.test.is_gpu_available())


class Ontology(object):
    def __init__(self, filename='/data/Train/go-basic.obo', with_rels=False):
        self.ont = self.load(filename, with_rels)
        self.ic = None

    def has_term(self, term_id):
        return term_id in self.ont

    def get_term(self, term_id):
        if self.has_term(term_id):
            return self.ont[term_id]
        return None

    def get_anchestors(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        q = deque()
        q.append(term_id)
        while (len(q) > 0):
            t_id = q.popleft()
            if t_id not in term_set:
                term_set.add(t_id)
                for parent_id in self.ont[t_id]['is_a']:
                    if parent_id in self.ont:
                        q.append(parent_id)
        return term_set

    def get_parents(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        for parent_id in self.ont[term_id]['is_a']:
            if parent_id in self.ont:
                term_set.add(parent_id)
        return term_set

    def get_namespace_terms(self, namespace):
        terms = set()
        for go_id, obj in self.ont.items():
            if obj['namespace'] == namespace:
                terms.add(go_id)
        return terms

    def get_namespace(self, term_id):
        return self.ont[term_id]['namespace']

    def get_term_set(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        q = deque()
        q.append(term_id)
        while len(q) > 0:
            t_id = q.popleft()
            if t_id not in term_set:
                term_set.add(t_id)
                for ch_id in self.ont[t_id]['children']:
                    q.append(ch_id)
        return term_set


class DFGenerator(Sequence):
    def __init__(self, df, terms_dict, nb_classes, batch_size):
        self.start = 0
        self.size = len(df)
        self.df = df
        self.batch_size = batch_size
        self.nb_classes = nb_classes
        self.terms_dict = terms_dict

    def __len__(self):
        return np.ceil(len(self.df) / float(self.batch_size)).astype(np.int32)

    def __getitem__(self, idx):
        batch_index = np.arange(idx * self.batch_size, min(self.size, (idx + 1) * self.batch_size))
        df = self.df.iloc[batch_index]
        labels = np.zeros((len(df), self.nb_classes), dtype=np.int32)
        feature_data = []
        protein_si = []
        for i, row in enumerate(df.itertuples()):
            feature_data.append(list(row.Promap_feature))
            protein_si.append(list(row.Protein_similary))
            data_onehot = np.array(feature_data)
            data_si = np.array(protein_si)
            for t_id in row.Annotations:
                if t_id in self.terms_dict:
                    labels[i, self.terms_dict[t_id]] = 1
        self.start += self.batch_size
        return ([data_onehot, data_si], labels)

    def __next__(self):
        return self.next()

    def reset(self):
        self.start = 0

    def next(self):
        if self.start < self.size:
            batch_index = np.arange(
                self.start, min(self.size, self.start + self.batch_size))
            df = self.df.iloc[batch_index]
            labels = np.zeros((len(df), self.nb_classes), dtype=np.int32)
            feature_data = []
            protein_si = []
            for i, row in enumerate(df.itertuples()):
                feature_data.append(list(row.Promap_feature))
                protein_si.append(list(row.Protein_similary))
                data_onehot = np.array(feature_data)
                data_si = np.array(protein_si)
                for t_id in row.Annotations:
                    if t_id in self.terms_dict:
                        labels[i, self.terms_dict[t_id]] = 1
            self.start += self.batch_size
            return ([data_onehot, data_si], labels)
        else:
            self.reset()
            return self.next()


def load_data(file):
    f = open(file, 'rb')
    data = pickle.load(f)
    return data


def load_weight(model_path1, model_path2):
    model = load_model(model_path1)
    loaded_model = load_model(model_path2)
    old_weights = loaded_model.get_weights()
    now_weights = model.get_weights()
    cnt = 0
    for i in range(len(old_weights)):
        if old_weights[cnt].shape == now_weights[i].shape:
            now_weights[i] = old_weights[cnt]
            cnt = cnt + 1

    print(f'{cnt} layers weights copied, total {len(now_weights)}')
    model.set_weights(now_weights)
    model.save(model_path1)


# def plot_curve(history):
#     plt.figure()
#     x_range = range(0, len(history.history['loss']))
#     plt.plot(x_range, history.history['loss'], 'bo', label='Training loss')
#     plt.plot(x_range, history.history['val_loss'], 'b', label='Validation loss')
#     plt.title('Training and validation loss')
#     plt.legend()


def init_evaluate(data_size, batch_size, model_file, data_path, term_path):
    with open(term_path, 'rb') as file:
        terms_df = pickle.load(file)
    with open(data_path, 'rb') as file:
        data_df = pickle.load(file)
    if len(data_df) > data_size:
        data_df = data_df.sample(n=data_size)
    model = load_model(f'/data/model_parma/annopro/{model_file}.h5')
    data_file = data_path.split('/')[-1].split('.')[0]
    terms = terms_df.index.values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    nb_classes = len(terms)
    labels = np.zeros((len(data_df), nb_classes), dtype=np.int32)
    for i, row in enumerate(data_df.itertuples()):
        for go_id in row.Annotations:
            if go_id in terms_dict:
                labels[i, terms_dict[go_id]] = 1
    print('predict……')
    data_generator = DFGenerator(data_df, terms_dict, nb_classes, batch_size)
    data_steps = int(math.ceil(len(data_df) / batch_size))
    preds = model.predict(data_generator, steps=data_steps)
    # preds=diamond_score("/data/Train/valid.txt",preds,data_path,term_path)
    return terms, labels, preds, data_file


def fmeasure(real_annots, pred_annots):
    cnt = 0
    precision = 0.0
    recall = 0.0
    p_total = 0
    for i in range(len(real_annots)):
        if len(real_annots[i]) == 0:
            continue
        tp = set(real_annots[i]).intersection(set(pred_annots[i]))
        fp = pred_annots[i] - tp
        fn = real_annots[i] - tp
        tpn = len(tp)
        fpn = len(fp)
        fnn = len(fn)
        cnt += 1
        recall += tpn / (1.0 * (tpn + fnn))
        if len(pred_annots[i]) > 0:
            p_total += 1
            precision_x = tpn / (1.0 * (tpn + fpn))
            precision += precision_x
    recall /= cnt
    if p_total > 0:
        precision /= p_total
    fscore = 0.0
    if precision + recall > 0:
        fscore = 2 * precision * recall / (precision + recall)
    return fscore, precision, recall


def evaluate_annotations(labels_np, preds_np, terms):
    fmax = 0.0
    tmax = 0.0
    precisions = []
    recalls = []
    labels = list(map(lambda x: set(terms[x == 1]), labels_np))
    for t in range(1, 101):
        threshold = t / 100.0
        preds = preds_np.copy()
        preds[preds >= threshold] = 1
        preds[preds != 1] = 0
        #         fscore, pr, rc = fmeasure(labels, prop_annotations(preds, terms))
        fscore, pr, rc = fmeasure(labels, list(map(lambda x: set(terms[x == 1]), preds)))
        precisions.append(pr)
        recalls.append(rc)
        if fmax < fscore:
            fmax = fscore
            tmax = t
    preds = preds_np.copy()
    preds[preds >= tmax / 100.0] = 1
    preds[preds != 1] = 0
    mcc = matthews_corrcoef(labels_np.flatten(), preds.flatten())
    precisions = np.array(precisions)
    recalls = np.array(recalls)
    sorted_index = np.argsort(recalls)
    recalls = recalls[sorted_index]
    precisions = precisions[sorted_index]
    return fmax, tmax, recalls, precisions, mcc


def evaluate(model_file, data_path, data_size=8000, batch_size=16,
             term_path='/data/Train/new_terms/new_terms_CC.pkl'):
    ont = ['GO:0003674', 'GO:0008150', 'GO:0005575']
    namespace = ['molecular_function', 'biological_process', 'cellular_component', 'all']
    terms, labels, preds, data_file = init_evaluate(data_size, batch_size, model_file, data_path, term_path)
    # with open("/home/zhengly/promap/data/go.pkl", 'rb') as file:
    #     go = pickle.loads(file.read())
    # plt.figure(1, figsize=(16, 3))
    evaluate_info = f'{model_file}, {data_file}:\n'
    print(f'evaluate ……')
    # if i == 3:
    #     chose = np.ones(len(terms), dtype=bool)
    # else:
    #     go_set = go.get_namespace_terms(namespace[i])
    #     go_set.remove(ont[i])
    #     chose = list(map(lambda x: x in go_set, terms))
    # _terms = terms[chose]
    # _labels = labels[:, chose]
    # _preds = preds[:, chose]
    roc_auc = roc_auc_score(labels.flatten(), preds.flatten())
    fmax, alpha, recalls, precisions, mcc = evaluate_annotations(labels, preds, terms)

    # plt.subplot(1, 4, i + 1)
    # plt.plot(recalls, precisions, color='darkorange', lw=1, label=f'Fmax={fmax:0.3f}')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.0])
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.title(f'P-R curve of {namespace[i]}')
    # plt.legend(loc="lower right")
    evaluate_info += f'\t, {len(terms)}: fmax={fmax:0.3f}, mcc={mcc:0.3f}, roc_auc={roc_auc:0.3f}, precision={precisions[alpha]:0.3f}, recall={recalls[alpha]:0.3f}, threshold={alpha}\n'
    # plt.show()
    aucli = []
    fs = []
    with open(term_path, 'rb') as file:
        terms_df = pickle.load(file)
    tags = terms_df['tag']
    for i in range(1, 10):
        tag_select = tags == i
        _terms = terms[tag_select]
        _labels = labels[:, tag_select]
        _preds = preds[:, tag_select]
        aucli.append(roc_auc_score(_labels.flatten(), _preds.flatten()))
    #         (res) = evaluate_annotations(_labels, _preds, _terms)
    #         fs.append(res[0])
    # plt.figure()
    # plt.plot(range(1, len(aucli) + 1), aucli, lw=1, label=f'STD of auc={np.std(aucli):0.5f}')
    # #     plt.plot(range(1,len(fs)+1), fs, lw=1, color='orange', label=f'STD of fmax={np.std(fs):0.5f}')
    # plt.xlabel('Depth')
    # plt.legend(loc="lower right")
    # plt.ylim([0.0, 1.0])
    # plt.show()
    #     evaluate_info += f'\tauc_std={np.std(aucli):0.5f}, fmax_std={np.std(fs):0.5f}\n'
    evaluate_info += f'\tauc_std={np.std(aucli):0.5f}\n'
    print(aucli)
    print(evaluate_info)
    with open("logfile.json", "a") as file:
        file.write(evaluate_info)


def train_model(model_path, data_path, valid_path='none', epochs=60, batch_size=32, data_size=200000,
                terms_file='/data/Train/new_terms/new_terms_CC.pkl', early=True):
    model = load_model(model_path)
    #     model.summary()
    checkpointer = ModelCheckpoint(filepath=model_path, verbose=1, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=3, verbose=1)
    tbCallBack = TensorBoard(log_dir="./model", histogram_freq=1, write_grads=True)
    logger = CSVLogger('/data/log/reCRNN_bi.log')
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df.index.values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    nb_classes = len(terms)
    t = open(data_path, 'rb')
    data_df = pickle.load(t)
    # t.close()
    if len(data_df) > data_size:
        data_df = data_df.sample(n=data_size, random_state=312)

    if valid_path == 'none':
        valid_df = data_df.sample(frac=0.2, random_state=312)
        train_df = data_df[~data_df.index.isin(valid_df.index)]
    else:
        train_df = data_df
        v = open(valid_path, 'rb')
        valid_df = pickle.load(v)
        # v.close()
    valid_steps = int(math.ceil(len(valid_df) / batch_size))
    train_steps = int(math.ceil(len(train_df) / batch_size))
    train_generator = DFGenerator(train_df, terms_dict, nb_classes, batch_size)
    valid_generator = DFGenerator(valid_df, terms_dict, nb_classes, batch_size)
    #     训练模型
    if early:
        his = model.fit(
            train_generator,
            steps_per_epoch=train_steps,
            epochs=epochs,
            validation_data=valid_generator,
            validation_steps=valid_steps,
            max_queue_size=batch_size,
            workers=12,
            callbacks=[logger, checkpointer, earlystopper])
    else:
        his = model.fit(
            train_generator,
            steps_per_epoch=train_steps,
            epochs=epochs,
            validation_data=valid_generator,
            validation_steps=valid_steps,
            max_queue_size=batch_size,
            workers=12,
            callbacks=[logger, checkpointer])
    # plot_curve(his)



def pretrain_model(model_file, term_file='/data/Train/new_terms/new_terms_CC.pkl'):
    terms_df = pd.read_pickle(term_file)
    terms = terms_df.index.values.flatten()
    batch_size = 32
    # promap通道
    params = {
        'max_kernel': 129,
        'initializer': 'glorot_normal',
        'dense_depth': 0,
        'nb_filters': 512,
        'optimizer': Adam(lr=2e-4),
        'loss': 'binary_crossentropy'
    }
    nb_classes = len(terms)
    inp_hot = Input(shape=(39, 39, 7), dtype=np.float32)
    cnn1 = Conv2D(64, (3, 3), activation='relu', input_shape=(39, 39, 7))(inp_hot)
    pool1 = MaxPool2D((2, 2))(cnn1)
    cnn2 = Conv2D(128, (3, 3), activation='relu')(pool1)
    pool2 = MaxPool2D((2, 2))(cnn2)
    cnn_out = Flatten()(pool2)


    # protein_similary
    inp_similary = Input(shape=(142246), dtype=np.float32)
    encoded1 = Dense(2048, activation='relu')(inp_similary)
    encoded2 = Dense(1024, activation='relu')(encoded1)
    encoded3 = Dense(512, activation='relu')(encoded2)
    decoded1 = Dense(1024, activation='relu')(encoded3)
    decoded2 = Dense(2048, activation='relu')(decoded1)

    # concenate
    concat = Concatenate()([cnn_out, decoded2])
    net = BatchNormalization()(concat)
    net = Dropout(0.5)(net)
    dense =Dense(nb_classes, activation='sigmoid')(net)
    model = Model(inputs=[inp_hot, inp_similary], outputs=dense)
    model.compile(optimizer=params['optimizer'], loss=params['loss'])
    model.save(model_file)
    model.summary()


def CRNN_model(model_file, frozen=False, load=True, origin_path='/data/model_parma/annopro/0725pre_annopro.h5',term_file='/data/Train/new_terms/new_terms_CC.pkl'):
    terms_df = pd.read_pickle(term_file)
    terms = terms_df.index.values.flatten()
    batch_size = 32
    # promap通道
    params = {
        'max_kernel': 129,
        'initializer': 'glorot_normal',
        'dense_depth': 0,
        'nb_filters': 512,
        'optimizer': Adam(lr=2e-4),
        'loss': 'binary_crossentropy'
    }
    nb_classes = len(terms)
    inp_hot = Input(shape=(39, 39, 7), dtype=np.float32)
    cnn1 = Conv2D(64, (3, 3), activation='relu', input_shape=(39, 39, 7),trainable=frozen)(inp_hot)
    pool1 = MaxPool2D((2, 2),trainable=frozen)(cnn1)
    cnn2 = Conv2D(128, (3, 3), activation='relu',trainable=frozen)(pool1)
    pool2 = MaxPool2D((2, 2),trainable=frozen)(cnn2)
    cnn_out = Flatten()(pool2)


    # protein_similary
    inp_similary = Input(shape=(142246), dtype=np.float32)
    encoded1 = Dense(2048, activation='relu',trainable=frozen)(inp_similary)
    encoded2 = Dense(1024, activation='relu',trainable=frozen)(encoded1)
    encoded3 = Dense(512, activation='relu',trainable=frozen)(encoded2)
    decoded1 = Dense(1024, activation='relu',trainable=frozen)(encoded3)
    decoded2 = Dense(2048, activation='relu',trainable=frozen)(decoded1)

    # concenate
    concat = Concatenate()([cnn_out, decoded2])
    net = BatchNormalization()(concat)
    net = Dropout(0.5)(net)
    out =Dense(nb_classes, activation='relu')(net)
    drop = Dropout(0.5)(out)
    repeat = RepeatVector(11)(drop)
    gru1 = LSTM(256, activation='tanh', return_sequences=True)(repeat)
    gru2 = LSTM(256, activation='tanh', return_sequences=True)(gru1)
    gru3 = LSTM(256, activation='tanh', return_sequences=True)(gru2)
    net = layers.Flatten()(gru3)
    classify = layers.Dense(nb_classes, activation='sigmoid')(net)
    model = Model(inputs=[inp_hot, inp_similary], outputs=classify)
    model.compile(optimizer=params['optimizer'], loss=params['loss'])
    if load:
        loaded_model = load_model(origin_path)
        old_weights = loaded_model.get_weights()
        now_weights = model.get_weights()
        cnt = 0
        for i in range(len(old_weights)):
            if old_weights[cnt].shape == now_weights[i].shape:
                now_weights[i] = old_weights[cnt]
                cnt = cnt + 1
                print(i, cnt, len(now_weights), len(now_weights), len(old_weights))
        print(f'{cnt} layers weights copied, total {len(now_weights)}')
        model.set_weights(now_weights)
    model.save(model_file)
    model.summary()


# pretrain_model('/data/model_parma/annopro/0725pre_annopro.h5')
# train_model('/data/model_parma/annopro/0725pre_annopro.h5', '/data/Train/cafa5/annopro/time_train_annopro.pkl',
#             valid_path='/data/Train/cafa5/annopro/time_valid_annopro.pkl', epochs=200, data_size=60000)
evaluate('0725pre_annopro', '/data/Train/cafa5/annopro/time_valid_annopro.pkl', data_size=10000, batch_size=32, term_path='/data/Train/new_terms/new_terms_CC.pkl')
CRNN_model('/data/model_parma/annopro/0725annopro.h5',frozen=False,load=True,origin_path='/data/model_parma/annopro/0725pre_annopro.h5',term_file='/data/Train/new_terms/new_terms_CC.pkl')
train_model('/data/model_parma/annopro/0725annopro.h5','/data/Train/cafa5/annopro/time_train_annopro.pkl',valid_path='/data/Train/cafa5/annopro/time_valid_annopro.pkl', batch_size=32,epochs=200)
evaluate('0725annopro', '/data/Train/cafa5/annopro/time_valid_annopro.pkl', data_size=10000, batch_size=32, term_path='/data/Train/new_terms/new_terms_CC.pkl')
