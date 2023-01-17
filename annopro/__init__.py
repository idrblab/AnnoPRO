import os
from tensorflow.keras.utils import Sequence
import matplotlib.pyplot as plt
from tensorflow.keras.models import load_model
import numpy as np
import pandas as pd
import math
import pickle
import annopro.data as data
import annopro.model_param as model_param
from annopro.focal_loss import BinaryFocalLoss
from annopro.data_procession.utils import NAMESPACES
from importlib import resources

__version__ = "0.0.1"


def main(file_path: str, used_gpu: str = None, with_diamond: bool = True):
    case_file = os.path.join(file_path, "case.txt")
    protein_file = os.path.join(file_path, "protein.pkl")
    if file_path == None:
        raise ValueError("Must provide the input fasta sequences.")
    if used_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = used_gpu
    for term_type in NAMESPACES.keys():
        init_evaluate(term_type=term_type,
                      data_path=protein_file,
                      case_file=case_file,
                      with_diamond=with_diamond,
                      output_dir=file_path)


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
        batch_index = np.arange(idx * self.batch_size,
                                min(self.size, (idx + 1) * self.batch_size))
        df = self.df.iloc[batch_index]
        labels = np.zeros((len(df), self.nb_classes), dtype=np.int32)
        feature_data = []
        protein_si = []
        for i, row in enumerate(df.itertuples()):
            feature_data.append(list(row.Promap_feature))
            protein_si.append(list(row.Protein_similary))
            data_onehot = np.array(feature_data)
            data_si = np.array(protein_si)
        self.start += self.batch_size
        return ([data_onehot, data_si])

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
            self.start += self.batch_size
            return ([data_onehot, data_si])
        else:
            self.reset()
            return self.next()


def diamond_score(diamond_scores_file, label, data_path, term_type):
    with resources.open_binary(data, "go.pkl") as file:
        go = pickle.loads(file.read())
    with resources.open_binary(data, "cafa_train.pkl") as file:
        train_df = pd.read_pickle(file)
    test_df = pd.read_pickle(data_path)
    annotations = train_df['Prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.Proteins] = i

    diamond_scores = {}
    with open(diamond_scores_file) as f:
        for line in f:
            it = line.strip().split("\t")
            if it[0] not in diamond_scores:
                diamond_scores[it[0]] = {}
            diamond_scores[it[0]][it[1]] = float(it[11])
    blast_preds = []

    for i, row in enumerate(test_df.itertuples()):
        annots = {}
        prot_id = row.Proteins
        # BlastKNN
        if prot_id in diamond_scores:
            sim_prots = diamond_scores[prot_id]
            allgos = set()
            total_score = 0.0
            for p_id, score in sim_prots.items():
                allgos |= annotations[prot_index[p_id]]
                total_score += score
            allgos = list(sorted(allgos))
            sim = np.zeros(len(allgos), dtype=np.float32)
            for j, go_id in enumerate(allgos):
                s = 0.0
                for p_id, score in sim_prots.items():
                    if go_id in annotations[prot_index[p_id]]:
                        s += score
                sim[j] = s / total_score
            ind = np.argsort(-sim)
            for go_id, score in zip(allgos, sim):
                annots[go_id] = score
        blast_preds.append(annots)
    with resources.open_binary(data, f"terms_{NAMESPACES[term_type]}.pkl") as term_path:
        terms = pd.read_pickle(term_path)
    terms = terms['terms'].values.flatten()
    alphas = {NAMESPACES['mf']: 0.55,
              NAMESPACES['bp']: 0.6, NAMESPACES['cc']: 0.4}

    for i in range(0, len(label)):
        annots_dict = blast_preds[i].copy()
        for go_id in annots_dict:
            annots_dict[go_id] *= alphas[go.get_namespace(go_id)]
        for j in range(0, len(label[0])):
            go_id = terms[j]
            label[i, j] = label[i, j]*(1 - alphas[go.get_namespace(go_id)])
            if go_id in annots_dict:
                label[i, j] = label[i, j] + annots_dict[go_id]
    return label


def plot_curve(history):
    plt.figure()
    x_range = range(0, len(history.history['loss']))
    plt.plot(x_range, history.history['loss'], 'bo', label='Training loss')
    plt.plot(x_range, history.history['val_loss'],
             'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()


def init_evaluate(term_type, data_path, case_file, output_dir:str, data_size=8000, batch_size=16,
                  with_diamond: bool = True):
    with resources.open_binary(data, f"terms_{NAMESPACES[term_type]}.pkl") as file:
        terms_df = pickle.load(file)
    with open(data_path, 'rb') as file:
        data_df = pickle.load(file)
    if len(data_df) > data_size:
        data_df = data_df.sample(n=data_size)
    data_df.index = range(len(data_df))
    with resources.path(model_param, f"{term_type}.h5") as model_file_path:
        model = load_model(
            model_file_path.absolute(),
            custom_objects={"focus_loss": BinaryFocalLoss})
    proteins = data_df["Proteins"]
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    nb_classes = len(terms)
    data_generator = DFGenerator(data_df, terms_dict, nb_classes, batch_size)
    data_steps = int(math.ceil(len(data_df) / batch_size))
    preds = model.predict(data_generator, steps=data_steps)
    if with_diamond:
        preds = diamond_score(case_file, preds, data_path, term_type=term_type)
    # label_di=defaultdict(list)
    protein = []
    go_terms = []
    score = []
    for i in range(len(preds)):
        for j in range(len(preds[i])):
            if preds[i][j] > 0:
                protein.append(proteins[i])
                go_terms.append(terms[j])
                score.append(preds[i][j])
    res = [protein, go_terms, score]
    res = pd.DataFrame(res)
    res = res.T
    res.columns = ['Proteins', 'GO-terms', 'Scores']
    res.sort_values(by='Scores', axis=0, ascending=False, inplace=True)
    result_file = os.path.join(output_dir, f"{term_type}_result.csv")
    res.to_csv(result_file, sep=',', index=False, header=True)
    return res