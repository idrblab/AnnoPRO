from annopro.data_procession.consine_distances import load_data, MinMaxScaleClip
import pandas as pd
import pickle

import numpy as np
from importlib import resources
from sklearn.metrics.pairwise import cosine_similarity
from annopro.data_procession.utils import Ontology
from annopro.data_procession.error import profeat_to_df
from fasta import FASTA


DATA_PACKAGE = "annopro.data"


class Data_process():

    def __init__(
            self,
            protein_file,
            proteins_fasta_file,
            save_file,
            num,
            grid_file="data_grid.pkl",
            assess_file="row_asses.pkl",
            prosim_file="cafa4_del.csv"):
        '''
        protein_file 是生成的蛋白特征文件，即profeat生成的文件
        split_file是需要包含所需信息的蛋白序列文件,
        save_file是生成的文件需要保存的位置,
        prosim_file是蛋白相似性比对的库
        grid_file，assess_file是使用的map的位置文件
        '''
        self.protein_file = protein_file
        self.split_file = proteins_fasta_file
        self.save_file = save_file
        self.grid_file = grid_file
        self.assess_file = assess_file
        self.prosim_file = prosim_file
        self.num = num
        self.__data__()

    def __data__(self):
        proteins_f = profeat_to_df(self.protein_file)
        proteins_f.dropna(axis=0, inplace=True)
        feature_data = proteins_f.iloc[:, :self.num]
        with resources.open_text(DATA_PACKAGE, "cafa4_del.csv") as cafa4_del:
            mia_data = load_data(cafa4_del, 1485)
        mia_data.columns = range(len(mia_data.columns))
        feature_data = (feature_data - mia_data.min()) / \
            ((mia_data.max() - mia_data.min()) + 1e-8)
        self.feature_data = feature_data
        self.proteins = list(proteins_f.index)
        with resources.open_binary(DATA_PACKAGE, self.grid_file) as gf:
            self.data_grid = pickle.load(gf)
        with resources.open_binary(DATA_PACKAGE, self.assess_file) as af:
            self.row_asses = pickle.load(af)
        with resources.open_text(DATA_PACKAGE, self.prosim_file) as pf:
            prosim_data = load_data(pf, 1485)
        prosim_standard = MinMaxScaleClip(prosim_data)
        prosim_map = np.mat(prosim_standard)
        prosim_map = MinMaxScaleClip(prosim_map)
        self.prosim_map = np.mat(prosim_map)
        self.go = Ontology()

    def calculate_feature(self, row_num, size):
        protein_seqs = FASTA(self.split_file).sequences
        class_labels = ['Composition', 'Autocorrelation', 'Physiochemical', 'Interaction',
                        'Quasi-sequence-order descriptors', 'PAAC for amino acid index set', 'Amphiphilic Pseudo amino acid composition']
        protein_all = []
        sequences_all = []
        feature_all = []
        prosim_all = []
        data_grid = self.data_grid
        row_asses = self.row_asses
        proteins = self.proteins
        feature_data = self.feature_data
        # 构建蛋白蛋白相似性
        result = np.mat(feature_data)
        consine_similarity = np.array(
            1)-cosine_similarity(result, self.prosim_map)
        # row, col = np.diag_indices_from(consine_similarity)
        # consine_similarity[row,col]=0
        print(proteins)
        print(protein_seqs)
        for i, protein in enumerate(proteins):
            if protein in protein_seqs:
                col_list = np.zeros(size)
                row = 0
                col = 0
                protein_all.append(protein)
                sequences_all.append(str(protein_seqs[protein].seq))
                prosim_all.append(consine_similarity[i])
                # 构建promap特征
                for j in range(len(data_grid['x'])):
                    channel = data_grid['subtype'][j]
                    index = class_labels.index(channel)
                    feature_index = row_asses[j]
                    row = j % row_num
                    col = j//row_num
                    col_list[col][row][index] = feature_data.iloc[i,
                                                                  feature_index]
                feature_all.append(col_list)
        data_t = [protein_all, sequences_all, feature_all, prosim_all]
        data_t = pd.DataFrame(data_t)
        data_t = data_t.T
        data_t.columns = ['Proteins', 'Sequence',
                          'Promap_feature', 'Protein_similary']
        data_t.to_pickle(self.save_file)
