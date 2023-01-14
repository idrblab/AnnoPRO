from consine_distances import load_data,MinMaxScaleClip
import pandas as pd
import pickle
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
from utils import Ontology
from erro import feature_to_csv
import os, sys

os.chdir(sys.path[0])


class Data_process():

    def __init__(self,protein_file,split_file,save_file,num,grid_file="../data/data_grid.pkl",assess_file="../data/row_asses.pkl",prosim_file="../data/cafa4_del.csv",go_file="../data/go.pkl"):
        '''
        protein_file 是生成的蛋白特征文件，即profeat生成的文件
        split_file是需要包含所需信息的蛋白序列文件,
        save_file是生成的文件需要保存的位置,
        prosim_file是蛋白相似性比对的库
        grid_file，assess_file是使用的map的位置文件
        '''
        self.protein_file=protein_file
        self.split_file=split_file
        self.save_file=save_file
        self.grid_file=grid_file
        self.assess_file=assess_file
        self.prosim_file=prosim_file
        self.go_file=go_file
        self.num=num
        self.__data__()

    def __data__(self):
        proteins_f=feature_to_csv(self.protein_file)
        proteins_f.dropna(axis=0,inplace=True)
        feature_data=proteins_f.iloc[:,:self.num]
        mia_data=load_data("../data/cafa4_del.csv",1485)
        mia_data.columns=range(len(mia_data.columns))
        feature_data=(feature_data - mia_data.min()) / ((mia_data.max() - mia_data.min()) + 1e-8)
        self.feature_data=feature_data
        self.proteins=list(proteins_f.index)
        self.data_grid=pickle.load(open(self.grid_file,"rb"))
        self.row_asses=pickle.load(open(self.assess_file,"rb"))
        prosim_data=load_data(self.prosim_file,1485)
        prosim_standard = MinMaxScaleClip(prosim_data)
        prosim_map=np.mat(prosim_standard)
        prosim_map = MinMaxScaleClip(prosim_map)
        self.prosim_map=np.mat(prosim_map)
        self.go=Ontology()


    def calculate_feature(self,row_num,size):
        train=pd.read_table(self.split_file,header=None)
        protein_t=[]
        sequence_t=[]
        # for i in train[train.index%2==0][0]:
        #     protein_t.append(i.strip(">"))
        # for j in train[train.index%2==1][0]:
        #     sequence_t.append(j)
        for i in train[0]:
            if "|" in i:
                protein_t.append(i.split("|")[1].strip("\n"+" "))
            else:
                sequence_t.append(i.strip("\n"))
        class_labels=['Composition','Autocorrelation','Physiochemical','Interaction','Quasi-sequence-order descriptors','PAAC for amino acid index set','Amphiphilic Pseudo amino acid composition']
        protein_all=[]
        sequences_all=[]
        feature_all=[]
        prosim_all=[]
        data_grid=self.data_grid
        row_asses=self.row_asses
        proteins=self.proteins
        feature_data=self.feature_data
        #构建蛋白蛋白相似性
        result=np.mat(feature_data)
        consine_similarity=np.array(1)-cosine_similarity(result,self.prosim_map)
        # row, col = np.diag_indices_from(consine_similarity) 
        # consine_similarity[row,col]=0
        for i in range(len(proteins)):
            if proteins[i] in protein_t:
                index_p=protein_t.index(proteins[i])
                col_list=np.zeros(size)
                row=0
                col=0
                protein_all.append(proteins[i])
                sequences_all.append(sequence_t[index_p])
                prosim_all.append(consine_similarity[i])
                #构建promap特征
                for j in range(len(data_grid['x'])):
                    channel=data_grid['subtype'][j]   
                    index=class_labels.index(channel)
                    feature_index=row_asses[j]
                    row=j%row_num
                    col=j//row_num
                    col_list[col][row][index]=feature_data.iloc[i,feature_index]
                feature_all.append(col_list)
        data_t=[protein_all,sequences_all,feature_all,prosim_all]
        data_t=pd.DataFrame(data_t)
        data_t=data_t.T
        data_t.columns=['Proteins','Sequence','Promap_feature','Protein_similary']
        f = open(self.save_file, 'wb')
        pickle.dump(data_t, f)
        f.close() 
