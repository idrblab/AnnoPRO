from consine_distances import load_data,MinMaxScaleClip
import pandas as pd
import pickle
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm

class Data_process():

    def __init__(self,protein_file,split_file,save_file,num,grid_file="/home/zhengly/PFmap/PFmap/data/data_grid.pkl",assess_file="/home/zhengly/PFmap/PFmap/data/row_asses.pkl",prosim_file="/home/zhengly/PFmap/PFmap/data/cafa4_del.csv",go_file="/home/zhengly/PFmap/PFmap/data/go.pkl"):
        '''
        protein_file 是生成的蛋白特征文件，
        split_file是需要包含所需信息的蛋白文件,
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
        df_data=load_data("/home/zhengly/PFmap/PFmap/data/cafa4_del.csv",1485)
        min_data=df_data.min()
        max_data=df_data.max()
        feature_data=load_data(self.protein_file,num=self.num)
        feature_data=(data - mean(data, axis=0)) / (std(data, axis=0)+1e-8)
        self.feature_data=feature_data
        proteins=pd.read_csv(self.protein_file,header=None)
        proteins.dropna(axis=0,inplace=True)
        proteins=proteins.iloc[:,0]
        self.proteins=list(proteins)
        self.data_grid=pickle.load(open(self.grid_file,"rb"))
        self.row_asses=pickle.load(open(self.assess_file,"rb"))
        prosim_data=load_data(self.prosim_file,1485)
        prosim_standard = MinMaxScaleClip(prosim_data)
        prosim_map=np.mat(prosim_standard)
        prosim_map = MinMaxScaleClip(prosim_map)
        self.prosim_map=np.mat(prosim_map)
        self.go=pickle.load(open(self.go_file,"rb"))



    def calculate_feature(self,row_num,size):
        F=open(self.split_file,'rb')
        train=pickle.load(F)
        protein_t=list(train['Proteins'])
        sequence_t=list(train['Sequence'])
        prop_t=list(train['Prop_annotations'])
        entry_t=list(train["Entry_name"])
        class_labels=['Composition','Autocorrelation','Physiochemical','Interaction','Quasi-sequence-order descriptors','PAAC for amino acid index set','Amphiphilic Pseudo amino acid composition']
        protein_all=[]
        sequences_all=[]
        prop_all=[]
        feature_all=[]
        prosim_all=[]
        entry_all=[]
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
                protein=proteins[i]
                row=0
                col=0
                protein_all.append(proteins[i])
                sequences_all.append(sequence_t[index_p])
                prop_all.append(prop_t[index_p])
                entry_all.append(entry_t[index_p])
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
        data_t=[protein_all,entry_all,sequences_all,feature_all,prosim_all,prop_all]
        data_t=pd.DataFrame(data_t)
        data_t=data_t.T
        data_t.columns=['Proteins','Entry_name','Sequence','Promap_feature','Protein_similary','Annotations']
        prop_annotations = []
        for i, row in tqdm(data_t.iterrows()):
            # Propagate annotations
            annot_set = set()
            annots = row.Annotations.split("; ")
            for go_id in annots:
                annot_set |= self.go.get_anchestors(go_id)
            annots = list(annot_set)
            prop_annotations.append(annots)
        data_t['Prop_annotations']=prop_annotations
        f = open(self.save_file, 'wb')
        pickle.dump(data_t, f)
        f.close() 
