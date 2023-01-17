import numpy as np
from scipy.sparse import data
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
# from numpy import *


def load_data(file, num):
    data = pd.read_csv(file, header=None)
    data.dropna(axis=0, inplace=True)
    data_noprotein = data.iloc[:, 1:num]
    return data_noprotein


def MinMaxScaleClip(data):
    data_standard = (data - data.min()) / ((data.max() - data.min()) + 1e-8)
    return data_standard


def StandardScaler(data):
    data_standard = (data - mean(data, axis=0)) / (std(data, axis=0)+1e-8)
    return data_standard


def calculate_consine(file, num, standard_type=['minmax', 'standard'], var_thr=np.e**(-8)):
    data = load_data(file, num)
    if standard_type == 'minmax':
        data = MinMaxScaleClip(data)
    elif standard_type == 'standard':
        data = StandardScaler(data)
    result = np.mat(data)
    result = np.transpose(result)
    consine_similarity = np.array(1)-cosine_similarity(result)
    for i in range(len(consine_similarity)):
        for j in range(len(consine_similarity[0])):
            if consine_similarity[i][j] < var_thr:
                consine_similarity[i][j] = 0
    return consine_similarity


# if __name__ == '_main_':
#     def calculate_consine(file, standard_type=['minmax', 'standard'], var_thr=e**(-8)):
#         data = load_data(file)
#         if standard_type == 'minmax':
#             data = MinMaxScaleClip(data)
#         elif standard_type == 'standard':
#             data = StandardScaler(data)
#         result = np.mat(data)
#         result = np.transpose(result)
#         consine_similarity = np.array(1)-cosine_similarity(result)
#         for i in range(len(consine_similarity)):
#             for j in range(len(consine_similarity[0])):
#                 if consine_similarity[i][j] < var_thr:
#                     consine_similarity[i][j] = 0
#         return consine_similarity
