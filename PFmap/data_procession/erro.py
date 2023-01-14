from collections import defaultdict
from data_process import Data_process
import pandas as pd
def feature_to_csv(data_path):
    with open(data_path,"r") as file:
        pro_list = list()
        feature_list = defaultdict(list)
        f = file.readlines()
        ind=-1
        jishu=0
        for i in f:
            jishu+=1
            if jishu==2:
                pass
            else:
                # if "T" in i:
                if "|" in i:
                    # if ind>=0:
                    #     print(ind,protein,len(feature_list[pro_list[ind]]))
                    jishu=1       
                    # protein = i.split(" ")[1]
                    protein=i.split("|")[1].strip("\n"+" ")
                    pro_list.append(protein)
                    ind+=1
                    # print(protein,len(feature_list[pro_list[ind-1]]))                
                else:
                    feature=i.strip('\n').split(" ")
                    while "" in feature:
                        feature.remove("")
                    for fea in feature:
                        try:
                            feature_list[pro_list[ind]].append(float(fea))
                        except:
                            print(fea,protein,jishu)
    data = pd.DataFrame(feature_list)
    data=data.T
    return data
# feature_list=feature_to_csv("/home/zhengly/promap/promap/example/gdf/output-protein.dat")
# data = pd.DataFrame(feature_list)
# data.T.to_csv("/home/zhengly/promap/promap/example/gdf/dgf.csv", header=None)
