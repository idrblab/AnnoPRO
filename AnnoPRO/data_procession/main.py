from data_predict import Data_process
# from absl import flags
# import tensorflow.compat.v1 as tf
import argparse
import os, sys

os.chdir(sys.path[0])
parser = argparse.ArgumentParser(description='Arguments for main.py')
parser.add_argument('--file_path', default=None, type=str)
args = parser.parse_args()


if args.file_path == None:
    raise ValueError("Must provide the input fasta sequences.")

protein_file=os.path.join(args.file_path,"output-protein.dat")
split_file=os.path.join(args.file_path,"input-protein.dat")
save_file=os.path.join(args.file_path,"protein.pkl")


data=Data_process(protein_file=protein_file,split_file=split_file,save_file=save_file,num=1484)
data.calculate_feature(row_num=39,size=(39,39,7))

