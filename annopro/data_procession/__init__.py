from annopro.data_procession.data_predict import Data_process
import os


def process(file_path: str):
    if file_path == None:
        raise ValueError("Must provide the input fasta sequences.")

    protein_file = os.path.join(file_path, "output-protein.dat")
    split_file = os.path.join(file_path, "input-protein.dat")
    save_file = os.path.join(file_path, "protein.pkl")
    data = Data_process(protein_file=protein_file,
                        split_file=split_file, save_file=save_file, num=1484)
    data.calculate_feature(row_num=39, size=(39, 39, 7))
