from annopro.data_procession.data_predict import Data_process
import os


def process(protein_fasta_file: str, profeat_file: str, output_dir: str):
    if protein_fasta_file == None:
        raise ValueError("Must provide the input fasta sequences.")

    save_file = os.path.join(output_dir, "protein.pkl")
    data = Data_process(protein_file=profeat_file,
                        split_file=protein_fasta_file, 
                        save_file=save_file, num=1484)
    data.calculate_feature(row_num=39, size=(39, 39, 7))
