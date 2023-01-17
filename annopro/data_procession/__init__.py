from annopro.data_procession.data_predict import Data_process


def process(proteins_fasta_file: str, profeat_file: str, save_file: str):
    if proteins_fasta_file == None:
        raise ValueError("Must provide the input fasta sequences.")

    data = Data_process(protein_file=profeat_file,
                        proteins_fasta_file=proteins_fasta_file, 
                        save_file=save_file, num=1484)
    data.calculate_feature(row_num=39, size=(39, 39, 7))
