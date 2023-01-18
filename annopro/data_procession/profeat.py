import annopro.data_procession._libprofeat as libprofeat
import annopro.data_procession._libprofeatconfig as libprofeatconfig
import os
from collections import defaultdict
import pandas as pd


def run(protein_fasta_file: str, output_dir: str):
    """
    run profeat
    
    """
    if not output_dir.endswith(os.sep):
        output_dir += os.sep
    config_dir = os.path.dirname(libprofeatconfig.__file__) + os.sep
    output_dir_len = len(output_dir)
    config_dir_len = len(config_dir)
    assert output_dir_len <= 100 and config_dir_len <= 300, "Too long directory nme"
    os.mkdir(output_dir)
    libprofeat.run(
        protein_fasta_file, 
        output_dir, 
        output_dir_len,
        config_dir,
        config_dir_len)


def profeat_to_df(data_path: str) -> pd.DataFrame:
    """
    Read profeat features as `pd.DataFrame`.
    """
    with open(data_path, "r") as file:
        protein_list = list()
        feature_list = defaultdict(list)
        protein_line_index = -1
        line: str
        for rol_index, line in enumerate(file, 1):
            line = line.strip()
            if protein_line_index == rol_index -1:
                continue

            if line.startswith(">"):
                protein_line_index = rol_index
                protein = line[1:]
                protein_list.append(protein)
            else:
                protein = protein_list[-1]
                for feature in line.split():
                    try:
                        feature_list[protein].append(float(feature))
                    except ValueError:
                        print(f"Invalid feature {feature} for {protein} at line {rol_index}")
    
    return pd.DataFrame(feature_list).T