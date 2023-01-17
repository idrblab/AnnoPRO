from collections import defaultdict
import pandas as pd


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