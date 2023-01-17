import annopro.data_procession._libprofeat as libprofeat
import annopro.data_procession._libprofeatconfig as libprofeatconfig
import os

def run(protein_fasta_file: str, output_dir: str):
    """
    run
    
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