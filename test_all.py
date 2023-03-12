import profeat as pf
import os
import annopro.resources as resources
import shutil
from annopro.data_procession.blast import blastp
from annopro.data_procession import process
from annopro.prediction import predict

protein_fasta_file = "test_proteins.fasta"
output_dir = "test_output"
if os.path.isdir(output_dir):
    shutil.rmtree(output_dir)

# profeat 步骤ok
pf.run(protein_fasta_file, output_dir)

# blastp ok
diaonds = f"{output_dir}/diamond_scores.txt"
blastp(resources.get_resource_path("cafa4.dmnd"), protein_fasta_file, diaonds)

# ok
promap_file = f"{output_dir}/promap_features.pkl"
process(protein_fasta_file, f"{output_dir}/output-protein.dat", promap_file)

# ok
predict(output_dir, promap_file, "0", diamond_scores_file=diaonds)


