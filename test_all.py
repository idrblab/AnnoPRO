import annopro.data_procession.profeat as pf
import os
import annopro.data as data
import shutil
from annopro.data_procession.blast import blastp
from annopro.data_procession import process
from annopro.prediction import predict
from importlib import resources

protein_fasta_file = "test_proteins.fasta"
output_dir = "test_output"
if os.path.isdir(output_dir):
    shutil.rmtree(output_dir)

# profeat 步骤ok
pf.run(protein_fasta_file, output_dir)

# blastp ok
with resources.path(data, "cafa4.dmnd") as path:
    blastp(path.absolute(), protein_fasta_file, f"{output_dir}/case.txt")

# ok
process(protein_fasta_file, f"{output_dir}/output-protein.dat", output_dir)

# ok
predict(output_dir, "0", True)


