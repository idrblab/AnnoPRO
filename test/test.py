import annopro.data_procession.profeat as pf
import os
import annopro.data as data
import shutil
from annopro.data_procession.blast import blastp
from annopro.data_procession import process
from annopro import main
from importlib import resources
os.chdir(os.path.dirname(__file__))
protein_fasta_file = "input-protein.fasta"
output_dir = "output1"
if os.path.isdir(output_dir):
    shutil.rmtree(output_dir)

# profeat 步骤ok
pf.run(protein_fasta_file, output_dir)

# blastp ok
with resources.path(data, "cafa4.dmnd") as path:
    blastp(path.absolute(), protein_fasta_file, f"{output_dir}/case.txt")

# TODO 报错，你的erro那个地方，请检查
process(protein_fasta_file, f"{output_dir}/output-protein.dat", output_dir)

# TODO 未能测试
main(output_dir, "0", True)


