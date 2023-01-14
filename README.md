# AnnoPRO
## Promap generation
* step 1: input proteins sequeces
* step 2: features extraction by Profeat
* step 3:  Feature pairwise distance calculation --> cosine, correlation, jaccard
* Step4: Feature 2D embedding --> umap, tsne, mds
* Step5: Feature grid arrangement --> grid, scatter
* Step5: Transform --> minmax, standard
![image](https://user-images.githubusercontent.com/76670356/204513203-2f0a430b-4b2c-4b1e-9587-3ee5a953150b.png)
## AnnoPRO architecture
* Encoding layers: Protein features was learned by CNNs and Protein similarity was learned by FCs.
* Decoding layers: LSTMs
![image](https://user-images.githubusercontent.com/76670356/204524869-31f558f0-0298-48c5-b4d2-3d5d087a2def.png)
## Installation
1. install gcc and  Profeat firstly.
```bash
git clone https://github.com/idrblab/AnnoPRO.git
cd AnnoPRO
conda create -n AnnoPRO python=3.8
conda activate AnnoPRO
pip install -r requirements.txt
unrar x profeat-new-version.rar
cd profeat-new-version
gfortran pro-des-35.f -o profeat
cp ../input-param.data input-param.data
```
The Profeat software source code uses Fortran language. It requires a related compilation environment (gcc) to run normally.<br /> 
2. generate proteins features <br />
upload your protein fasta file into the input-protein.dat
```bash
profeat
```
Then create a new file folder such as protein_A
```bash
cp input-protein.dat out-protein.dat protein_A
```
3. get database
```bash
cd ../
./AnnoPRO/get_data.sh
```
4. predict proteins functions
```bash
cd AnnoPRO
predict.sh protein_A
```
The result is displayed in the ./protein_A/bp(cc,mf)_result.csv
## Dependencies
AnnoPRO is tested to work under Python 3.8. and  cuda 11.2.
