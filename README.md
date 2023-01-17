# AnnoPRO
![AUR](https://img.shields.io/badge/license-MIT-blue.svg)
![python](https://img.shields.io/badge/python->=3.8-success.svg)
![keras](https://img.shields.io/badge/keras-2.5.0-success.svg)
![PMID](https://img.shields.io/badge/PMID-Not%20available-red.svg)
## AnnoPRO generation
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
1. install diamond and compilers

dependency `lapjv` requires `g++` or other Cpp compiler, and annopro contains fortran extensional module and require `gfortran` or other fortran compiler. `diamond` will be invoked by annopro for blast. Here is an example of installing them on Ubuntu.

```bash
sudo apt install diamond-aligner
sudo apt install gcc g++ gfortran
```

2. install annopro

```bash
git clone https://github.com/idrblab/AnnoPRO.git
cd AnnoPRO
conda create -n annopro python=3.8
conda activate annopro
pip install -r requirements.txt
python setup.py install
```

## Usage

```
annopro input-protein.fasta -o output
```
The result is displayed in the `./output/bp(cc,mf)_result.csv`

## Note
AnnoPRO is tested to work under Python 3.8. and  cuda 11.2.

## Contact
If any questions, please create an [issue](https://github.com/idrblab/AnnoPRO/issues/new/choose) on this repo, we will deal with it as soon as possible.
