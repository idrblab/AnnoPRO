# PFmap
## Promap generation
* step 1: input proteins sequeces
* step 2: features extraction by Profeat
* step 3:  Feature pairwise distance calculation --> cosine, correlation, jaccard
* Step4: Feature 2D embedding --> umap, tsne, mds
* Step5: Feature grid arrangement --> grid, scatter
* Step5: Transform --> minmax, standard
![image](https://user-images.githubusercontent.com/76670356/204513203-2f0a430b-4b2c-4b1e-9587-3ee5a953150b.png)
## PFmap architecture
* Encoding layers: Protein features was learned by CNNs and Protein similarity was learned by FCs.
* Decoding layers: LSTMs
![image](https://user-images.githubusercontent.com/76670356/204524869-31f558f0-0298-48c5-b4d2-3d5d087a2def.png)
## Installation
1. install gcc and  Profeat firstly.
```bash
git clone https://github.com/idrblab/PFmap.git
cd PFmap
conda create -n PFmap python=3.8
conda activate PFmap
pip install -r requirements.txt
unrar x ./profeat-new-version.rar
cd ./profeat-new-version
gfortran pro-des-35.f -o profeat
```
