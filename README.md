# AnnoPRO
![AUR](https://img.shields.io/badge/license-MIT-blue.svg)
![python](https://img.shields.io/badge/python->=3.8-success.svg)
[![pypi](https://github.com/idrblab/AnnoPRO/actions/workflows/pypi.yml/badge.svg)](https://pypi.org/project/annopro/)
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
1. install compilers

dependency `lapjv` requires `g++` or other Cpp compiler, and annopro contains fortran extensional module and require `gfortran` or other fortran compiler. Here is an example of installing them on Ubuntu.

```bash
sudo apt install gcc g++ gfortran
# or you can install by conda in your virtual env
# command name is like 
# gcc: x86_64-conda_cos6-linux-gnu-cc
# g++: x86_64-conda_cos6-linux-gnu-c++
# gfortran: x86_64-conda_cos6-linux-gnu-gfortran
conda install gcc_linux-64 gxx_linux-64 gfortran_linux-64
```

2. install annopro

You can install it directly by `pip install annopro` or install from source code as following steps.
But you should install numpy first if you install it from source code because we need `numpy.f2py` to help us build fortran extension submodule.
```bash
git clone https://github.com/idrblab/AnnoPRO.git
cd AnnoPRO
conda create -n annopro python=3.8
conda activate annopro
pip install .
```

## Usage
- Use it as a terminal command. For all parameters, type `annopro -h`.
```bash
annopro -i test_proteins.fasta -o output
```
- Use it as a python executable package

```bash
python -m annopro -i test_proteins.fasta -o output
```

- Use it as a library to integrated with your project.
```python
from annopro import main
main("test_proteins.fasta", "output")
```

The result is displayed in the `./output/bp(cc,mf)_result.csv`.

**Notice**: if you use annopro for the first time, annopro will
automatically download required resources when they are used
(lazy download mechanism)

## Possible problems
1. pip is looking at multiple versions of XXX to determine which version is compatible with other requirements. this could take a while.

Your pip is latest, back to old version such as 20.2, or just add `--use-deprecated=legacy-resolver` param.

2. Argument mismatch when building source code.

Because your gfortran is latest and imcompatible,
edit setup.py and uncomment `-fallow-argument-mismatch` or 
just use a earlier version of gfortran such as 4.8.5, 8.4

## Contact
If any questions, please create an [issue](https://github.com/idrblab/AnnoPRO/issues/new/choose) on this repo, we will deal with it as soon as possible.
