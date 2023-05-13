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

![image](https://github.com/idrblab/AnnoPRO/assets/76670356/220c8fd7-6d81-4f4f-bf79-ce53d0a7a2bd)

## AnnoPRO architecture
* Encoding layers: Protein features was learned by CNNs and Protein similarity was learned by FCs.
* Decoding layers: LSTMs

![image](https://github.com/idrblab/AnnoPRO/assets/76670356/251065df-0c8a-4232-8855-cf7c40e6f3d6)

## Installation

You can install it directly by `pip install annopro` or install from source code as following steps.
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


## Contact
If any questions, please create an [issue](https://github.com/idrblab/AnnoPRO/issues/new/choose) on this repo, we will deal with it as soon as possible.
