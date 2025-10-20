# ITTF
Repo for Genetic algorithm approach to assess correlation in TF-binding sites


## Overview
This repository houses the scripts used to conduct the in-silico experiments for the paper *Why do transcription factors target linear binding motifs? Mutational robustness and evolvability of information encoding strategies on DNA*

All scripts used to generate the results visualized in the paper are in `src`. For an overview on how
to use the scripts and how to regenerate the paper's results, please review the notebook: 📓 [View the Jupyter notebook on GitHub](https://github.com/ErillLab/ITTF/blob/main/tutorial.ipynb)


## Installation
### 1. Clone the repository

```bash
git clone https://github.com/ErillLab/ITTF.git
cd ITTF
```

### 2. Create and activate the Conda environment
```bash
conda env create -f environment.yml
conda activate ITTF
pip install dipwmsearch
```

## Usage
Type of execution: "parallel":
```bash
mpiexec -np <num_processes> python fitMotifPopulation.py
```
### Example
```bash
mpiexec -np 5 python fitMotifPopulation.py
```

Type of execution: "serial":
```bash
python fitMotifPopulation.py
```
