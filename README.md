# ITTF

Code and data for the paper:

> Mekasha, E., Mascolo, E., & Erill, I. (2026). *Why don't transcription factors exploit
> inter-nucleotide correlations in binding motifs? Evolution of information encoding
> strategies on DNA.*

## Overview

A DNA sequence motif can encode binding information in two ways. One, through nucleotide
conservation at individual positions, or two, through correlations between positions. Real
transcription factor binding motifs appear to store only a small fraction of their
information as correlations. This repository features the in-silico  experiments
used to investigate why.

Populations of simulated organisms are generated using a genetic algorithm with a certain target
total information (encoded through either PI/MI). They are then compared via mutational
robustness, genome-search performance, direct competition, and the structure of their
genotype-phenotype maps.

For a guided walkthrough that reproduces the paper's results, start with the notebook:
📓 [tutorial.ipynb](https://github.com/ErillLab/ITTF/blob/main/tutorial.ipynb)

## Repository map

| Path | What it holds |
|---|---|
| [`src/`](src/) | All simulation and figure-generation code. See the two sections below. |
| [`src/Objects/`](src/Objects/) | Core classes: `MotifPopulation`, `MotifModel`, `Genome`, `coevOrganism`/`coevPopulation`, `GeneticCode`, `MI_VS_SI`, `MotifFLA`. |
| [`src/settings_*.json`](src/) | Run parameters, one file per simulation entry point. |
| [`Notebook/`](Notebook/) | Small curated inputs for the tutorial. The exact evolved populations (`Populations/`) and per-generation trajectories (`Results/`) behind the paper's figures are stored in the repository|
| [`Output/`](Output/) | Raw simulation output, written by the scripts. Also holds the generated figure PDFs in [`Output/Figures/`](Output/Figures/).. |
| [`FigureImports/`](FigureImports/) | Drawn schematic panels that are composited into the final figure plates. |
| `tutorial.ipynb` | Annotated walkthrough of every experiment and figure. |



## Terminology: paper vs. code

The code predates the paper's notation, so the names differ:

| Paper | Code | Meaning |
|---|---|---|
| SI | `IC`, `PI`, `propInformation` | Single-locus information — conservation at individual positions |

So a population named `MIpop` is one evolved to encode its information as correlations,
and `ICpop` is one evolved to encode it positionally.

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

`dipwmsearch` is installed separately because it is not available on the conda

## Running the simulations

**All commands run from the repository root.** The scripts inputs are detailed in the settings files
as `./src/settings_*.json` and import themselves as `src.<module>`, so they must be
called with `-m` from the root rather than executed as files from inside `src/`.

Each script has its own settings file (edit the JSON, then run):

| Script | Settings file | What it does |
|---|---|---|
| `python -m src.fitMotifPopulation` | `src/settings_fitMotifPop.json` | Evolves a motif population to a target information content under a chosen encoding strategy. |
| `python -m src.fitCoevPopulation` | `src/settings_fitCoevPop.json` | Coevolutionary simulation of TF and genome, used for the competition and grid experiments. |
| `python -m src.genomeSearch` | `src/settings_GenomeExperiment.json` | Scores an evolved motif model against a genome and reports AUPRC. |
| `python -m src.informationDecay` | `src/settings_informationDecay.json` | Applies mutations to a population and tracks information loss. |
| `python -m src.main_removeResidual` | `src/settings_residual.json` | Equalizes residual information between two populations. |

### Parallel execution

`fitCoevPopulation` is the only MPI-enabled entry point. Set `"run_mode": "parallel"` in
`src/settings_fitCoevPop.json` and run with `mpiexec`

```bash
mpiexec -np 5 python -m src.fitCoevPopulation
```

With `"run_mode": "serial"` the same script runs single-process. 

## Regenerating the figures


| SubFig | Script | Role |
|---|---|---|
| Panel | [`src/figureFunctions.py`](src/figureFunctions.py) | One function per panel. Returns a figure and saves nothing. Most functions take an `axis` so the panel can be drawn into a plate |
| Whole Figure | [`src/mainFigures.py`](src/mainFigures.py) | One function per *paper* figure. Lays out panels with gridspec, adds the A/B/C labels, loads schematic PNGs, and saves to `Output/Figures/`. Contains no plotting logic of its own. |



## Citation

If you use this code, please cite:



## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).
