# ITTF
Genetic algorithm approach to assess correlation in TF-binding sites
- Under Development -

## Basic go-to guide:
Scripts related to MI vs IC motif mutational resillience (*Figures 1-4*):
- main_Figure1.py
- Divergence_Experiment1.py
- main_removeResidual.py 
- objects/organism.py
- objects/population.py

Scripts related to genome search experiment(Figure 5):
- genomeSearch.py
- removeResidualInformation_genomeSearch.py
- objects/genome.py

Scripts related to Coevolutionary Simulation:
- generatePerfectClassifier_CoevSimulation.py
- runCompetition_CoevSimulation.py
- objects/coevGenome.py
- objects/coevOrganism.py

Output contains saved output to 
generate figures using figure functions in src/Visualization. 
These output files will likely be 
transitioned from pickle objects to 
json files as the development of the 
jupyter notebook begins
