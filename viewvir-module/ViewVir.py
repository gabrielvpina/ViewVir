import os
import pandas as pd
from modules.tblfmt import process_diamondTbl
from modules.findorf import findorf

vvfolder = "ViewVir-results"

# Importando dados
ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])

# Merge de todas as tabelas
allVirus = pd.concat([ncbiNames, ncbiSpecie])

# Processando dados
print("Enter the diamond table")
dmndtable = input()
process_diamondTbl(dmndtable,vvfolder)

# Criar Orfs
findorf(vvfolder)

