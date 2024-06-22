import os
import pandas as pd
import argparse
import re
from modules.tblfmt import process_diamondTbl
from modules.findorf import findorf
from modules.processBED import bedgc1


parser = argparse.ArgumentParser()
parser.add_argument("-o","--outdir", type=str, help="Output directory")
args = parser.parse_args()

vvfolder = str(args.outdir)

# Importando dados
ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])

# Merge de todas as tabelas
allVirus = pd.concat([ncbiNames, ncbiSpecie])

# Processando dados
#print("Enter the diamond table")
#dmndtable = input()
dmndtable = "diamond.tsv"
process_diamondTbl(dmndtable,vvfolder)

# Criar Orfs
findorf(vvfolder)

# Processing ORFs
bedgc1(vvfolder)




