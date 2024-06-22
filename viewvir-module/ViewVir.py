import os
import pandas as pd
import argparse
import re
from modules.tblfmt import process_diamondTbl
from modules.findorf import findorf
from modules.newORF import gc1_ORFs,gc5_ORFs,gc11_ORFs


parser = argparse.ArgumentParser()
parser.add_argument("-o","--outdir", type=str, help="Output directory")
args = parser.parse_args()

vvfolder = str(args.outdir)

if vvfolder == "None":
    vvfolder = "ViewVir-results"


# Importando dados
ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])

# Merge de todas as tabelas
allVirus = pd.concat([ncbiNames, ncbiSpecie])

# Processando dados do diamond
dmndtable = "diamond.tsv"
process_diamondTbl(dmndtable,vvfolder)

# Criar Orfs
findorf(vvfolder)

# Processing ORFs
gc1_ORFs(vvfolder)
gc5_ORFs(vvfolder)
gc11_ORFs(vvfolder)




