import os
import argparse
import pandas as pd
from modules.tblfmt import process_diamondTbl,renameFasta
from modules.findorf import findorf
from modules.newORF import gc1_ORFs,gc5_ORFs,gc11_ORFs
from modules.plots import scatterPlot
from modules.IntProCD import interpro

# Importando dados
ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])
# Merge de todas as tabelas
allVirus = pd.concat([ncbiNames, ncbiSpecie])


parser = argparse.ArgumentParser()
parser.add_argument("-out","--outdir", type=str, help="Output directory name")
parser.add_argument("-scan","--interproscan", type=str, help="Interproscan executable path")
parser.add_argument("-cpu","--cpu", type=int, help="CPU for interproscan")
args = parser.parse_args()

########################### INPUT ###########################
# Output dir
vvfolder = str(args.outdir)
if vvfolder == "None":
    vvfolder = "ViewVir-results"

# Input contig
# Interproscan path
interpro_path = str(args.interproscan)

CPU = str(args.cpu)
if CPU == 0:
    CPU = 1


######################################## Processando contigs ##############################
# Renomear arquivo fasta original
# renameFasta(arquivo_input)

###########################################################################################

######################################## Processando diamond ##############################
dmndtable = "diamond.tsv"
process_diamondTbl(dmndtable,vvfolder)
renameFasta(vvfolder)

# Criar Orfs
findorf(vvfolder)

# Processing ORFs
gc1_ORFs(vvfolder)
gc5_ORFs(vvfolder)
gc11_ORFs(vvfolder)

# Scatter Plot
scatterPlot(vvfolder)

################################### CONSERVED DOMAINS ######################################

if interpro != "None":
    interpro(vvfolder,interpro_path,CPU)






