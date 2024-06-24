import os
import argparse
import pandas as pd
from modules.tblfmt import viralFilter,renameFasta
from modules.findorf import findorf
from modules.newORF import gc1_ORFs,gc5_ORFs,gc11_ORFs
from modules.plots import scatterPlot
from modules.IntProCD import interpro
from modules.contigProcess import cap3,diamondTable,processDmndOut

# Importando dados
ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])
# Merge de todas as tabelas
allVirus = pd.concat([ncbiNames, ncbiSpecie])


parser = argparse.ArgumentParser()
parser.add_argument("-in","--input", type=str, help="Fasta with contigs")
parser.add_argument("-out","--outdir", type=str, help="Output directory name")
parser.add_argument("-vir","--viralDB",type=str, help=".dmnd file for diamond blastx")
parser.add_argument("-scan","--interproscan", type=str, help="Interproscan executable path")
parser.add_argument("-cpu","--cpu", type=int, help="CPU for interproscan")
args = parser.parse_args()

########################### INPUT ###########################

# Output dir
vvfolder = str(args.outdir)
if vvfolder == "None":
    vvfolder = "ViewVir-results"

# Input contig
inputContig = str(args.input)
if inputContig == "None":
    print("Please select fasta file")


# Interproscan path
interpro_path = str(args.interproscan)

# CPU
CPU = str(args.cpu)
if CPU == 0:
    CPU = 1

# Viral diamond database
viralDB = str(args.viralDB)
if viralDB == "None":
    print("Please select .dmnd file")



######################################## Processando contigs ##############################
# Assembly contigs
cap3(inputContig,vvfolder)

# Renomear arquivo fasta original
renameFasta(vvfolder)

######################################## Processando diamond ##############################

diamondTable(viralDB,vvfolder)

processDmndOut(vvfolder)

viralFilter(vvfolder)

######################################## Criar Orfs ########################################
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

