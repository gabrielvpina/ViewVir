import os
import argparse
import pandas as pd
from modules.tblfmt import viralFilter,renameFasta
from modules.findorf import findorf
from modules.newORF import gc1_ORFs,gc5_ORFs,gc11_ORFs
from modules.IntProCD import interpro
from modules.contigProcess import cap3,diamondTable,processDmndOut
from modules.makeblast import blastn,blastx
from modules.makePlot import scatterPlot, scatterPlotBLAST, generate_orf_plots, combine_html
    

parser = argparse.ArgumentParser()
parser.add_argument("-in","--input", type=str, help="Fasta with non-host contigs")
parser.add_argument("-out","--outdir", type=str, help="Output directory name")
parser.add_argument("-vir","--viralDB",type=str, help="Diamond database (.dmnd)")
parser.add_argument("-scan","--interproscan", type=str, help="Interproscan executable path -> /path/to/interproscan/./interproscan.sh")
parser.add_argument("-N","--blastn",type=str, help="BLASTn database path")
parser.add_argument("-X","--blastx",type=str, help="BLASTx database path")
parser.add_argument("-cpu","--cpu", type=int, help="CPU usage <int>")
parser.add_argument("-norf","--numORFs",type=int, help="Number of biggest ORFs selected")


args = parser.parse_args()

########################### INPUT ###########################
# Input contig
inputContig = str(args.input)
if inputContig == "None":
    print("Please select fasta file")

# Viral diamond database
viralDB = str(args.viralDB)
if viralDB == "None":
    print("Please select .dmnd file")

# Output dir
vvfolder = str(args.outdir)
if vvfolder == "None":
    vvfolder = "ViewVir-results"

# Interproscan path
interpro_path = str(args.interproscan)

# CPU
CPU = str(args.cpu)
if CPU == 0:
    CPU = 1

# Number of orfs
nORF = str(args.numORFs)
if nORF == "None":
    nORF = 2

# BLASTn
blastn_database = str(args.blastn)

# BLASTx
blastx_database = str(args.blastx)


######################################## Processando contigs ##############################
# Assembly contigs
cap3(inputContig,vvfolder)

# Renomear arquivo fasta original
renameFasta(vvfolder)

######################################## Processando diamond ##############################

diamondTable(viralDB,vvfolder,CPU)

processDmndOut(vvfolder)

viralFilter(vvfolder)

######################################## Criar Orfs ########################################
findorf(vvfolder)

# Processing ORFs
gc1_ORFs(vvfolder,nORF)
gc5_ORFs(vvfolder,nORF)
gc11_ORFs(vvfolder,nORF)


################################### CONSERVED DOMAINS ######################################

if interpro != "None":
    interpro(vvfolder,interpro_path,CPU)


################################# All Plots #################################
if blastn != "None":
    blastn(vvfolder,blastn_database,CPU)
    if blastx != "None":
        blastx(vvfolder,blastx_database,CPU)
        # Scatter Plot
        scatterPlotBLAST(vvfolder)
    else:
        scatterPlot(vvfolder)



suffixes = ['_ORFgc1.fasta', '_ORFgc5.fasta', '_ORFgc11.fasta']
output_file = 'orf_plots.html'

html_ORF_file = os.path.join(vvfolder,"orf_plots.html")
html_Combined_file = os.path.join(vvfolder,"ViewVir.html")

generate_orf_plots(vvfolder,html_ORF_file, suffixes)
scatterplot_html = scatterPlotBLAST(vvfolder)

# Combine HTMLs
script_dir = os.path.dirname(os.path.abspath(__file__))
image_path = os.path.join(script_dir, 'data', 'roundlogo.png')
combine_html(scatterplot_html ,html_ORF_file, html_Combined_file, image_path)



os.rmdir("temp")
