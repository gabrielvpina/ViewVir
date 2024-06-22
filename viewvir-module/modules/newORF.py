from Bio import SeqIO
from orffinder import orffinder

def findorf(outputFolder):

    outputFolder = "ViewVir-results"
    
    fileFasta = os.path.join(outputFolder, "*_nonDNA.fasta")
    nonDNA_files = glob.glob(fileFasta)

    sequence = SeqIO.read(nonDNA_files, "fasta")
    



