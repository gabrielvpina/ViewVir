import os, glob, subprocess
from Bio import SeqIO, Seq
from Bio.Seq import Seq


def findorf(outputFolder):
    
    fileFasta = os.path.join(outputFolder, "*_nonDNA.fasta")
    nonDNA_files = glob.glob(fileFasta)

    for nonDNA_file in nonDNA_files:
        entrada = nonDNA_file
        nome = os.path.basename(nonDNA_file).replace("_nonDNA.fasta", "")
        saida = outputFolder

        # Executa orfipy - genetic code 01
        orfipy_gc1 = ["orfipy","--partial-3","--partial-5","--table", "1", "--outdir", saida, entrada, "--pep", f"{nome}_ORFgc1.fasta"]

        subprocess.run(orfipy_gc1, check=True)

        # Executa orfipy - genetic code 05
        orfipy_gc5 = ["orfipy","--partial-3","--partial-5", "--table", "5", "--outdir", saida, entrada, "--pep",  f"{nome}_ORFgc5.fasta"]

        subprocess.run(orfipy_gc5, check=True)

        # Executa orfipy - genetic code 11
        orfipy_gc11 = ["orfipy","--partial-3","--partial-5", "--table", "11", "--outdir", saida, entrada, "--pep",  f"{nome}_ORFgc11.fasta"]

        subprocess.run(orfipy_gc11, check=True)

        # Remove arquivos .log
        log_files = [os.path.join(outputFolder, file) for file in os.listdir(outputFolder) if file.endswith(".log")]
        for log_file in log_files:
            os.remove(log_file)
