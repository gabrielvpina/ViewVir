import os,subprocess

def interpro(outputFolder,ipFolder,CPU):

    for orf01 in os.listdir(outputFolder):
        if orf01.endswith("_ORFgc1.fasta"):
            orf01_path = os.path.join(outputFolder, orf01)

            interpro01 = [ipFolder, "--cpu", CPU, "--input", orf01_path, "-o", orf01_path.replace(".fasta","_CDD.tsv"), "-f", "TSV"]

            subprocess.run(interpro01, check=True)

    for orf05 in os.listdir(outputFolder):
        if orf05.endswith("_ORFgc5.fasta"):
            orf05_path = os.path.join(outputFolder, orf05)

            interpro05 = [ipFolder, "--cpu", CPU, "--input", orf05_path, "-o", orf05_path.replace(".fasta","_CDD.tsv"), "-f", "TSV"]

            subprocess.run(interpro05, check=True)

    for orf11 in os.listdir(outputFolder):
        if orf11.endswith("_ORFgc11.fasta"):
            orf11_path = os.path.join(outputFolder, orf11)

            interpro11 = [ipFolder, "--cpu", CPU, "--input", orf11_path, "-o", orf11_path.replace(".fasta","_CDD.tsv"), "-f", "TSV"]

            subprocess.run(interpro11, check=True)

    




