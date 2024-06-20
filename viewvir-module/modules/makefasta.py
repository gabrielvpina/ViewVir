import os
import pandas as pd

def makefasta(pre_file):

    inputBasename = os.path.basename(pre_file)

    output_pre = os.path.join(pre_file, inputBasename.replace("_nonDNA.tsv", "_pre.tsv"))
    arquivo_PRE = tableRNA[["QuerySeq", "FullQueryLength"]]
    arquivo_PRE.to_csv(output_pre, sep="\t", index=False)

    # Processar _pre.tsv para gerar RNA-virus.fasta
    samp = os.path.basename(output_pre).replace("_pre.tsv", "")
    with open(output_pre, "r") as infile, open(f"{pre_file}/{samp}_nonDNA.fasta", "w") as outfile:
        lines = infile.readlines()
        for line in lines[1:]:  # Ignorar a primeira linha
            line = line.strip().replace('\t', '\n')
            outfile.write(f">{line}\n")
    
    os.remove(output_pre)
