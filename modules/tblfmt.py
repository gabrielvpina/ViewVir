import os
import pandas as pd
from Bio import SeqIO

ncbiSpecie = pd.read_csv("data/NCBI_virSpecies.csv", names=["Species", "Genome.composition"])
ncbiNames = pd.read_csv("data/NCBI_virName.csv", names=["Species", "Genome.composition"])
allVirus = pd.concat([ncbiNames, ncbiSpecie])

def viralFilter(vvfolder):

    # result directory
    viewvirFolder = vvfolder
    if not os.path.exists(viewvirFolder):
        os.makedirs(viewvirFolder)

    for file_path in os.listdir(vvfolder): 
        if file_path.endswith("_proc.tsv"):
          
            
            file = pd.read_csv(os.path.join(vvfolder,file_path), sep="\t")

            # join
            file = pd.merge(file, allVirus, on="Species", how="left")

            # change order
            outrasCols = list(set(file.columns) - {"Genome.composition"})
            file = file[["QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "QCover",
                      "SubjTitle", "Species", "Genome.composition", "FullQueryLength"]]

            # filter length
            file = file[file["QseqLength"] >= 500]

            # remove duplicates
            file = file.drop_duplicates()
            
            inputBasename = os.path.basename(file_path)
            # out path
            output_file = os.path.join(viewvirFolder, inputBasename.replace("_proc.tsv", "_processed.tsv"))
    
            # write table
            file.to_csv(output_file, sep="\t", index=False)


            # non-dna filter
            tableRNA = file[file["Genome.composition"].isin(["ssRNA(-)", "ssRNA(+/-)", "ssRNA(+)", "ssRNA-RT", "RNA", "unknown", "NA"])]

            # non-dna table
            nonDNA_table = os.path.join(viewvirFolder, inputBasename.replace("_proc.tsv", "_nonDNA.tsv"))
            tableRNA.to_csv(nonDNA_table, sep="\t", index=False)


            #### Generate fasta file
            #vvfolder = "ViewVir-results"
            inputBasename = os.path.basename(vvfolder)
            files = os.listdir(vvfolder)

            output_pre = os.path.join(vvfolder, inputBasename.replace("_nonDNA.tsv", "_pre.tsv"))
            
            arquivo_PRE = tableRNA[["QuerySeq", "FullQueryLength"]]
            arquivo_PRE.to_csv(output_pre, sep="\t", index=False)

            # generate fasta 
            samp = os.path.basename(output_pre).replace("_pre.tsv", "")
            with open(output_pre, "r") as infile, open(f"{vvfolder}/{samp}_nonDNA.fasta", "w") as outfile:
                lines = infile.readlines()
                for line in lines[1:]:  # Ignorar a primeira linha
                    line = line.strip().replace('\t', '\n')
                    outfile.write(f">{line}\n")
    
            os.remove(output_pre)
            os.remove(os.path.join(vvfolder,file_path))



def renameFasta(outputFolder):

    for fasta in os.listdir(outputFolder):
        if fasta.endswith(".fasta"):
            input_fasta = os.path.join(outputFolder, fasta)
    
    registros_modificados = []

    
    for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), start=1):
        
        novo_id = f"contig_{i:02d}"
        record.id = novo_id
        record.description = novo_id  # Atualizar a descrição também

        
        registros_modificados.append(record)

    fastafile = os.path.join
    
    with open(input_fasta, "w") as output_handle:
        SeqIO.write(registros_modificados, output_handle, "fasta")

        
