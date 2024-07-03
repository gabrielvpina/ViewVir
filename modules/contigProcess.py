import os
import subprocess
import pandas as pd

def cap3(inputContig,outputFolder):
    # directory w/ results
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

    # copy original fasta in result directory
    copy_command = ["cp", inputContig, outputFolder]
    subprocess.run(copy_command, check=True)

    for cpContig in os.listdir(outputFolder):
        if cpContig.endswith(".fasta"):
            # run CAP3
            cap3_command = ["cap3", os.path.join(outputFolder,cpContig)]
            subprocess.run(cap3_command, check=True)

    for contigs in os.listdir(outputFolder):
        if contigs.endswith(".cap.contigs"):
            sample = contigs.replace(".fasta.cap.contigs","")
            print(f"Processando amostra {sample}")

            for singlets in os.listdir(outputFolder):
                if singlets.endswith(".cap.singlets"):

                    cat_command = f"cat {os.path.join(outputFolder,contigs)} {os.path.join(outputFolder,singlets)} >> {os.path.join(outputFolder, f'{sample}_merged.fasta')}"
                    print(cat_command)
                    subprocess.run(cat_command, shell=True, check=True)
                    # merge the outputs of cap3
    
    
    log_files = [os.path.join(outputFolder, file) for file in os.listdir(outputFolder)
                 if file.endswith((".links",".qual",".info",".contigs",".singlets",".ace"))]
    for log_file in log_files:
        os.remove(log_file)
    # remove original fasta file
    os.remove(os.path.join(outputFolder,inputContig))


def diamondTable(viralDB, outputFolder,CPU):
    # debugs
    if not os.path.exists(outputFolder):
        raise FileNotFoundError(f"O diretório {outputFolder} não existe.")

    for inputContig in os.listdir(outputFolder):
        if inputContig.endswith("_merged.fasta"):
            sample = inputContig.replace("_merged.fasta", "")

            command_dmnd = f"diamond blastx -d {viralDB} -q {os.path.join(outputFolder, inputContig)} --threads {CPU} --outfmt 6 qseqid sseqid qlen slen pident evalue qcovhsp stitle full_qseq --max-target-seqs 1 --out {os.path.join(outputFolder, sample + '_diamond.tsv')}"
            
            # run diamond
            try:
                subprocess.run(command_dmnd, shell=True, check=True)
                print(f"Processado: {sample}")
            except subprocess.CalledProcessError as e:
                print(f"Erro ao processar {sample}: {e}")

def processDmndOut(outputFolder):
    for table in os.listdir(outputFolder):
        if table.endswith("_diamond.tsv"):
            sample = table.replace('_diamond.tsv', '_proc.tsv')

            entrada = os.path.join(outputFolder, table)
            saida = os.path.join(outputFolder, sample)
            
            # replace w/ sed - important
            sed_com = f"sed 's/\\[/\\t/g; s/\\]//g' {entrada}"
            
            # headers of table
            newcols = "QuerySeq\tSubjectSeq\tQseqLength\tSseqLength\tPident\tEvalue\tQCover\tSubjTitle\tSpecies\tFullQueryLength\n"
            
            try:
                # run sed
                result = subprocess.run(sed_com, shell=True, check=True, capture_output=True, text=True)
                
                # write headers in final table
                with open(saida, 'w') as out_file:
                    out_file.write(newcols)
                    out_file.write(result.stdout)
                
                # remove original table
                os.remove(entrada)
            except subprocess.CalledProcessError as e:
                print(f"Erro ao processar o arquivo {entrada}: {e}")
            except Exception as e:
                print(f"Erro inesperado ao processar o arquivo {entrada}: {e}")

