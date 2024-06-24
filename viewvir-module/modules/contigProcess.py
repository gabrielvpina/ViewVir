import os
import subprocess
import pandas as pd

def cap3(inputContig,outputFolder):
    # Pasta para os resultados
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

    # Copiar o arquivo de entrada para a pasta de saída
    copy_command = ["cp", inputContig, outputFolder]
    subprocess.run(copy_command, check=True)

    for cpContig in os.listdir(outputFolder):
        if cpContig.endswith(".fasta"):
            # Executar CAP3
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
                    # O resultado é o arquivo sample_merged.fasta + o arquivo .fasta original
    
    
    log_files = [os.path.join(outputFolder, file) for file in os.listdir(outputFolder)
                 if file.endswith((".links",".qual",".info",".contigs",".singlets",".ace"))]
    for log_file in log_files:
        os.remove(log_file)
    # Remover arquivo original da pasta
    os.remove(os.path.join(outputFolder,inputContig))


def diamondTable(viralDB, outputFolder):
    # Verificar se o diretório existe
    if not os.path.exists(outputFolder):
        raise FileNotFoundError(f"O diretório {outputFolder} não existe.")

    for inputContig in os.listdir(outputFolder):
        if inputContig.endswith("_merged.fasta"):
            sample = inputContig.replace("_merged.fasta", "")

            command_dmnd = f"diamond blastx -d {viralDB} -q {os.path.join(outputFolder, inputContig)} --outfmt 6 qseqid sseqid qlen slen pident evalue qcovhsp stitle full_qseq --max-target-seqs 1 --out {os.path.join(outputFolder, sample + '_diamond.tsv')}"
            
            # Executar o comando diamond
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
            
            # Comando sed para substituir colchetes por tabulações e remover colchetes de fechamento
            sed_com = f"sed 's/\\[/\\t/g; s/\\]//g' {entrada}"
            
            # Cabeçalho para adicionar ao arquivo processado
            newcols = "QuerySeq\tSubjectSeq\tQseqLength\tSseqLength\tPident\tEvalue\tQCover\tSubjTitle\tSpecies\tFullQueryLength\n"
            
            try:
                # Executar comando sed e capturar a saída
                result = subprocess.run(sed_com, shell=True, check=True, capture_output=True, text=True)
                
                # Escrever o cabeçalho e a saída do sed no arquivo de saída
                with open(saida, 'w') as out_file:
                    out_file.write(newcols)
                    out_file.write(result.stdout)
                
                # Remover a tabela original
                os.remove(entrada)
            except subprocess.CalledProcessError as e:
                print(f"Erro ao processar o arquivo {entrada}: {e}")
            except Exception as e:
                print(f"Erro inesperado ao processar o arquivo {entrada}: {e}")

