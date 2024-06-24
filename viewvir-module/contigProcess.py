import os
import subprocess

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
    
#def diamond(): 

#def processDmndOut():

