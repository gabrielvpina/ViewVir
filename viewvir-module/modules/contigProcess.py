import os,subprocess

def cap3(inputContig,outputFolder):

    cap3 = ["cap3", inputContig]
    subprocess.run(cap3, check=True)

    for contigs in os.listdir(outputFolder):
        if contigs.endswith(.cap.contigs):
            sample = os.path.splitext(contigs)[0]
            print(f"Processando amostra {sample}")

            singlets = f"{sample}.cap.singlets"

            if os.path.exists(singlets):
                command = f"cat {contigs} {singlets} > {os.path.join(outputFolder, f'{sample}_merged.fasta')}"
                # O resultado é o arquivo sample_merged.fasta + o arquivo .fasta original
    
    log_files = [os.path.join(outputFolder, file) for file in os.listdir(outputFolder) if file.endswith(".links",".qual",".info",".contigs",".singlets",".ace")]
    for log_file in log_files:
            os.remove(log_file)

seq = "mycontigs.fasta"
folder = "test"
cap3(seq,folder)

#def diamond(): 

#def processDmndOut():

