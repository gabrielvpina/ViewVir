#!/bin/bash

# Criar um novo diretório
dir_name="diamond-processed"
mkdir -p "$dir_name"

# Obter a lista de arquivos na pasta diamond-output/
lista_arquivos=$(ls diamond-output/)

# Iterar sobre os arquivos na lista
for arquivo in $lista_arquivos; do

    # Remover o sufixo ".fasta" do nome do arquivo de entrada
    nome_sem_sufixo=$(echo "$arquivo" | sed 's/\.fasta$//')
    
    # Caminho de entrada e saída para o Fastp
    entrada="diamond-output/$arquivo"
    saida="$dir_name/$nome_sem_sufixo"


    sed -i 's/\[/\t/g; s/\]//g' $entrada


echo "QuerySeq	SubjectSeq	QseqLength	SseqLength	Pident	Evalue	QCover	SubjTitle	Species	FullQueryLength" | cat - $entrada > $saida 
# && mv $saida $entrada

done 

rename 's/.fasta//' diamond-processed/*

rm -r diamond-output
