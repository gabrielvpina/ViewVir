#!/bin/bash

# Criar um novo diretório
dir_name="diamond-output"
mkdir -p "$dir_name"

# Obter a lista de arquivos na pasta STAR_unmapped/
lista_arquivos=$(ls STAR_unmapped/)

# Iterar sobre os arquivos na lista
for arquivo in $lista_arquivos; do
    # Caminho de entrada e saída para o Fastp
    entrada="STAR_unmapped/$arquivo"
    
    saida="$dir_name/$arquivo"

./diamond blastx -d viralDB.dmnd -q "$entrada" --outfmt '6' qseqid sseqid qlen slen pident evalue bitscore stitle full_qseq --max-target-seqs '1' --out "$saida".tsv 
    
done



