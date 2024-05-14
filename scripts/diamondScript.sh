#!/bin/bash

# Criar um novo diretório
c

# Obter a lista de arquivos na pasta STAR_unmapped/
lista_arquivos=$(ls "cap3-processed/")

# Iterar sobre os arquivos na lista
for arquivo in $lista_arquivos; do
    # Caminho de entrada e saída para o Fastp
    entrada="cap3-processed/$arquivo"
    
    saida="$dir_name/$arquivo"

diamond blastx -d viralDB.dmnd -q "$entrada" --outfmt '6' qseqid sseqid qlen slen pident evalue qcovhsp stitle full_qseq --max-target-seqs 1 --out "$saida".tsv 
    
done

rm -r cap3-processed



