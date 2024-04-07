#!/bin/bash

# Criar um novo diretório
dir_name="diamond-processed"
mkdir -p "$dir_name"

# Obter a lista de arquivos na pasta diamond-output/
lista_arquivos=$(ls diamond-output/)

# Iterar sobre os arquivos na lista
for arquivo in $lista_arquivos; do
    # Caminho de entrada e saída para o Fastp
    entrada="diamond-output/$arquivo"
    saida="$dir_name/$arquivo"


sed -i 's/\[/\t/g; s/\]//g' $entrada


echo "QuerySeq	SubjectSeq	QseqLength	SseqLength	Pident	Evalue	BitScore	SubjTitle	Species	FullQueryLength" | cat - $entrada > $saida 
# && mv $saida $entrada

done 
