#!/bin/bash

# Criar um novo diretório
dir_name="cap3-output"
mkdir -p "$dir_name"

# Obter a lista de arquivos na pasta STAR_unmapped/
lista_arquivos=$(ls "nonHost_sequences/")

# Iterar sobre os arquivos na lista
for arquivo in $lista_arquivos; do
    # Caminho de entrada e saída para o Fastp
    entrada="nonHost_sequences/$arquivo"
    
    saida="$dir_name/$arquivo"

cap3 $entrada >> $saida
    
done

cd cap3-output 

rm *[^.fasta]

cd ..


