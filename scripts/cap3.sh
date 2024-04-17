#!/bin/bash

# Criar um novo diretório
dir_name="cap3-output"
mkdir -p "$dir_name"

cp nonHost_sequences/* cap3-output/

# Obter a lista de arquivos na pasta STAR_unmapped/
lista_arquivos=$(ls "cap3-output/")

# Iterar sobre os arquivos na lista
for arquivo in $lista_arquivos; do

    cap3 "cap3-output/$arquivo"

done


