#!/bin/bash

dir_name="ORFs"
mkdir -p "ViewVir-results/$dir_name"

lista_arquivos=$(ls "ViewVir-results"/*.fasta)

for arquivo in $lista_arquivos; do

entrada=$arquivo    
saida="ViewVir-results/$dir_name"

nome=$(basename "$arquivo" .fasta)

orfipy --partial-3 --partial-5 --outdir $saida $entrada --bed "${nome}.bed"

done

rm ViewVir-results/ORFs/*.log
