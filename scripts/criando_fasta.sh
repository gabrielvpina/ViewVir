#!/bin/bash

lista_arquivos=$(ls "ViewVir-results"/*/*pre_fasta.tsv)

for arquivo in $lista_arquivos; do

samp=$(basename "$arquivo" _pre_fasta.tsv)

   # Removendo 1 linha | inserindo ">" em todas as linhas | substituir tab por \n
  cat $arquivo | sed '1d' | sed 's/^/>/' | sed 's/\t/\n/' | sed 's/"//' | sed 's/"//' >> "ViewVir-results/${samp}_RNA-virus.fasta"
   

done
   
   
   




