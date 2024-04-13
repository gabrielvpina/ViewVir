#!/bin/bash

# Criar um novo diretório
dir_name="cap3-processed"
mkdir -p "$dir_name"

for arquivo in "cap3-output"/*.cap.contigs; do
   
    if [[ -f "$arquivo" ]]; then
       # Extrair nome do arquivo
       samp=$(basename "$arquivo" .cap.contigs)

       echo "Processing sample $samp"

           
           arquivo_2="cap3-output/${samp}.cap.singlets"

           if [[ -f "$arquivo_2" ]]; then

                 cat "$arquivo" "$arquivo_2" >> "$dir_name/${samp}_merged.fasta" 
        
            fi

    fi   

done    

rm -r cap3-output
