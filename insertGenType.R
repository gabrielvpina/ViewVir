# Insert the directory of diamond-processed file
WorkingDir = "/home/gabriel/bioinfo/bioinfo_scripts/scripts_viral"
setwd(WorkingDir)

install.packages("dplyr")

library(dplyr)

# import all tables from ICTV
ictvMSL <- read.csv("ICTV_data/ICTV_MSL.csv")
colnames(ictvMSL) <- c("Species","Genome.composition")

ictvVMR_species <- read.csv("ICTV_data/ICTV_VMR-species.csv")
ictvVMR_virNames <- read.csv("ICTV_data/ICTV_VMR-virNames.csv")

# merge all tables with rbind
allVirus <- rbind(ictvMSL, ictvVMR_species, ictvVMR_virNames)

arquivos <- zlist.files(path = "/home/gabriel/bioinfo/bioinfo_scripts/scripts_viral/diamond-processed", pattern = "\\.tsv$", full.names = TRUE)

# Crie a pasta para os resultados
dir.create("diamond-results")

# Loop para cada tabela .tsv
for (arquivo_path in arquivos) {
  # Leia o arquivo
  arquivo <- read.table(arquivo_path, header = TRUE, sep = "\t")
  
  # Realize a operação de junção
  arquivo <- left_join(arquivo, allVirus, by = c("Species" = "Species"))
  
  # Mude a ordem das colunas
  outrasCols <- setdiff(names(arquivo), "Genome.composition")
  arquivo <- arquivo[, c("QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "BitScore", "SubjTitle",
                         "Species", "Genome.composition", "FullQueryLength")]
  
  # Obtenha o nome do arquivo de entrada
  nome_arquivo_input <- basename(arquivo_path)
  
  # Defina o caminho completo do arquivo de saída, incluindo o nome da pasta
  caminho_arquivo_output <- file.path("diamond-results", sub(".tsv", "_output.tsv", nome_arquivo_input))
  
  # Escreva os dados em um arquivo de saída
  write.table(arquivo, file = caminho_arquivo_output, sep = "\t")

}






















}




