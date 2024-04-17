
# Required packages
if (!require("dplyr", quietly = TRUE))
install.packages("dplyr")

if (!require("plotly", quietly = TRUE))
install.packages("plotly")

if (!require("gapminder", quietly = TRUE))
install.packages("gapminder")

if (!require("htmltools", quietly = TRUE))
install.packages("htmltools")

if (!require("ggplot2", quietly = TRUE))
install.packages("ggplot2")

if (!require("DT", quietly = TRUE))
  install.packages("DT")

library(dplyr)
library(ggplot2)
library(plotly)
library(gapminder)
library(htmlwidgets)
library(DT)


# import all tables from ICTV
ictvMSL <- read.csv("data/ICTV_MSL.csv")
colnames(ictvMSL) <- c("Species","Genome.composition")

ictvVMR_species <- read.csv("data/ICTV_VMR-species.csv")

ictvVMR_virNames <- read.csv("data/ICTV_VMR-virNames.csv")

ncbiSpecie <- read.csv("data/NCBI_virSpecies.csv")
colnames(ncbiSpecie) <- c("Species","Genome.composition")

ncbiNames <- read.csv("data/NCBI_virName.csv")
colnames(ncbiNames) <- c("Species","Genome.composition")

# merge all tables with rbind
allVirus <- rbind(ncbiNames, ncbiSpecie, ictvMSL, ictvVMR_species, ictvVMR_virNames)

arquivos <- list.files(path = "diamond-processed/", pattern = "\\.tsv$", full.names = TRUE)

# Crie a pasta para os resultados
if (!dir.exists("ViewVir-results")) {
  dir.create("ViewVir-results")
}

# Loop para cada tabela .tsv
for (arquivo_path in arquivos) {

# Processando a Tabela    
  # Leia o arquivo
  arquivo <- read.table(arquivo_path, header = TRUE, sep = "\t")
  
  # Realize a operação de junção
  arquivo <- left_join(arquivo, allVirus, by = c("Species" = "Species"))
  
  # Mude a ordem das colunas
  outrasCols <- setdiff(names(arquivo), "Genome.composition")
  arquivo <- arquivo[, c("QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "QCover", "SubjTitle",
                         "Species", "Genome.composition", "FullQueryLength")]
  
  # Filtro de sequencias
  arquivo <- arquivo[arquivo$QseqLength >= 450,]
  
  arquivo  <- dplyr::distinct(arquivo)
  
  
  

############################## OLD VERSION OUTPUT ###########################################  
  # Obtenha o nome do arquivo de entrada
  #nome_arquivo_input <- basename(arquivo_path)
  
  
  # Defina o caminho completo do arquivo de saída, incluindo o nome da pasta
  #caminho_arquivo_output <- file.path("ViewVir-results/", sub(".tsv", "_output.tsv", nome_arquivo_input))
  #plot_output <- file.path("ViewVir-results", sub(".tsv","_bubblePlt.html", nome_arquivo_input))
  
  # Escreva os dados em um arquivo de saída
  #write.table(arquivo, file = caminho_arquivo_output, sep = "\t")
################################################################################  
  
  
#################### função subpastas (teste) ###################################  
  # Obtenha o nome do arquivo de entrada
  nome_arquivo_input <- basename(arquivo_path)
  
  # Crie o diretório com o nome do arquivo dentro de "diamond-results"
  pasta_nome_arquivo <- file.path("ViewVir-results", gsub(".fasta.tsv", "", nome_arquivo_input))
  
  dir.create(pasta_nome_arquivo, showWarnings = FALSE)
  
  # Defina o caminho completo do arquivo de saída dentro do diretório criado
  caminho_arquivo_output <- file.path(pasta_nome_arquivo, paste0(sub(".tsv", "_output.tsv", nome_arquivo_input)))
  
  plot_output_bubble <- file.path(pasta_nome_arquivo, paste0(sub(".tsv", "_bubblePlt.html", nome_arquivo_input)))
  plot_output_table <- file.path(pasta_nome_arquivo, paste0(sub(".tsv", "_table.html", nome_arquivo_input)))
  
  
  # Escreva os dados em um arquivo de saída dentro do diretório criado
  write.table(arquivo, file = caminho_arquivo_output, sep = "\t", row.names = FALSE)
  
  
################################### RNA Virus Table ##################################################################  
  
  arquivoRNA <- arquivo[arquivo$Genome.composition %in% c("ssRNA(-)","ssRNA(+/-)","ssRNA(+)","ssRNA-RT","RNA"),]
  
  arquivoRNA  <- dplyr::distinct(arquivoRNA)
  
  output_RNA_table <- file.path(pasta_nome_arquivo, paste0(sub(".tsv", "_RNA-virus.tsv", nome_arquivo_input)))
  
  write.table(arquivoRNA, file = output_RNA_table, sep = "\t", row.names = FALSE)
  
  arquivoPRE <- arquivo[arquivo$Genome.composition %in% c("ssRNA(-)","ssRNA(+/-)","ssRNA(+)","ssRNA-RT","RNA"),]
  arquivoPRE <- arquivoPRE[,c(1,11)]
 
  
  output_pre <- file.path(pasta_nome_arquivo, paste0(sub(".tsv", "_pre_fasta.tsv", nome_arquivo_input)))
  
  write.table(arquivoPRE, file = output_pre, sep = "\t", row.names = FALSE)
  
  
###################################### PRE FASTA TABLE ##################################################################
  
  
  
  
  


  
# Grafico Interativo
  
  arquivo$MatchSequence <- paste(arquivo$Species, arquivo$SubjTitle, sep = ' - ')
  
  plot <- arquivo %>%
    ggplot(aes(QseqLength, MatchSequence, size = QCover, fill=Pident, color=Genome.composition)) +
    geom_point() +
    scale_y_discrete(labels = function(x) rep("", length(x))) +
    theme_bw()
  
  ggplotly(plot)
  
  plot <- plotly::ggplotly(plot)
  
  table <- datatable(arquivo)
  
  htmlwidgets::saveWidget(plot, file = plot_output_bubble, selfcontained = FALSE)
  htmlwidgets::saveWidget(table, file = plot_output_table, selfcontained = FALSE)

}
















