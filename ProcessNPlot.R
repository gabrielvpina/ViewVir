if (!require("here", quietly = TRUE))
install.packages("here")
library(here)

# Insert the directory of diamond-processed file
WorkingDir = here()
setwd(WorkingDir)

# Required packages
if (!require("dplyr", quietly = TRUE))
install.packages("dplyr")

if (!require("pandoc", quietly = TRUE))
install.packages("pandoc")

if (!require("plotly", quietly = TRUE))
install.packages("plotly")

if (!require("gapminder", quietly = TRUE))
install.packages("gapminder")

if (!require("htmltools", quietly = TRUE))
install.packages("htmltools")

if (!require("ggplot2", quietly = TRUE))
install.packages("ggplot2")

library(dplyr)
library(ggplot2)
library(plotly)
library(gapminder)
library(htmlwidgets)
library(pandoc)



# import all tables from ICTV
ictvMSL <- read.csv("ICTV_data/ICTV_MSL.csv")
colnames(ictvMSL) <- c("Species","Genome.composition")

ictvVMR_species <- read.csv("ICTV_data/ICTV_VMR-species.csv")
ictvVMR_virNames <- read.csv("ICTV_data/ICTV_VMR-virNames.csv")

# merge all tables with rbind
allVirus <- rbind(ictvMSL, ictvVMR_species, ictvVMR_virNames)

arquivos <- list.files(path = here("diamond-processed"), pattern = "\\.tsv$", full.names = TRUE)

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
  arquivo <- arquivo[, c("QuerySeq", "SubjectSeq", "QseqLength", "SseqLength", "Pident", "Evalue", "BitScore", "SubjTitle",
                         "Species", "Genome.composition", "FullQueryLength")]

  
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
  plot_output <- file.path(pasta_nome_arquivo, paste0(sub(".tsv", "_bubblePlt.html", nome_arquivo_input)))
  
  # Escreva os dados em um arquivo de saída dentro do diretório criado
  write.table(arquivo, file = caminho_arquivo_output, sep = "\t")
##################################################################################
  
# Grafico Interativo
  
  plot <- arquivo %>%
    ggplot( aes(QseqLength, Species, size = Pident, color=Genome.composition)) +
    geom_point() +
    scale_y_discrete(labels = function(x) rep("", length(x))) +
    theme_bw()
  
  ggplotly(plot)
  
  plot <- plotly::ggplotly(plot)
  
  htmlwidgets::saveWidget(plot, file = plot_output, selfcontained = TRUE)

}

