############################################################ 
#         Examples for how to run the Metabolizer          #
#                                                          #
#         cankutcubuk [at] {gmail} [dot] {com}             #
#                     2016-2019                            #
#             @ Sevilla and Valencia, Spain                #
############################################################ 



##################################################
#             Example of Homo sapiens            #
##################################################


rm(list=ls())

#library(devtools)
source("https://raw.githubusercontent.com/babelomics/metabolizer/develop/metabolizer.git/script/utils/functions.main.14032019.R")
load("./metabolizer.git/data//KEGGfiles/hsa/hsa_module_data_March2019.RData")

### just to create a fake gene expression file: an example input file ###
# KEGG or NCBI (Entrez) gene ids are used in the Metabolizer tool. Check which one is used for your model organism
# You can use xxx_geneIDs.RData file to convert your gene symbols into gene ids that are required by the Metabolizer.
  all_module_genes_vec <- get.All.module.genes(hsa_module_data)
# simulate some missing genes
  all_module_genes_vec <- all_module_genes_vec[-sample(1:length(all_module_genes_vec),40)] 
  fakeExp <- matrix(data=rexp(length(all_module_genes_vec)*20, rate = 0.1), nrow = length(all_module_genes_vec), ncol = 20)
  dimnames(fakeExp) <- list(all_module_genes_vec, LETTERS[1:20])
### end of fake gene expression file ###

results <- metabolizer(hsa_module_data, 
                       combat.vals=fakeExp, 
                       output_folder=NULL, saveName=NULL, 
                       onesample=F, 
                       expbased=T, fluxbased=F, 
                       default_value=0.5, 
                       moduleinfo="./metabolizer.git/data/KEGGfiles/hsa/hsa_moduleinfo.RData")


##################################################
# Example of Drosophila melanogaster (fruit fly) #
##################################################

rm(list=ls())

source("./metabolizer.git/script/utils/functions.main.14032019.R")
load("./metabolizer.git/data//KEGGfiles/dme/dme_module_data_March2019.RData")

### just to create a fake gene expression file: an example input file ###
# KEGG or NCBI (Entrez) gene ids are used in the Metabolizer tool. Check which one is used for your model organism
# You can use xxx_geneIDs.RData file to convert your gene symbols into gene ids that are required by the Metabolizer.
all_module_genes_vec <- get.All.module.genes(hsa_module_data)
# simulate some missing genes
all_module_genes_vec <- all_module_genes_vec[-sample(1:length(all_module_genes_vec),40)] 
fakeExp <- matrix(data=rexp(length(all_module_genes_vec)*20, rate = 0.1), nrow = length(all_module_genes_vec), ncol = 20)
dimnames(fakeExp) <- list(all_module_genes_vec, LETTERS[1:20])
### end of fake gene expression file ###

results <- metabolizer(hsa_module_data, 
                       combat.vals=fakeExp, 
                       output_folder=NULL, saveName=NULL, 
                       onesample=F, 
                       expbased=T, fluxbased=F, 
                       default_value=0.5, 
                       moduleinfo="./metabolizer.git/data/KEGGfiles/dme/dme_moduleinfo.RData")

