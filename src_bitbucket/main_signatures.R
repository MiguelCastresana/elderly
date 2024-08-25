sysname = Sys.getenv("USERNAME")

if(sysname == "migcas") {
  setwd('//wsl.localhost/Ubuntu/home/miguecas/github_miguel/elderly/') 
  place = "migcas"
} else if(sysname == "Miguel") {
  setwd("C:/Users/Miguel/OneDrive - Karolinska Institutet/Desktop/Postdoc/") 
  place = "Miguel"
}


# Package loading
require(dplyr)
require(stringi)
require(stringr)
require(genefu)
require(MetaGxBreast)
require(ggplot2)
require(genefu)
require(data.table)
require(readxl)

require(survival)
require(survminer)
require(lubridate)


'%!in%' <- function(x,y)!('%in%'(x,y))


###############################       data loading              ###############################
###############################       data annotation           ###############################
###############################           data formatting         ###############################
###############################                                  ###############################


# Load annotation

load("data_bitbucket/esets_1_37")
load("data_bitbucket/esets_38_39")



# all annotation
total_annotation = read.delim("data_bitbucket/total_annotation_all_genes")
source("src_bitbucket//annotations.R")

# curate datasets due to format requirements
# source("src_bitbucket/data_curation.R")

data_clean = readRDS("data_bitbucket/dataclean.rds")

groups = names(data_clean)




###############################       Oncotype DX / 21-gene      ###############################
###############################                                  ###############################
###############################                                  ###############################
###############################                                  ###############################


# Load complete mapping

total_annotation = read.delim("data_bitbucket/oncotype_dx/total_annotation_oncotype")

# LOAD ORIGINAL oncotype DX
data("sig.oncotypedx")

NCBI70 = sig.oncotypedx

# Load missing dataset annotation
load("data_bitbucket/oncotype_dx/gpl_annotation_missing_datasets_oncotype")



# Run Oncotype: Output named "oncotype_results"
source("src_bitbucket/oncotype_function.R")

# 
# 
# 
# 
lista = list()
i = 1
for(i in 1:length(oncotype_results)){
  
  
  
  
  lista[[i]] = oncotype_results[[i]]$result
}




###############################       Mammaprint / 70-gene       ###############################
###############################                                  ###############################
###############################                                  ###############################
###############################                                  ###############################


#load original annotation file
load(paste("data_bitbucket/70_gene/NKIAnnotation.RData",sep=""))
NCBI70<-read.table("data_bitbucket/70_gene/mammaprint_nick.txt",sep="\t",header=T)

# 70-gene correlations
correlations_70gene = read.delim("data_bitbucket/70_gene/correlationsnadel.txt")


# Load complete mapping
total_annotation = read.delim("data_bitbucket/70_gene/total_annotation_mammaprint")

# Load missing dataset annotation
load("data_bitbucket/70_gene/gpl_annotation_missing_datasets_mammaprint")


# Run Oncotype: Output named "mammaprint_results"
source("src_bitbucket/mammaprint_function.R")


lista = list()
i = 1
for(i in 1:length(mammaprint_results)){
  
  
  lista[[i]] = mammaprint_results[[i]]$result
}



###############################       GGI            ###############################
###############################                      ###############################
###############################                      ###############################
###############################                      ###############################


# Load complete mapping
total_annotation = read.delim("data_bitbucket/ggi/total_annotation_ggi")



# LOAD ORIGINAL ggi
sig.ggi = read.delim("data_bitbucket/ggi/ggi_signature_complete")

NCBI70 = sig.ggi



# Load missing dataset annotation

load("data_bitbucket/ggi/gpl_annotation_missing_datasets_ggi")

# Load annotation
source("src_bitbucket/annotations.R")

# Run Oncotype: Output named "ggi_results"
source("src_bitbucket/ggi_function.R")

ggi_genefu = ggi_results

lista = list()
i = 1
for(i in 1:length(ggi_genefu)){
  
  col1 = names(ggi_genefu[[i]]$risk)
  col2 = as.numeric(as.vector(ggi_genefu[[i]]$risk))
  
  dat = as.data.frame(cbind(col1,col2))
  lista[[i]] = dat
}



###############################       CCS            ###############################
###############################                      ###############################
###############################                      ###############################
###############################                      ###############################



# Load CCS genes
ccs_genes = read_excel("data_bitbucket/Cell_cycle/CCS_genes.xlsx")


# Load complete mapping

total_annotation = read.delim("data_bitbucket/total_annotation_all_genes")



# Run Cell cycle: Output named "cell_cycle_results"
source("src_bitbucket/cell_cycle_function.R")




###############################       PAM50            ###############################
###############################                      ###############################
###############################        ROR-P              ###############################
###############################                      ###############################








# Load complete mapping

total_annotation = read.delim("data_bitbucket/PAM50_ROR_P//total_annotation_pam50")

# total_annotation$gene[which(total_annotation$gene == "PALM2AKAP2")] <- "PALM2-AKAP2"



data("pam50")

NCBI70<-pam50[[7]]



# Load missing dataset annotation

load("data_bitbucket/PAM50_ROR_P/gpl_annotation_missing_datasets_pam50")






paramDir <- paste("data_bitbucket/PAM50_ROR_P/bioclassifier_R",sep = "") # the location of unchanging files such as the function library and main program
inputDir <- paste("data_bitbucket/PAM50_ROR_P/",sep = "")  # the location of the data matrix, and where output will be located

inputFile <- "MedCenForPAM50_Genes.txt" # the input data matrix as a tab delimited text file
short <- "MedCenForPAM50_Genes" # short name that will be used for output files


calibrationParameters <- -1 # the column of the "mediansPerDataset.txt" file to use for calibration; 
# NA will force centering within the test set & -1 will not do any 
# adjustment (when adjustment performed by used)

hasClinical <- FALSE  # may include tumor size as second row, with 'T' as the gene name, 
# and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
# set this variable to FALSE if tumor size is not available

collapseMethod <- "mean"  # can be mean or iqr (probe with max iqr is selected)
# typically, mean is preferred for long oligo and
# iqr is preferred for short oligo platforms
# (not applied if dataet is already collapsed, i.e. 1 probe per gene)


source("src_bitbucket/PAM50_RORP_function.R")





###############################       Merge results            ###############################
###############################                      ###############################
###############################                     ###############################
###############################                      ###############################





# Load gene signature results
load("data_bitbucket/final_results/mammaprint_results")
load("data_bitbucket/final_results/ggi_results")
load("data_bitbucket/final_results/oncotype_results")
load("data_bitbucket/final_results/cell_cycle_results")


source("src_bitbucket/final_analysis_all_above70.R")





###############################      Kaplan Meier           ###############################
###############################                      ###############################
###############################                     ###############################
###############################                      ###############################



source("src_bitbucket/survival_analysis.R")








###############################      55-65 analysis         ###############################
###############################                      ###############################
###############################                     ###############################
###############################                      ###############################





# Load gene signature results
load("data_bitbucket/final_results/mammaprint_results")
load("data_bitbucket/final_results/ggi_results")
load("data_bitbucket/final_results/oncotype_results")
load("data_bitbucket/final_results/cell_cycle_results")


source("C:/Users/migcas/OneDrive - Karolinska Institutet/Desktop/Postdoc/src_bitbucket/final_analysis_all_55_65.R")





###############################      Kaplan Meier           ###############################
###############################                      ###############################




source("src_bitbucket/survival_analysis.R")





