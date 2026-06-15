#########################################################
## PAM50 subtyping preprocessing script
## Author: Miguel Castresana-Aguirre / Annelie Johansson / Christina Yau
## Date: 20th April 2022
## 1. Calculates median of medians of all probes in a 
##    random sample of 50% ER+ and 50% ER- patients
## 2. Center data on median of medians
## 3. Save data in correct format for PAM50 subtyping
#########################################################






#####        TIME FOR THE ANALYSIS ####################################



i = 1
for(i in 1:length(data_clean)){
  
  ## Load data ----
  # setwd("")
  # - normalized gene expression data (recommend collapsed data, i.e. 1 probe per gene)
  # - matching patient data (with information on ER-status)
  # - annotation file with ProbeIDs and GeneNames
  if(i==17){
    next
  }
  GEX_data <- data_clean[[i]][[1]]
  PAT_data <- data_clean[[i]][[3]]
  annot <- data_clean[[i]][[2]]
  
  if("ORC6"%in%rownames(GEX_data) & "ORC6L"%in%rownames(GEX_data)){
    
    GEX_data = GEX_data[which(rownames(GEX_data)%!in%"ORC6"),]
    
  }
  if("NUF2"%in%rownames(GEX_data) & "CDCA1"%in%rownames(GEX_data)){
    
    GEX_data = GEX_data[which(rownames(GEX_data)%!in%"NUF2"),]
    
  }
  if("NDC80"%in%rownames(GEX_data) & "KNTC2"%in%rownames(GEX_data)){
    
    GEX_data = GEX_data[which(rownames(GEX_data)%!in%"NDC80"),]
    
  }
  
  # For example, genes can have several different gene names:
  annot$symbol[which(annot$symbol == "ORC6")] <- "ORC6L"
  annot$symbol[which(annot$symbol == "NUF2")] <- "CDCA1"
  annot$symbol[which(annot$symbol == "NDC80")] <- "KNTC2"
  
  
  # some genes here have the same entrez and differnt symbol
  
  if(ncol(annot)<3){
    
    annot = annot[!duplicated(annot[,2]),]
    
  }else{
    
    annot = annot[!duplicated(annot[,3]),]
    
  }
  
  rownames(GEX_data)[which(rownames(GEX_data)=="ORC6")]<- "ORC6L"
  rownames(GEX_data)[which(rownames(GEX_data)=="NUF2")]<- "CDCA1"
  rownames(GEX_data)[which(rownames(GEX_data)=="NDC80")]<- "KNTC2"
  
  GEX_data = GEX_data[which(rownames(GEX_data)%in%annot$symbol),]
  
  
  
  
  
  # ## Get genes used in the PAM50 subtyping 
  pam50_genes <- NCBI70$probe
  
  # 
  colnames(GEX_data) = PAT_data$sample_name
  
  PAT_data1 <- PAT_data[complete.cases(PAT_data),]
  
  
  if(nrow(PAT_data1)<1){
    next
  }
  
  
  GEX_data <- GEX_data[,which(colnames(GEX_data)%in%PAT_data$sample_name)]
  
  # Step 1: Get the order of PAT_data's SampleID based on the column names of GEX_data
  order_index <- match(PAT_data$sample_name, colnames(GEX_data))
  
  # Step 2: Reorder PAT_data based on the obtained order
  PAT_data <- PAT_data[order_index, ]
  
  # PAT_data = PAT_data[!is.na(PAT_data$er),]
  
  GEX_data = GEX_data[,which(colnames(GEX_data)%in%PAT_data$sample_name)]
  
  
  ## Extract the 50 genes form your data
  # annot <- annot[which(annot$symbol%in%pam50_genes),]
  # GEX_data <- GEX_data[which(rownames(GEX_data)%in%annot$symbol),]
  
  ## Function to calculate median of all probes ----
  # Function to random sub-sample ER+ cases to make 50% and 50% composition
  # Output: Median of dataset for centering
  
  # Lets set the seeds so this can be replicated:
  
  
  
  
  library(parallel)
  
  # Example data (Replace with your actual data)
  # GEX_data <- matrix(rnorm(1000), ncol = 10) # Example gene expression data
  # PAT_data <- data.frame(er = c(rep("positive", 5), rep("negative", 5))) # Example patient data
  
  # The makemed function adjusted for direct access to GEX_data and PAT_data
  makemed <- function() {
    neg <- sum(PAT_data$er == "negative", na.rm = TRUE)
    pos <- sum(PAT_data$er == "positive", na.rm = TRUE)
    
    indneg <- which(PAT_data$er == "negative")
    indpos <- which(PAT_data$er == "positive")
    
    if (neg >= pos) {
      nb <- pos
      wantneg <- sample(indneg, nb, replace = FALSE)
      temp <- cbind(GEX_data[, wantneg], GEX_data[, indpos])
    } else {
      nb <- neg
      wantpos <- sample(indpos, nb, replace = FALSE)
      temp <- cbind(GEX_data[, indneg], GEX_data[, wantpos])
    }
    
    out <- apply(temp, 1, function(x) { median(x, na.rm = TRUE) })
    return(out)
  }
  
  
  
  # Number of Monte Carlo simulations
  n_simulations <- 100
  
  # Detect number of cores, reserving one for system stability
  no_cores <- 4
  
  # Run simulations in parallel
  # Initialize a cluster
  cl <- makeCluster(no_cores)
  
  # Set a single seed for reproducibility
  clusterSetRNGStream(cl, 1556)
  
  # Export the makemed function and data to the cluster
  clusterExport(cl, varlist = c("makemed", "GEX_data", "PAT_data"))
  
  # Use parLapply to perform the simulations
  Medians <- parLapply(cl, X = 1:n_simulations, fun = function(x) makemed())
  
  # Stop the cluster
  stopCluster(cl)
  
  # Convert list of medians to a matrix (if necessary)
  Medians_matrix <- do.call(cbind, Medians)
  
  
  
  Medians <- as.data.frame(Medians_matrix)
  
  # Optionally, set the column names to reflect the simulation number
  colnames(Medians) <- paste0("Simulation_", 1:n_simulations)
  
  # Ensure the row names match those of GEX_data, if GEX_data has meaningful row names
  rownames(Medians) <- rownames(GEX_data)
  
  
  
  ## Calculate average of all medians
  # av_med <- apply(Medians, 1, mean)
  # summary(av_med)
  
  
  ## Calculate median of all medians
  med_med <- apply(Medians, 1, median)
  summary(med_med)
  
  
  ## The difference between using the mean or median is neglibile ##
  ## Going with the theme - let's do median of medians ##
  
  ## Save medians of medians 
  # write.table(med_med, "ER_EqualSplit_MedianOfMedians.txt", sep="\t", row.names=T, quote=F)
  # med_med <- read.table("ER_EqualSplit_MedianOfMedians.txt", sep = "\t", stringsAsFactors = FALSE)
  
  #### 2. Center your data to median of medians ----
  GEX_data_med_med <- sweep(GEX_data, 1, med_med)
  
  ## Save centered gene expression data
  # save(GEX_data_med_med, file = "MedCenForPAM50.RData")
  # load("MedCenForPAM50.RData")
  
  #### 3. Get data in right format for PAM50 subtyping ----
  
  
  
  if(length(annot$symbol[which(annot$symbol%!in%rownames(GEX_data))])<1){
    annot = annot[match(rownames(GEX_data_med_med), annot$symbol),]
  }
  
  
  ## make sure identical ProbeNames in gene expression data and annotation data
  identical(as.vector(annot$symbol), rownames(GEX_data_med_med)) # TRUE
  
  ## Get genes used in the PAM50 subtyping 
  
  ## If less than 50 genes are matching to your data: explore why?
  # For example, here only 48 genes are matching to our data
  length(which(pam50_genes %in% annot$symbol)) # 48
  # pam50_genes[which(pam50_genes %!in% annot$symbol)] # 48
  
  
  
  
  
  
  ## Select gene expression data with only PAM50 genes
  annot_pam50 <- annot
  GEX_pam50 <- GEX_data_med_med
  
  
  
  ## Extract the 50 genes form your data
  annot_pam50 <- annot_pam50[which(annot_pam50$symbol%in%pam50_genes),]
  GEX_pam50 <- GEX_pam50[which(rownames(GEX_pam50)%in%annot_pam50$symbol),]
  # 
  identical(rownames(GEX_pam50), as.vector(annot_pam50$symbol)) # TRUE
  
  ## Make sure we have the data in the correct format for PAM50 subtyping ##
  
  ## Sample names first row
  GEX_pam50 <- rbind(colnames(GEX_pam50), GEX_pam50)
  
  ## Gene names first column
  GEX_pam50 <- cbind(GeneName = c("X", as.vector(annot_pam50$symbol)), GEX_pam50)
  
  
  dir.create("data_bitbucket/PAM50_ROR_P/PAM50_results/")
  
  dir.create(paste("data_bitbucket/PAM50_ROR_P/PAM50_results/",groups[i],sep = ""))
  
  
  ## Save
  write.table(GEX_pam50, paste("data_bitbucket/PAM50_ROR_P/PAM50_results/",
                               groups[i],"/MedCenForPAM50_Genes.txt",sep = ""), sep="\t", row.names=F, col.names = F, quote=F)
  
  
  
  paramDir <- paste("data_bitbucket/PAM50_ROR_P/bioclassifier_R",sep = "") # the location of unchanging files such as the function library and main program
  inputDir <- paste("data_bitbucket/PAM50_ROR_P/PAM50_results/",groups[i],sep = "")  # the location of the data matrix, and where output will be located
  
  inputFile <- "MedCenForPAM50_Genes.txt" # the input data matrix as a tab delimited text file
  short <- "MedCenForPAM50_Genes" # short name that will be used for output files
  
  
  ####
  # run the assignment algorithm
  ####
  
  # load functions 
  source(paste(paramDir,"subtypePrediction_functions.R",sep="/")) # original file
  
  # run scripts
  
  source(paste(paramDir,"subtypePrediction_distributed_AJ.R",sep="/")) # I have fixed script a little from the original file!
  
  print(i)
  
  # if(i==1){break}
  
}

