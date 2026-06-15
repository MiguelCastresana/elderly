

merge_function = function(){
  
  
  
  
  mammaprint_adaptative_newmap = mammaprint_results
  ggi_genefu = ggi_results
  # Load summary of datasets
  load("data_bitbucket/SUMMARY_datasets")
  
  load("data_bitbucket/datasets_samples")
  
  
  # Summary of survival endpoints per dataset
  g = final_summary$dataset
  
  

  
  
  
  
  
  
  # Filter samples with mappings lower than 75% for Mammaprint and that dont have all non reference oncotype genes and that could be ran in PAM50
  
  sub = final_summary[which(final_summary$ratio_mapped_mammaprint_miguel>0.75 & final_summary$mapped_oncotype_miguel==1 & final_summary$mapped_genes_pam50>0),] # 23 datasets
  
  sum(sub$nb_patients) # 6082
  
  sum(sub$nb_patients_above70,na.rm = TRUE) # 1020
  
  
  
  
  
  
  # Load METABRIC survival data
  
  rfs_metabric = read.delim("data_bitbucket/metabric_survival.txt")
  
  
  
  # Apparently we have treatment data but we need to label it properpy as either: chemo.plus.hormono,chemotherapy,hormonotherapy, untreated
  
  treatment = vector()
  i = 1
  for(i in 1:nrow(rfs_metabric)){
    
    if(rfs_metabric$CT[i]%!in%"NO/NA" & rfs_metabric$RT[i]%!in%"NO/NA" & rfs_metabric$HT[i]%!in%"NO/NA"){
      
      treatment[i] = "chemo.plus.hormono.plus.radio"
    }else if(rfs_metabric$CT[i]%in%"NO/NA" & rfs_metabric$RT[i]%in%"NO/NA" & rfs_metabric$HT[i]%in%"NO/NA"){
      
      treatment[i] = "untreated"
    }else if(rfs_metabric$CT[i]%in%"NO/NA" & rfs_metabric$RT[i]%!in%"NO/NA" & rfs_metabric$HT[i]%!in%"NO/NA"){
      
      treatment[i] = "hormono.plus.radio"
    }else if(rfs_metabric$CT[i]%in%"NO/NA" & rfs_metabric$RT[i]%in%"NO/NA" & rfs_metabric$HT[i]%!in%"NO/NA"){
      
      treatment[i] = "hormonotherapy"
    }else if(rfs_metabric$CT[i]%in%"NO/NA" & rfs_metabric$RT[i]%!in%"NO/NA" & rfs_metabric$HT[i]%in%"NO/NA"){
      
      treatment[i] = "radiotherapy"
    }else if(rfs_metabric$CT[i]%!in%"NO/NA" & rfs_metabric$RT[i]%!in%"NO/NA" & rfs_metabric$HT[i]%in%"NO/NA"){
      
      treatment[i] = "chemo.plus.radio"
    }else if(rfs_metabric$CT[i]%!in%"NO/NA" & rfs_metabric$RT[i]%in%"NO/NA" & rfs_metabric$HT[i]%!in%"NO/NA"){
      
      treatment[i] = "chemo.plus.hormono"
    }else if(rfs_metabric$CT[i]%!in%"NO/NA" & rfs_metabric$RT[i]%in%"NO/NA" & rfs_metabric$HT[i]%in%"NO/NA"){
      
      treatment[i] = "chemotherapy"
    }
    
  }
  
  rfs_metabric$treatment = treatment
  
  
  
  
  
  pheno = rfs_metabric[,c("METABRIC.ID","Age.At.Diagnosis","ER.Status","Lymph.Nodes.Positive","PR.Expr","Her2.Expr","Pam50Subtype","Grade", "Size","T","Death","TLR","LR","TDR","DR","treatment")]
  
  
  # Samples ID change
  pheno$METABRIC.ID = gsub("-", "_", pheno$METABRIC.ID)
  
  
  
  # Adapt RFS data based on TLR (loco regional relapse), TDR (distance relapse), and T (death)
  RFS_days = numeric()
  RFS_status = numeric()
  i = 1
  for(i in 1:nrow(pheno)){
    
    
    if(pheno$LR[i] ==0 & pheno$DR[i] ==0){
      
      RFS_days = c(RFS_days,min(pheno$TLR[i],pheno$TDR[i],pheno$T[i]))
      RFS_status = c(RFS_status,0)
    }else if(pheno$LR[i]==1 | pheno$DR[i] ==1){
      
      RFS_days = c(RFS_days,min(pheno$TLR[i],pheno$TDR[i]))
      RFS_status = c(RFS_status,1)
    }else{
      
      RFS_days = c(RFS_days,NA)
      RFS_status = c(RFS_status,NA)
    }
  }
  
  pheno_1 = as.data.frame(cbind(pheno$METABRIC.ID,RFS_status,RFS_days))
  
  names(pheno_1) = c("METABRIC.ID","dthbcstat2","T_dth")
  
  pheno_1 = merge(pheno,pheno_1,by = "METABRIC.ID")
  
  
  
  pheno_1 = pheno_1[,-c(10:15)]
  
  names(pheno_1) = c("sample_name","age_at_initial_pathologic_diagnosis","er","N","pgr","her2","PAM50_metabric","grade","tumor_size","treatment","dthbcstat2","T_dth")
  
  
  
  pheno_1$dataset = rep("METABRIC",nrow(pheno_1))
  
  
  # From mm to cm
  pheno_1$tumor_size = pheno_1$tumor_size/10
  
  
  
  final_dat = final_dat[which(final_dat$dataset%!in%"METABRIC"),]
  
  
  # Remove known duplicates
  dups = c("MB_0025", "MB_0196", "MB_0326", "MB_0329", "MB_0330", "MB_0335", "MB_0355",
           "MB_0407", "MB_0433", "MB_0547", "MB_2720", "MB_6206")
  
  pheno_1 = pheno_1[which(pheno_1$sample_name%!in%dups),]
  
  
  # Merge both
  final_dat = rbindlist(list(final_dat, pheno_1), fill = TRUE)
  
  
  
  
  # Convert lymph node count into binary >0 positive
  
  final_dat = dplyr::mutate(final_dat,
                            N = ifelse(N >= 1, 1, N))
  
  
  # Lets create binary variables for chemo, hormono and chemo+hormono
  
  i = 1
  chemo_hormono = vector()
  chemo = vector()
  hormono = vector()
  for(i in 1:nrow(final_dat)){
    
    if(final_dat$treatment[i]%in%"chemo.plus.hormono" | final_dat$treatment[i]%in%"chemo.plus.hormono.plus.radio" 
       | final_dat$treatment[i]%in%"chemo.plus.radio" | final_dat$treatment[i]%in%"chemotherapy"){
      
      chemo[i] = "YES"
      hormono[i] = "NO"
      chemo_hormono[i] = "NO"
    }
    if(final_dat$treatment[i]%in%"chemo.plus.hormono" | final_dat$treatment[i]%in%"chemo.plus.hormono.plus.radio" 
       | final_dat$treatment[i]%in%"hormono.plus.radio" | final_dat$treatment[i]%in%"hormonotherapy"){
      
      hormono[i] = "YES"
      chemo[i] = "NO"
      chemo_hormono[i] = "NO"
    }
    if(final_dat$treatment[i]%in%"chemo.plus.hormono" | final_dat$treatment[i]%in%"chemo.plus.hormono.plus.radio"){
      
      chemo_hormono[i] = "YES"
      hormono[i] = "YES"
      chemo[i] = "YES"
    }
    
    if(final_dat$treatment[i]%in%"radiotherapy" | final_dat$treatment[i]%in%"untreated"){
      
      chemo[i] = "NO"
      hormono[i] = "NO"
      chemo_hormono[i] = "NO"
    }
  }
  
  
  final_dat$chemo = chemo
  final_dat$hormono = hormono
  final_dat$chemo_hormono = chemo_hormono
  
  

  
  final_summary[is.na(final_summary)] <- 0
  
  
  
  # Add metabric to the rest
  pos = which(final_summary$dataset%in%"METABRIC")
  vec1 = rep(0,23)
  vec2 = rep(0,15)
  dthbcstat2 = c(vec1,length(final_dat$dthbcstat2[!is.na(final_dat$dthbcstat2)]),vec2)
  
  final_summary[,17] = dthbcstat2
  colnames(final_summary)[17] = "dthbcstat2"
  
  
  
  
  # Select patients above 70, that have ER status, that have Oncotype coverage = 1; mammaprint > 0.75, pam50>0, that have an endpoint that is not days to death. Also we remove duplicates
  # and we add the metabric survival data
  # We will remove the end point days_to_death since it is too general and it may involve commorbidities
  
  
  datasets_to_select = final_summary$dataset[which(final_summary$ratio_mapped_mammaprint_miguel>0.75 & final_summary$mapped_oncotype_miguel==1 & final_summary$mapped_genes_pam50>0)]
  
  
  
  sub = final_dat[which(!is.na(final_dat[,5])  & final_dat$dataset%in%datasets_to_select & final_dat$age_at_initial_pathologic_diagnosis>=70 
                        & (!is.na(final_dat$recurrence_status) | 
                             !is.na(final_dat$dmfs_status) | !is.na(final_dat$dthbcstat2))),]
  
  
  
  duplicates = as.vector(unlist(esetsAndDups$duplicates))
  
  sub = sub[which(sub$sample_name%!in%duplicates),]
  
  
  
  ######################################         Preparation of the final file      ########################################################################
  
  # PAM50 
  
  files = list.files("data_bitbucket/PAM50_ROR_P/PAM50_results//")
  


  
  i = 1
  lista = list()
  for(i in 1:length(files)){
    
    lista[[i]] = read.delim(paste("data_bitbucket/PAM50_ROR_P/PAM50_results//",files[i],"/MedCenForPAM50_Genes_pam50scores.txt",sep = ""))
  }
  
  pam50 = do.call(rbind.data.frame,lista)
  
  
  pam50 = pam50[which(pam50$X%in%sub$sample_name),]
  colnames(pam50)[1] = "samples"
  
  
  # GGI
  
  
  load("data_bitbucket/final_results/ggi_results")
  
  
  lista = list()
  i = 1
  for(i in 1:length(ggi_genefu)){
    
    col1 = names(ggi_genefu[[i]]$risk)
    col2 = as.numeric(as.vector(ggi_genefu[[i]]$risk))
    
    dat = as.data.frame(cbind(col1,col2))
    lista[[i]] = dat
  }
  
  ggi = do.call(rbind.data.frame,lista)
  
  ggi = ggi[which(ggi$col1%in%sub$sample_name),]
  colnames(ggi)[1] = "samples"
  colnames(ggi)[2] = "ggi_risk"
  
  # MAMMAPRINT
  
  load("data_bitbucket/final_results/mammaprint_results")
  
  
  lista = list()
  i = 1
  for(i in 1:length(mammaprint_adaptative_newmap)){
    
    
    
    
    lista[[i]] = mammaprint_adaptative_newmap[[i]]$result
  }
  
  mammaprint = do.call(rbind.data.frame,lista)
  mammaprint[,4] = rownames(mammaprint)
  mammaprint = mammaprint[,c(4,1,2,3)]
  colnames(mammaprint)[1] = "samples"
  mammaprint = mammaprint[which(mammaprint$samples%in%sub$sample_name),]
  
  
  
  
  # ONCOTYPE
  
load("data_bitbucket/final_results//oncotype_results")
  
  
  lista = list()
  i = 1
  for(i in 1:length(oncotype_results)){
    
    
    
    
    lista[[i]] = oncotype_results[[i]]$result
  }
  
  oncotype = do.call(rbind.data.frame,lista)
  oncotype[,8] = rownames(oncotype)
  oncotype = oncotype[,c(8,1:7)]
  colnames(oncotype)[1] = "samples"
  oncotype = oncotype[which(oncotype$samples%in%sub$sample_name),]
  colnames(oncotype)[8] = "oncotype_class"
  
  
  
  
  dat = left_join(pam50,mammaprint,by = "samples")
  
  dat = left_join(dat,ggi,by = "samples")
  
  dat = left_join(dat,oncotype,by = "samples")
  
  megatable = merge(sub,dat, by.x = "sample_name",by.y = "samples")
  
  megatable = as.data.frame(megatable)
  
  # We will combine all the survival endpoints into one. For that we will prioritize from highest group to smallest. The biggest is dthbcstat2, then dmfs and then days to tumor recurrence
  unified_survival_days = numeric()
  unified_survival_status = vector()
  endpoint = vector()
  i = 1
  for(i in 1:nrow(megatable)){
    
    if(!is.na(megatable$dthbcstat2[i])==TRUE){
      
      unified_survival_days[i] = megatable$T_dth[i]
      unified_survival_status[i] = megatable$dthbcstat2[i]
      endpoint[i] = "dthbc"
      next
    }
    if(!is.na(megatable$dmfs_days[i])==TRUE){
      
      unified_survival_days[i] = megatable$dmfs_days[i]
      unified_survival_status[i] = megatable$dmfs_status[i]
      endpoint[i] = "dmfs"
      next
    }
    if(!is.na(megatable$days_to_tumor_recurrence[i])==TRUE){
      
      unified_survival_days[i] = megatable$days_to_tumor_recurrence[i]
      unified_survival_status[i] = megatable$recurrence_status[i]
      endpoint[i] = "tumor_recurrence"
      next
      
    }
  }
  
  
  
  megatable$unified_survival_days = unified_survival_days
  megatable$unified_survival_status = unified_survival_status
  megatable$endpoint = endpoint
  
  
  megatable[megatable=="norecurrence"]<-FALSE
  megatable[megatable=="recurrence"]<-TRUE
  
  megatable$unified_survival_days = as.numeric(as.vector(megatable$unified_survival_days))
  
  megatable[megatable=="0"]<-FALSE
  megatable[megatable=="1"]<-TRUE
  
  
  megatable[megatable=="neg"]<-"negative"
  megatable[megatable=="pos"]<-"positive"
  
  
  # Merge cell cycle results
  
  
  load("data_bitbucket/final_results/cell_cycle_results")
  
  
  sub = merge(megatable,cell_cycle_results[,c(1,2)],by = "sample_name")
  
  
  write.table(sub,"data_bitbucket/final_results/all_results", sep="\t", row.names=F, col.names = T, quote=F)
  
  return(sub)
  
}

results_final = merge_function()



