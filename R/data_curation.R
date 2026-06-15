




# We curate all datasets accordingly for further use

i = 1
c = 1

data_clean = list()


for(i in 1:length(groups)){
  
  
  if(i ==17){data_clean[[i]] = NA
  next}
  
  # Load data
  exp = exprs(alldatasets[[i]])
  
  
  data = exp
  
  
  # data = data[colSums(is.na(data)) != ncol(data), ]
  
  
  
  # Filter healthy people
  
  pheno = pData(alldatasets[[i]])
  
  pos = which(pheno$sample_type%in%"tumor")
  
  pheno = pheno[pos,]
  data = data[,pos]
  
  
  # Mean normalization of the probes
  rm=rowMeans(data)
  
  data=data-rm
  
  # Select platform for the dataset
  
  plat = as.vector(final_summary$biomart[which(final_summary$dataset%in%groups[i])])
  
  plat = plat[which(plat%in%total_annotation$platform)]
  
  # Extra platform information based on maximum number of probes mapped to the annotation provided by Genefu
  plat_added = as.vector(plat_new[which(plat_new$datasets_notplatform%in%groups[i]),2])
  
  plat = c(plat,plat_added)
  plat = plat[!is.na(plat)]
  
  
  if(length(plat)>=1){
    
    # Extract the annotation file from Nick?s approach (Total annotation is the file that gathers all the platforms that I could find)
    NCBI70Affy = total_annotation[which(total_annotation[,9]%in%plat),]
    
    
    
    safe = NCBI70Affy[,c(1:3)]
    
    
    safe = safe[complete.cases(safe),]
    
    safe = safe[!duplicated(safe),]
    
    
    
    var = featureData(alldatasets[[i]])
    annot = pData(var)
    names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
    
    
    annot = split_annotation_genefu(annot)
    
    annot = annot[,c(1,2,3)]
    
    names(annot) = c("probe","symbol","entrez")
    
    
    
    
    
    # dat = left_join(NCBI70Affy,annot, by = "symbol")
    # 
    # 
    # res = dat[,c(3,1)]
    # 
    # names(res) = c("probe","symbol")
    
    safe = safe[which(safe$probe%in%rownames(data)),]
    
    res = rbind(safe, annot)
    
    
    annot = res
    
    annot = annot[which(annot$probe%in%rownames(data)),]
    
    annot$probe = trimws(annot$probe,"both")
    annot$symbol = trimws(annot$symbol,"both")
    
    # annot = annot[which(as.vector(annot$symbol)%in%NCBI70$probe),]
    
    
    
    # Collapse probes for single gene-probe
    annot = annot[complete.cases(annot),]
    
    annot = annot[!duplicated(annot),]
    
    data = data[which(rownames(data)%in%annot$probe),]
    meann = function(sig_gene){
      
      ok = mean(sig_gene,na.rm = TRUE)
      return(ok)
    }
    
    con=merge(data,annot, by.x="row.names", by.y="probe", all.x=T)
    
    
    
    final=aggregate(con,list(con$symbol), meann)
    
    values=final[,-c(2,ncol(final),(ncol(final)-1))]
    
    rownames(values)<-values[,"Group.1"]
    
    
    
    # objeto = run_mammaprint(annot)
    # 
    # mammaprint_adaptative_newmap[[i]] = objeto
  }else{
    
    
    # We coul not find any platform information for these datasets so we need to rely on the annotation from Genefu
    
    var = featureData(alldatasets[[i]])
    
    
    annot = pData(var)
    names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
    
    annot = annot[!duplicated(annot[,3]),]
    
    
    # It can happen that some of the annotations are merge with ///. We split them:
    
    annot = split_annotation_genefu(annot)
    
    
    annot$EntrezGene.ID = as.numeric(as.vector(annot$EntrezGene.ID))
    annot = annot[which(!is.na(annot$EntrezGene.ID)),c(1:2)]
    
    
    
    
    names(annot) = c("probe","symbol")
    
    
    if(groups[i]%in%names(gpl_datasets)){
      
      pos = which(names(gpl_datasets)%in%groups[i])
      
      annot_extra = gpl_datasets[[pos]]
      
      names(annot_extra) = c("probe","symbol")
      
      annot = rbind(annot,annot_extra)
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
    }
    
    annot = annot[which(annot$probe%in%rownames(exp)),]
    
    
    
    
    data = data[which(rownames(data)%in%annot$probe),]
    meann = function(sig_gene){
      
      ok = mean(sig_gene,na.rm = TRUE)
      return(ok)
    }
    
    
    con=merge(data,annot, by.x="row.names", by.y="probe", all.x=T)
    
    final=aggregate(con,list(con$symbol), meann)
    
    values=final[,-c(2,ncol(final))]
    
    rownames(values)<-values[,"Group.1"]
    
    
  }
  
  pheno = pheno[,c("sample_name","er")]
  
  values = values[,-1]
  
  ok = list(values,annot,pheno)
  
  
  data_clean[[i]] = ok
  
  
  print(i)
}
