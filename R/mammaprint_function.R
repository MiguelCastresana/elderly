
mammaprint_function = function(){
  
  defaultW <- getOption("warn") 
  
  options(warn = -1) 
  
  run_mammaprint = function(annot){
    
    
    uppANKI=rownames(data) %in% annot$probe
    a=data[uppANKI,]
    
    
    rmA=rowMeans(a,na.rm=TRUE)
    mcA=a-rmA
    
    
    con=merge(mcA,annot, by.x="row.names", by.y="probe", all.x=T)
    dim(con)
    
    meann = function(sig_gene){
      
      ok = mean(sig_gene,na.rm = TRUE)
      return(ok)
    }
    final=aggregate(con,list(con$ORIGINALID), meann)
    
    
    
    values=final[,-c(2,ncol(final))]
    
    
    cv=tapply(oriannot$Averagegood, oriannot$ORIGINALID,function(x) x[1])
    
    
    rownames(values)<-values[,"Group.1"]
    
    x=data.matrix(cv)
    
    
    correlations=x[rownames(values),]
    length(correlations)
    
    
    
    
    correl=data.matrix(correlations)
    
    z=data.matrix(values)
    p=cor(correl, y = z, use = "complete.obs", method = c("pearson"))
    k=cor(correl, y = z, use = "complete.obs", method = c("kendall"))
    s=cor(correl, y = z, use = "complete.obs", method = c("spearman"))
    
    pp=t(p)
    kk=t(k)
    ss=t(s)
    
    ptypes=cbind(pp,kk,ss)
    colnames(ptypes)=c("pearson","kendall", "spearman")
    ptypes=ptypes[-1,]
    
    uppcorr=data.matrix(ptypes[,"pearson"])
    
    
    y1=as.matrix(correlations_70gene[,2])
    
    scaledNKIcor=rescale(y1, q=0.05,na.rm = TRUE)
    scaledUppcor=rescale(uppcorr, q=0.05,na.rm = TRUE)
    
    
    # Also, can fit a linear model and look at the equation of the line
    cof1 = lm(scaledNKIcor ~ y1)$coeff[1]
    
    cof2 = lm(scaledNKIcor ~ y1)$coeff[2]
    
    correlations_70gene = (cof2 * 0.3) + cof1
    
    
    cof2_1 = lm(scaledUppcor ~ uppcorr)$coeff[1]
    
    cof2_2 = lm(scaledUppcor ~ uppcorr)$coeff[2]
    
    x = (correlations_70gene - cof2_1)/cof2_2
    
    
    
    Upp_corr = data.frame(uppcorr)
    Upp_corr[,"Prognosis"]=ifelse(Upp_corr[,1] >= x, "good", "poor")
    Upp_corr[,"Prognosis_binary"] = ifelse(Upp_corr[,1] >= x, 0, 1)
    mammaprint_UPP_calls = Upp_corr
    table(mammaprint_UPP_calls$Prognosis_binary)
    colnames(mammaprint_UPP_calls) = c("mammaprint.pearson", "mammaprint.prognosis",
                                       "mammaprint.prognosis_binary")
    
    
    objeto = list(result = mammaprint_UPP_calls,ratio_mapped = length(correlations)/70,total_mapped = length(correlations),
                  missed_genes = c(NCBI70$ORIGINALID[which(NCBI70$ORIGINALID%!in%names(correlations))]))
    
    return(objeto)
  }
  
  
  
  
  
  
  i = 1
  c = 1
  
  mammaprint_results = list()
  
  
  for(i in 1:length(groups)){
    
    
    
    
    if(i ==17){
      
      
      objeto = list(result = NA,ratio_mapped = 0,total_mapped = 0)
      
      mammaprint_results[[i]] = objeto
      next
      
    }
    # Load data
    exp = exprs(alldatasets[[i]])
    
    
    data = exp
    
    
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
      NCBI70Affy = total_annotation[which(total_annotation[,16]%in%plat),]
      
      
      
      safe = NCBI70Affy[,c(11,1)]
      
      names(safe) = c("probe","ORIGINALID")
      
      # Some datasets have the original IDs as probe names
      
      safe1 = as.data.frame(cbind(NCBI70[,1],NCBI70[,1]))
      names(safe1) = c("probe","ORIGINALID")
      
      
      
      safe = as.data.frame(rbind(safe,safe1))
      
      safe = safe[complete.cases(safe),]
      
      safe = safe[!duplicated(safe),]
      
      NCBI70Affy = NCBI70[,c(1,6)]
      names(NCBI70Affy) = c("ORIGINALID","entrez")
      # We have NICKS. Now we will also combine it with the extra information we may get from merging genefu annot
      
      
      
      var = featureData(alldatasets[[i]])
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      annot = split_annotation_genefu(annot)
      
      annot = annot[,c(1,3)]
      
      names(annot) = c("probe","entrez")
      
      annot$entrez = as.numeric(as.vector(annot$entrez))
      NCBI70Affy$entrez = as.numeric(as.vector(NCBI70Affy$entrez))
      
      annot = annot[which(!is.na(annot$entrez)),]
      NCBI70Affy = NCBI70Affy[which(!is.na(NCBI70Affy$entrez)),]
      
      
      
      
      dat = left_join(NCBI70Affy,annot, by = "entrez")
      
      
      res = dat[,c(3,1)]
      
      res = rbind(safe, res)
      
      res = res[complete.cases(res),]
      
      res = res[!duplicated(res),]
      
      
      annot = res
      
      objeto = run_mammaprint(annot)
      
      mammaprint_results[[i]] = objeto
    }else{
      
      
      # We coul not find any platform information for these datasets so we need to rely on the annotation from Genefu
      
      var = featureData(alldatasets[[i]])
      
      
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      # It can happen that some of the annotations are merge with ///. We split them:
      
      annot = split_annotation_genefu(annot)
      
      
      annot$EntrezGene.ID = as.numeric(as.vector(annot$EntrezGene.ID))
      annot = annot[which(!is.na(annot$EntrezGene.ID)),]
      
      # We map this annotation to the original mammaprint annotation file so we can link our genes/probes to Original IDs
      merge1 = merge(annot, NCBI70, by.x="EntrezGene.ID", by.y="LLID", all.x=T)
      
      merge1 = merge1[which(!is.na(merge1$probe) & !is.na(merge1$ORIGINALID)),]
      
      annot = merge1[,c(2,5)]
      
      if(groups[i]%in%names(gpl_datasets)){
        
        pos = which(names(gpl_datasets)%in%groups[i])
        
        annot_extra = gpl_datasets[[pos]]
        
        names(annot_extra) = c("probe","ORIGINALID")
        
        annot = rbind(annot,annot_extra)
        
        annot = annot[complete.cases(annot),]
        annot = annot[!duplicated(annot),]
      }
      
      
      # DATASET 17 has no good mappings
      if(nrow(annot)<1){
        
        
        objeto = list(result = NA,ratio_mapped = 0,total_mapped = 0)
        
        mammaprint_adaptative_newmap[[i]] = objeto
      }else{
        
        
        objeto = run_mammaprint(annot)
        mammaprint_results[[i]] = objeto
        
        
      }
      
      
    }
    
    
    
    
    print(i)
  }
  
  
  
  
  return(mammaprint_results)
  options(warn = defaultW)
  

}

mammaprint_results = mammaprint_function()



save(mammaprint_results,file = "data_bitbucket/final_results/mammaprint_results")






