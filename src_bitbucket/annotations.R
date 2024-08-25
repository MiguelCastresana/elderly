

'%!in%' <- function(x,y)!('%in%'(x,y))

#################################################################################################
# Load esets from either your own directory or load them from MetaGxBreast
# require(MetaGxBreast)
# 
# esetsAndDups <- loadBreastEsets(loadString =  "majority",removeDuplicates = TRUE)
# 
# extra = loadBreastEsets(loadString = c("METABRIC", "TCGA"),removeDuplicates = TRUE)
# 
#################################################################################################

# Load
# load("data_bitbucket/esets_1_37")
# load("data_bitbucket/esets_38_39")

# Load first 37 datasets
if(!exists("esetsAndDups")) {
  load("data_bitbucket/esets_1_37")
}
# Load datasets 38 and 39 (TCGA and METABRIC)
if(!exists("extra")) {
  load("data_bitbucket/esets_38_39")
}

# Join all the esets since we cannot have them in the same file due to metagxbreast loading approach
i = 1
alldatasets = list()



groups = c(names(esetsAndDups$esets),names(extra$esets))

i = 1
c= 1
for(i in 1:length(groups)){
  
  
  
  if(i<38){
    
    alldatasets = c(alldatasets,esetsAndDups$esets[[i]])
  }else{
    
    alldatasets = c(alldatasets,extra$esets[[c]])
    c = c + 1
  }
  
  
}




split_annotation_genefu = function(annot){
  
  # It can happen that some of the annotations are merge with ///. We split them:
  
  mm<-apply(annot,1,paste,collapse=" ")
  pos = grep("///", mm)
  
  if(length(pos)>0){
    
    
    add = list()
    kk = 1
    for(kk in 1:length(pos)){
      
      separar = list()
      l = numeric()
      j = 1
      for(j in 1:ncol(annot)){
        
        
        
        out = str_split(annot[pos[kk],j], "///")
        clean = lapply(out, function(x) x[x != ""])
        if(length(clean[[1]])<1){
          clean[[1]] = NA
        }
        
        l[j] = sapply(clean,function(x) length(x))
        separar[[j]] = clean
        
      }
      
      l_ordered = l[order(-l)]
      first_col = rep(as.vector(annot[pos[kk],1]),(max(l)*l_ordered[2])/l[1])
      second_col = rep(unlist(separar[[2]]),(max(l)*l_ordered[2])/l[2])
      third_col = rep(unlist(separar[[3]]),(max(l)*l_ordered[2])/l[3])
      fourth_col = rep(annot[pos[kk],4],(max(l)*l_ordered[2])/l[4])
      
      
      
      new = as.data.frame(cbind(first_col,second_col,third_col,fourth_col))
      names(new) = names(annot)
      
      add[[kk]] = new
      
      
    }
    
    
    add = do.call(rbind.data.frame, add)
    add = add[complete.cases(add),]
    
    add = add[!(is.na(add[,2]) | add[,2]==""), ]
    add = add[!(is.na(add[,3]) | add[,3]==""), ]
    
    annot = annot[-pos,]
    
    annot = rbind(annot,add)
    
    
  }
  colnames(annot)[1] = "probe"
  
  return(annot)
  
}



split_annotation_gpl = function(annot){
  
  # It can happen that some of the annotations are merge with ///. We split them:
  
  mm<-apply(annot,1,paste,collapse=" ")
  pos = grep("///", mm)
  
  if(length(pos)>0){
    
    
    add = list()
    kk = 1
    for(kk in 1:length(pos)){
      
      separar = list()
      l = numeric()
      j = 1
      for(j in 1:ncol(annot)){
        
        l[j] = length(str_split(annot[pos[kk],j], "///")[[1]])
        separar[[j]] = str_split(annot[pos[kk],j], "///")
      }
      
      
      first_col = rep(as.vector(annot[pos[kk],1]),max(l))
      second_col = unlist(separar[[2]])
      
      
      new = as.data.frame(cbind(first_col,second_col))
      names(new) = names(annot)
      
      add[[kk]] = new
      
      
    }
    
    
    add = do.call(rbind.data.frame, add)
    add = add[complete.cases(add),]
    
    add = add[!(is.na(add[,2]) | add[,2]==""), ]
    
    
    annot = annot[-pos,]
    
    annot = rbind(annot,add)
    
    
  }
  colnames(annot)[1] = "probe"
  
  return(annot)
  
}



split_annotation_gpl_comma = function(annot){
  
  # It can happen that some of the annotations are merge with ///. We split them:
  
  mm<-apply(annot,1,paste,collapse=" ")
  pos = grep(",", mm)
  
  if(length(pos)>0){
    
    
    add = list()
    kk = 1
    for(kk in 1:length(pos)){
      
      separar = list()
      l = numeric()
      j = 1
      for(j in 1:ncol(annot)){
        
        l[j] = length(str_split(annot[pos[kk],j], ",")[[1]])
        separar[[j]] = str_split(annot[pos[kk],j], ",")
      }
      
      
      first_col = rep(as.vector(annot[pos[kk],1]),max(l))
      second_col = unlist(separar[[2]])
      
      
      new = as.data.frame(cbind(first_col,second_col))
      names(new) = names(annot)
      
      add[[kk]] = new
      
      
    }
    
    
    add = do.call(rbind.data.frame, add)
    add = add[complete.cases(add),]
    
    add = add[!(is.na(add[,2]) | add[,2]==""), ]
    
    
    annot = annot[-pos,]
    
    annot = rbind(annot,add)
    
    
  }
  colnames(annot)[1] = "probe"
  
  return(annot)
  
}




# summary datasets

# final_summary = read.delim(paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/SUMMARY",sep = ""))

load("data_bitbucket/SUMMARY_datasets")

datasets_notplatform = final_summary$dataset[which(final_summary$biomart%!in%total_annotation$platform)]


positions = match(datasets_notplatform, groups)


i = 1
plat_new = numeric()
for(i in 1:length(datasets_notplatform)){
  
  
  if(positions[i]>=38){
    
    
    var = featureData(extra$esets[[positions[[i]]-37]])
  }else{
    
    var = featureData(esetsAndDups$esets[[positions[i]]])
  }
  
  annot = pData(var)
  names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
  
  # data = data[,which(colnames(data)%in%annot$probe)]
  
  
  annot = split_annotation_genefu(annot)
  
  annot$EntrezGene.ID = as.numeric(as.vector(annot$EntrezGene.ID))
  

  if("ORIGINALID"%in%colnames(total_annotation)){
    
    dat = total_annotation[which(total_annotation$ProbeSet%in%annot$probe | total_annotation$ORIGINALID%in%annot$probe),]
    
  }else{
    
    dat = total_annotation[which(total_annotation$probe%in%annot$probe),]
    
  }
  

  pos = which.max(as.numeric(table(dat$platform)))
  
  if(length(pos)>0){
    
    
    plat_new[i] = names(table(dat$platform))[pos]
    
  }else{
    
    plat_new[i] = NA
    
  }
  
  
  
  print(i)
}


plat_new = as.data.frame(cbind(datasets_notplatform,plat_new))














