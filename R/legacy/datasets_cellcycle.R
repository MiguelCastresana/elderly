



'%!in%' <- function(x,y)!('%in%'(x,y))


# Package loading
require(dplyr)
require(stringi)
require(stringr)
require(genefu)
require(MetaGxBreast)
require(ggplot2)





### LOAD THE DATA AND COMBINE IT

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

db_annotation = list()

files = list.files("C:/Users/migcas//OneDrive - Karolinska Institutet/Desktop/Postdoc/data/annotation_package_ALL_genes//")

dbs = sub(".RData.*", "", files) 
vec = vector()
i = 1
for(i in 1:length(files)){
  
  
  d = loadRData(paste("C:/Users/migcas//OneDrive - Karolinska Institutet/Desktop/Postdoc/data/annotation_package_ALL_genes//",files[i],sep=""))
  db_annotation[[i]] = d
  
  vec = c(vec,rep(dbs[i],nrow(d)))
}


sapply(db_annotation, function(x) ncol(x))


db_annotation = do.call(rbind.data.frame, db_annotation)



db_annotation[,9] = vec






colnames(db_annotation)[1] = "probe"
colnames(db_annotation)[2] = "entrez"
colnames(db_annotation)[3] = "symbol"
colnames(db_annotation)[9] = "platform"

db_annotation = db_annotation[complete.cases(db_annotation[,2]),]
db_annotation = db_annotation[!duplicated(db_annotation),]




total_annotation = db_annotation


write.table(total_annotation,"C:/Users/migcas//OneDrive - Karolinska Institutet/Desktop/Postdoc/data/total_annotation_all_genes", sep = "\t", quote = FALSE, row.names = FALSE)







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




########## GPL information       ###################################




# Load GPLs and parse the data

files = list.files(paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/data/GPL/",sep=""))

gpl_datasets = list()

i = 1
for(i in 1:length(files)){

  annot = read.delim(paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/data/GPL/",files[i],sep=""))


  g = sub("\\.txt.*", "", files[i])
  pos = which(groups%in%g)

  exp = exprs(alldatasets[[pos]])

  var = featureData(alldatasets[[pos]])
  annot_genefu = pData(var)
  names(annot_genefu) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")


  if(colnames(annot)[2]=="Gene.title"){

    annot = annot[,c(1,3)]
    annot = annot[complete.cases(annot),]
    annot = annot[!(is.na(annot[,2]) | annot[,2]==""), ]

    annot = split_annotation_gpl(annot)

   

    annot = annot[complete.cases(annot),]
    annot = annot[!duplicated(annot),]
    gpl_datasets[[i]] = annot


  }else if(colnames(annot)[2] == "Comment.OLIGO_ID."){


    annot = annot[,c(2,3)]
    annot = annot[complete.cases(annot),]
    annot = annot[!(is.na(annot[,2]) | annot[,2]==""), ]

    annot = split_annotation_gpl(annot)



    names(annot) = c("probe","symbol")

   
    annot = annot[complete.cases(annot),]
    annot = annot[!duplicated(annot),]
    gpl_datasets[[i]] = annot

  }else if(colnames(annot)[2]=="Search_key"){

    annot[,1] = paste("probe_", annot$ID, sep="")

    # take only the cases where ID is NM*

    annot = annot[annot$ID%in%rownames(exp),]

    annot = annot[,c(1,7)]

    names(annot) = c("probe","symbol")

    
    gpl_datasets[[i]] = annot
  }else if(colnames(annot)[2]=="PLATE"){

    annot = annot[which(annot[,1]%in%rownames(exp)),]

    annot = annot[,c(1,7)]

    annot = split_annotation_gpl_comma(annot)

   
    gpl_datasets[[i]] = annot

  }else{

    gpl_datasets[[i]] = NA

  }

  print(i)

}

names(gpl_datasets) = sub("\\.txt.*", "", files)






save(gpl_datasets, file = paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/data/gpl_annotation_missing_datasets_all_genes",sep = ""))









# Load complete mapping

total_annotation = read.delim(paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/data/total_annotation_all_genes",sep = ""))

# total_annotation$gene[which(total_annotation$gene == "PALM2AKAP2")] <- "PALM2-AKAP2"


# summary datasets

final_summary = read.delim(paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/SUMMARY",sep=""))




# Load missing dataset annotation

load(paste("C:/Users/",place,"/OneDrive - Karolinska Institutet/Desktop/Postdoc/data/gpl_annotation_missing_datasets_all_genes",sep = ""))











# Giving platform annotation based on most mapped platform to those annotation without platform



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
  
  dat = total_annotation[which(total_annotation$probe%in%annot$probe ),]
  
  pos = which.max(as.numeric(table(dat$platform)))
  
  if(length(pos)>0){
    
    
    plat_new[i] = names(table(dat$platform))[pos]
    
  }else{
    
    plat_new[i] = NA
    
  }
  
  # if(i >3){break}
  
  
  print(i)
}


plat_new = as.data.frame(cbind(datasets_notplatform,plat_new))









# Home
place = "Miguel"

# Work
place = "migcas"




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
    
    # Extract the annotation file from Nick´s approach (Total annotation is the file that gathers all the platforms that I could find)
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


save(data_clean,file = "C:/Users/migcas/OneDrive - Karolinska Institutet/Desktop/Postdoc/merged_data/data_clean_alldatasets_merged")


load("C:/Users/Miguel//OneDrive - Karolinska Institutet/Desktop/Postdoc/merged_data/data_clean_alldatasets_merged")

# Which datasets dont have platform available? Since it is giving me issues running PCA due to NAs






i = 1
c = 1
cases = numeric()
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
  
  if(length(plat)<1){
    cases = c(cases,i)
  }
  
  print(i)
}




# Taking all of them

grupo = vector()      
GEX_data = data.frame()
i = 1
for(i in 1:length(data_clean)){
  

  if(i == 17){next}
  grupo = c(grupo,rep(groups[i],ncol(data_clean[[i]][[1]])))
  GEX_data <- merge(GEX_data, data_clean[[i]][[1]],
                    by = 'row.names', all = TRUE)
  
  rownames(GEX_data) = GEX_data[,1]
  GEX_data = GEX_data[,-1]
  
  print(i)
  
}

GEX_data = GEX_data[-1,]



# Filtering some datasets


'%!in%' <- function(x,y)!('%in%'(x,y))


# We will try to pick only these platforms
# selection = final_summary$dataset[which(!is.na(final_summary$biomart))]
selection = final_summary$dataset[which(final_summary$biomart%in%c("affy_hg_u133a","affy_hg_u95av2","affy_hg_u133_plus_2","affy_hg_u133a///affy_hg_u133_plus_2",
                                                                   "affy_hg_u133b///affy_hg_u133a","affy_hugene_1_0_st_v1","affy_u133_x3p"))]

pos = which(groups%in%selection)

grupo = vector()      
GEX_data = data.frame()
i = 1
for(i in 1:length(data_clean)){
  
  if(i%!in%pos){next}
  if(i == 17){next}
  grupo = c(grupo,rep(groups[i],ncol(data_clean[[i]][[1]])))
  GEX_data <- merge(GEX_data, data_clean[[i]][[1]],
                    by = 'row.names', all = TRUE)
  
  rownames(GEX_data) = GEX_data[,1]
  GEX_data = GEX_data[,-1]
  
  print(i)
  
}

GEX_data = GEX_data[-1,]


save(GEX_data,file = "C:/Users/Miguel//OneDrive - Karolinska Institutet/Desktop/Postdoc/merged_data/GEX_data_alldatasets_merged_datasets_all_affymetrix")

load("C:/Users/Miguel//OneDrive - Karolinska Institutet/Desktop/Postdoc/merged_data/GEX_data_alldatasets_merged")




grupo = vector()      
i = 1
for(i in 1:length(data_clean)){
  
  if(i%!in%pos){next}
  if(i == 17){next}
  grupo = c(grupo,rep(groups[i],ncol(data_clean[[i]][[1]])))
 
  
  print(i)
  
}



a = as.data.frame(rownames(GEX_data))



# valor = colSums(is.na(data))



# 
# 
# valor = rowSums(is.na(data))
# 
# valor = which(valor<8000)
# 
# select = names(valor)
# 
# 
# data = data[,which(colnames(data)%!in%select)]



library(tidyr)
GEX_data = GEX_data %>% drop_na()


# Perform the PCA



data = GEX_data
data = t(data)
sample_info = as.data.frame(cbind(rownames(data),grupo))




# data = as_tibble(data)
rownames(data) = colnames(GEX_data)


sample_pca <- prcomp(na.omit(data), scale=TRUE)





pc_eigenvalues <- sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues


pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")


pc_scores <- sample_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores



pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()



# 
a = pc_scores$sample
sample_info = sample_info[which(sample_info$V1%in%a),]

names(sample_info) = c("sample","study")

pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(study))) +
  geom_point(size=3) +  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlim(-20,20) +ylim(-50,50)

# print the result (in this case a ggplot)
pca_plot









