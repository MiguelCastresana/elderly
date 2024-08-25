
###########  ANALYSIS    #######################################



oncotype_function = function(x){
  
  # Make a recurrence score function
  oneToFifteen <- function(x) {
    x <- x - min(x)
    x <- x * (15/max(x))
    x
  }
  GRB7groupscore <- function(arr) {
    arr["GRB7", ] <- arr["GRB7", ] * 0.9
    arr["ERBB2", ] <- arr["ERBB2", ] * 0.1
    gs <- apply(arr[c("GRB7", "ERBB2"), ], 2, sum)
    gs[gs < 8] <- 8
    gs
  }
  ERgroupscore <- function(arr) {
    
    if(length(arr["ESR1", ])>0){arr["ESR1", ] <- arr["ESR1", ] * 0.2}
    if(length(arr["PGR", ])>0){arr["PGR", ] <- arr["PGR", ] * 0.3}
    if(length(arr["BCL2", ])>0){arr["BCL2", ] <- arr["BCL2", ] * 0.25}
    if(length(arr["SCUBE2", ])>0){arr["SCUBE2", ] <- arr["SCUBE2", ] * 0.25}
    
    
    gs <- apply(arr[c("ESR1", "PGR", "BCL2", "SCUBE2"), ], 2, sum)
    gs
  }
  proliferationgroupscore <- function(arr) {
    arr["AURKA", ] <- arr["AURKA", ] * 0.2
    arr["MKI67", ] <- arr["MKI67", ] * 0.2
    arr["MYBL2", ] <- arr["MYBL2", ] * 0.2
    arr["CCNB1", ] <- arr["CCNB1", ] * 0.2
    arr["BIRC5", ] <- arr["BIRC5", ] * 0.2
    gs <- apply(arr[c("AURKA", "MKI67", "MYBL2", "CCNB1", "BIRC5"), ], 2,
                sum)
    gs[gs < 6.5] <- 6.5
    gs
  }
  invasiongroupscore <- function(arr) {
    arr["CTSL2", ] <- arr["CTSL2", ] * 0.5
    arr["MMP11", ] <- arr["MMP11", ] * 0.5
    gs <- apply(arr[c("CTSL2", "MMP11"), ], 2, sum)
    gs
  }
  RS <- function(arr) {
    arr["GRB7gs", ] <- arr["GRB7gs", ] * 0.47
    arr["ERgs", ] <- arr["ERgs", ] * -0.34
    arr["proliferationgs", ] <- arr["proliferationgs", ] * 1.04
    arr["invasiongs", ] <- arr["invasiongs", ] * 0.1
    arr["CD68", ] <- arr["CD68", ] * 0.05
    arr["GSTM1", ] <- arr["GSTM1", ] * -0.08
    arr["BAG1", ] <- arr["BAG1", ] * -0.07
    RSu <- apply(arr, 2, sum)
    RS <- 20 * (RSu - 6.7)
    RS[RS < 0] <- 0
    RS[RS > 100] <- 100
    class <- rep("intermediate", length(RS))
    class[RS > 31] <- "high"
    class[RS < 18] <- "low"
    rbind(RSu, RS, class)
  }
  paikfunction <- function(mat, group) {
    
    
    
    meann = function(sig_gene){
      
      ok = mean(sig_gene,na.rm = TRUE)
      return(ok)
    }
    
    
    
    reference <- apply(mat[group, ], 2, meann)
    matr <- t(apply(mat, 1, function(x) x - reference))
    # print(matr[,1:5])
    matr <- apply(matr, 1, function(x) oneToFifteen(x))
    matr <- t(matr)
    # print(matr[,1:25])
    GRB7gs <- GRB7groupscore(matr)
    ERgs <- ERgroupscore(matr)
    proliferationgs <- proliferationgroupscore(matr)
    invasiongs <- invasiongroupscore(matr)
    mat2 <- rbind(GRB7gs, ERgs, proliferationgs, invasiongs, CD68 = matr["CD68",
    ], GSTM1 = matr["GSTM1", ], BAG1 = matr["BAG1", ])
    
    o <- RS(mat2)
    rbind(mat2[c("GRB7gs", "ERgs", "proliferationgs", "invasiongs"), ], o)
  }
  
  
  onco_data_input = function(annot){
    
    data = data[which(rownames(data)%in%annot$probe),]
    
    
    annot = annot[match(rownames(data), annot$probe),]
    
    
    paikRMA <- data[as.vector(annot$probe), ]
    
    rownames(paikRMA) = annot$symbol
    
    
    
    # We select the symbols we have information about so then we can do the mean expression of those
    select = as.vector(annot$symbol[which(annot$symbol%in%rownames(paikRMA))])
    
    
    meann = function(sig_gene){
      
      ok = mean(sig_gene,na.rm = TRUE)
      return(ok)
    }
    
    paikRMA <- apply(paikRMA, 2, function(x) tapply(x, select, meann))
    
    
    
    # Where are the reference genes
    gr <- annot[match(rownames(paikRMA), annot$symbol), "group"] == "reference"
    
    # Check if we are missing a gene that is not a reference gene because then we cannot perform the analysis
    ll = sig.oncotypedx[which(sig.oncotypedx$symbol%!in%rownames(paikRMA)),3]
    missing = sig.oncotypedx[which(sig.oncotypedx$symbol%!in%rownames(paikRMA)),c(1,4)]
    ref_missing = missing[which(missing$group=="reference"),]
    ll = table(missing$group)
    not_miss = c("undefined","proliferation","invasion","her2","estrogen")
    not = not_miss[which(not_miss%in%names(ll))]
    
    
    
    if(length(not)>0){
      
      
      objeto = list(result = NA,ratio_mapped = (nrow(paikRMA)-(5-nrow(ref_missing))) /16, total_mapped_without_ref = (nrow(paikRMA)-(5-nrow(ref_missing))),
                    missing = as.vector(sig.oncotypedx[which(sig.oncotypedx$symbol%!in%rownames(paikRMA)),1]),class = not,reference_genes_missing = ref_missing[,1])
      
      
      
    }else{
      
      objeto = list(result = as.data.frame(t(data.frame(paikfunction(paikRMA, gr)))),ratio_mapped = (nrow(paikRMA)-(5-nrow(ref_missing))) /16, total_mapped_without_ref = (nrow(paikRMA)-(5-nrow(ref_missing))),
                    missing = as.vector(sig.oncotypedx[which(sig.oncotypedx$symbol%!in%rownames(paikRMA)),1]), class = not,reference_genes_missing = ref_missing[,1])
      
    }
    
    return(objeto)
  }
  
  
  
  
  
  
  
  i = 1
  c = 1
  
  oncotype_results = list()
  problematic_datasets = vector()
  
  
  for(i in 1:length(groups)){
    
    
    # Load data
    exp = exprs(alldatasets[[i]])
    
    
    data = exp
    
    
    
    # Filter healthy people
    
    pheno = pData(alldatasets[[i]])
    
    pos = which(pheno$sample_type%in%"tumor")
    
    pheno = pheno[pos,]
    data = data[,pos]
    
    
    # Select platform for the dataset
    
    plat = as.vector(final_summary$biomart[which(final_summary$dataset%in%groups[i])])
    
    plat = plat[which(plat%in%total_annotation$platform)]
    
    # Extra platform information based on maximum number of probes mapped to the annotation provided by Genefu
    plat_added = as.vector(plat_new[which(plat_new$datasets_notplatform%in%groups[i]),2])
    
    plat = c(plat,plat_added)
    plat = plat[!is.na(plat)]
    
    
    if(length(plat)>=1){
      
      # Extract the annotation file from Nick?s approach (Total annotation is the file that gathers all the platforms that I could find)
      NCBI70Affy = total_annotation[which(total_annotation[,5]%in%plat),]
      
      
      safe = NCBI70Affy[,c(1:4)]
      
      names(safe) = c("probe","entrez","symbol","group")
      
      
      NCBI70Affy = NCBI70[,c(3,1,2,4)]
      names(NCBI70Affy) = c("entrez","symbol","probe","group")
      
      # We have NICKS. Now we will also combine it with the extra information we may get from merging genefu annot
      
      
      
      var = featureData(alldatasets[[i]])
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      annot = split_annotation_genefu(annot)
      
      
      annot = annot[,c(1,3,2)]
      
      names(annot) = c("probe","entrez","symbol")
      
      annot$entrez = as.numeric(as.vector(annot$entrez))
      
      annot = annot[which(!is.na(annot$entrez)),]
      
      
      dat = merge(NCBI70Affy,annot, by = "entrez",all.x= TRUE)
      
      
      res = dat[,c(5,1,2,4)]
      
      names(res) = c("probe","entrez","symbol","group")
      
      res1 = dat[,c(3,1,2,4)]
      
      names(res1) = names(res)
      
      res = rbind(safe, res,res1)
      
      res = res[complete.cases(res),]
      
      res = res[!duplicated(res),]
      
      
      annot = res
      
      annot = annot[which(annot$symbol%in%NCBI70$symbol),]
      
      
      oncotype_results[[i]] = onco_data_input(annot)
      
      
    }else{
      
      
      # We could not find any platform information for these datasets so we need to rely on the annotation from Genefu
      
      var = featureData(alldatasets[[i]])
      
      
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      # It can happen that some of the annotations are merge with ///. We split them:
      
      annot = split_annotation_genefu(annot)
      
      
      annot$EntrezGene.ID = as.numeric(as.vector(annot$EntrezGene.ID))
      annot = annot[which(!is.na(annot$EntrezGene.ID)),]
      
      # We map this annotation to the original mammaprint annotation file so we can link our genes/probes to Original IDs
      merge1 = merge(annot, NCBI70, by.x="EntrezGene.ID", by.y="EntrezGene.ID", all.x=T)
      
      merge1 = merge1[which(!is.na(merge1$probe) & !is.na(merge1$NCBI.gene.symbol)),]
      
      annot = merge1[,c(2,1,3,7)]
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
      
      names(annot) = c("probe","entrez","symbol","group")
      
      annot = annot[which(annot$symbol%in%NCBI70$symbol),]
      
      if(groups[i]%in%names(gpl_datasets)){
        
        pos = which(names(gpl_datasets)%in%groups[i])
        
        annot_extra = gpl_datasets[[pos]]
        
        names(annot_extra) = c("probe","entrez","symbol","group")
        
        annot = rbind(annot,annot_extra)
        
        annot = annot[complete.cases(annot),]
        annot = annot[!duplicated(annot),]
        
        # With this we get rid of the other symbols for the same gene which are not used
        annot = annot[which(annot$symbol%in%NCBI70$symbol),]
      }
      
      
      if(nrow(annot)<1){
        
        
        oncotype_results[[i]] = onco_data_input(annot)
        
      }else{
        
        
        annot = annot[which(annot$symbol%!in%"CTSV"),]
        
        oncotype_results[[i]] = onco_data_input(annot)
        
        
      }
      
      
    }
    
    
    
    
    print(i)
  }
  
  return(oncotype_results)
}

oncotype_results = oncotype_function()


save(oncotype_results,file = "data_bitbucket/final_results/oncotype_results")

