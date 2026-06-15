


ggi_function = function(){
  
  
  # Information required to run the rankclass method for GGI (compute grades if we dont have that information)
  
  
  # the following defines a vector of 91 HGU-133A probe sets used for simplified grading
  
  u133a<-c("202954_at","222077_s_at","201088_at","203554_x_at","218355_at","210052_s_at","202580_x_at","204092_s_at","218755_at","201584_s_at",
           "204825_at","202095_s_at","201710_at","211762_s_at","209680_s_at","209408_at","202870_s_at","38158_at","202107_s_at","204767_s_at",
           "203046_s_at","210559_s_at","221520_s_at","214710_s_at","210821_x_at","218726_at","201475_x_at","204033_at","202705_at","204649_at",  
           "209836_x_at","203276_at","203755_at","203214_x_at","212949_at","219510_at","204026_s_at","203432_at","204768_s_at","209773_s_at",
           "214431_at","218883_s_at","205733_at","222039_at","214096_s_at","211072_x_at","202779_s_at","218447_at","213911_s_at","212141_at",
           "221591_s_at","209251_x_at","217835_x_at","201890_at","213671_s_at","218009_s_at","207828_s_at","204695_at","212021_s_at","201090_x_at",
           "218039_at","204603_at","222036_s_at","204252_at","201524_x_at","204817_at","221436_s_at","201195_s_at","208696_at","218556_at",  
           "211058_x_at","212723_at","203022_at","210334_x_at","203744_at","213103_at","204703_at","218346_s_at","218471_s_at","205898_at",  
           "204072_s_at","219455_at","217889_s_at","219238_at","216264_s_at","221562_s_at","216520_s_at","220141_at","218483_s_at","221771_s_at",
           "220917_s_at")
  
  # for other than HGU-133A platforms, a corresponding vector of Entrez Gene identifiers is utilized
  
  entrezgene<-c("11065","29127","3838","9232","24137","22974","2305","6790","10112","10212","9833","332","4605","3838","3833","11004", 
                "991","9700","4171","2237","8914","983","55143","891","1058","55355","4141","9319","9133","10024","552900","4001",
                "701","983","23397","10721","11130","7112","2237","6241","8833","79682","641","146909","56901","10376","27338","56942", 
                "3015","4173","54478","84790","55969","6241","4141","9055","1063","993","4288","10376","51203","9156","4173","1017",  
                "7334","9700","83461","8140","22948","29095","10376","23210","10535","332","3149","90627","8100","27244","582","1524",  
                "10129","79846","79901","55650","3913","23410","7178","79864","56912","54737","57728")
  
  # a vector of 1 and -1 describes which genes are upregulated in grade 3 and grade 1 tumours, respectively
  
  weights<-c(rep(1,75),rep(-1,16))
  
  
  
  
  rankclass = function(arraydata,annotation){
    
    
    arraydata = prueba
    
    rows <- rownames(arraydata)
    annotation="entrezgene"
    
    
    
    sub = annot[which(annot$probe%in%rows),c(1,2)]
    
    
    
    
    rownames(arraydata) = sub[,2]
    rows <- rownames(arraydata)
    
    if(annotation == "u133a"){signature<-u133a}
    if(annotation == "entrezgene"){signature<-entrezgene}
    
    use<-signature %in% rows
    signature<-signature[use]
    weights<-weights[use]
    
    upgenes1<-signature[weights==-1]
    upgenes3<-signature[weights==1]
    
    rc<-function(expvec,rows,upgenes1,upgenes3){
      
      rank_expvec<-rank(expvec)
      GGLrank<-mean(rank_expvec[rows %in% upgenes1])
      GGHrank<-mean(rank_expvec[rows %in% upgenes3])
      GGR<-1
      if(GGHrank > GGLrank){GGR<-3}
      GGR
    }
    
    if(is.vector(arraydata)){G<-rc(arraydata,rows,upgenes1,upgenes3)}
    if(is.array(arraydata)){G<-apply(arraydata,2,function(x) rc(x,rows,upgenes1,upgenes3))}
    G
  }
  
  
  # I modified the genefu GGI function to allow us to collapse the probes by genes, in case we do not want to use their mapping system
  
  
  ggi_modified = function (data, annot, do.mapping = FALSE, mapping, hg, verbose = FALSE) 
  {
    scale.raw.ggi <- function(ggi, hg) {
      if (length(hg) != length(ggi)) {
        stop("bad length of hg!")
      }
      mhg1 <- mean(ggi[hg == 1], na.rm = TRUE)
      mhg3 <- mean(ggi[hg == 3], na.rm = TRUE)
      mm <- mhg1 + (mhg3 - mhg1)/2
      res.scaled <- ((ggi - mm)/(mhg3 - mhg1)) * 2
      return(res.scaled)
    }
    ggi.gl <- cbind(sig.ggi[, c("probe", "EntrezGene.ID")], 
                    coefficient = ifelse(sig.ggi[, "grade"] == 1, -1, 
                                         1))
    tt <- sig.score_modified(x = ggi.gl, data = data, annot = annot, do.mapping = do.mapping, 
                             mapping = mapping, signed = TRUE, verbose = verbose)
    
    myprobe <- tt$probe
    mymapping <- tt$mapping
    res <- tt$score
    if (!missing(hg)) {
      if (length(hg) != nrow(data)) {
        stop("hg must have the same length nrow(data)!")
      }
      mhg1 <- mean(res[hg == 1], na.rm = TRUE)
      mhg3 <- mean(res[hg == 3], na.rm = TRUE)
      mm <- mhg1 + (mhg3 - mhg1)/2
      res.scaled <- ((res - mm)/(mhg3 - mhg1)) * 2
      res <- list(score = res.scaled, risk = ifelse(res.scaled >= 
                                                      0, 1, 0), mapping = mymapping, probe = myprobe)
    }
    else {
      riskt <- rep(NA, length(res))
      names(riskt) <- names(res)
      res <- list(score = res, risk = riskt, mapping = mymapping, 
                  probe = myprobe)
    }
    return(res)
  }
  
  
  size = 0
  
  sig.score_modified = function(x = ggi.gl, data = data, annot = annot, do.mapping = do.mapping,cutoff = NA, 
                                mapping = mapping, signed = TRUE, verbose = verbose){
    
    if (missing(data) || missing(annot)) {
      stop("data and annot parameters must be specified")
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (nrow(x) == 0) {
      stop("empty gene list!")
    }
    myprobe <- as.character(x[, "probe"])
    mygid <- as.character(x[, "EntrezGene.ID"])
    mycoef <- as.numeric(x[, "coefficient"])
    names(mycoef) <- names(mygid) <- names(myprobe) <- myprobe
    nix <- order(abs(mycoef), decreasing = TRUE, na.last = NA)
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]
    if (do.mapping) {
      gid1 <- mygid
      gid2 <- as.character(annot[, "EntrezGene.ID"])
      names(gid2) <- dimnames(annot)[[1]]
      rm.ix <- is.na(gid1) | duplicated(gid1)
      gid1 <- gid1[!rm.ix]
      rr <- geneid.map(geneid1 = gid2, data1 = data, geneid2 = gid1, 
                       verbose = FALSE)
      if (is.na(rr$geneid1[1])) {
        res <- rep(NA, nrow(data))
        names(res) <- dimnames(data)[[1]]
        gf <- c(mapped = 0, total = nrow(x))
        if (verbose) {
          message(sprintf("probe candidates: 0/%i", 
                          nrow(x)))
        }
        return(list(score = res, mapping = gf, probe = cbind(probe = NA, 
                                                             EntrezGene.ID = NA, new.probe = NA)))
      }
      nix <- match(rr$geneid2, mygid)
      myprobe <- myprobe[nix]
      mygid <- mygid[nix]
      mycoef <- mycoef[nix]
      gid1 <- rr$geneid2
      if (is.null(names(gid1))) {
        stop("problem with annotations!")
      }
      gid2 <- rr$geneid1
      if (is.null(names(gid2))) {
        stop("problem with annotations!")
      }
      data <- rr$data1
      names(mycoef) <- names(mygid) <- mygid <- names(myprobe) <- myprobe <- as.character(gid1)
      dimnames(data)[[2]] <- as.character(gid2)
    }
    else {
      nix <- is.element(myprobe, dimnames(data)[[2]])
      myprobe <- myprobe[nix]
      mygid <- mygid[nix]
      mycoef <- mycoef[nix]
      gid1 <- gid2 <- mygid
      
      
      ##### ADDED CODE TO COLLAPSE PROBES FOR THE SAME GENE #########################
      meann = function(sig_gene){
        
        ok = mean(sig_gene,na.rm = TRUE)
        return(ok)
      }
      
      
      final=aggregate(t(data),list(colnames(data)), meann)
      # dim(final)
      
      
      final = t(final)
      
      colnames(final) = final[1,]
      
      data = final[-1,]
      
      class(data) <- "numeric"
      
      ###############################################################################  
      
      data <- data[, myprobe, drop = FALSE]
    }
    if (length(myprobe) == 0) {
      if (verbose) {
        message(sprintf("probe candidates: 0/%i", size))
      }
      tt <- rep(NA, nrow(data))
      names(tt) <- dimnames(data)[[1]]
      return(list(score = tt, mapping = c(mapped = 0, total = nrow(x)), 
                  probe = cbind(probe = names(gid1), EntrezGene.ID = gid1, 
                                new.probe = names(gid2))))
    }
    if (size == 0 || size > nrow(x)) {
      size <- length(myprobe)
    }
    nix <- 1:size
    myprobe <- myprobe[nix]
    mygid <- mygid[nix]
    mycoef <- mycoef[nix]
    gid1 <- gid1[nix]
    gid2 <- gid2[nix]
    if (!is.na(cutoff)) {
      nix <- abs(mycoef) > cutoff
      myprobe <- myprobe[nix]
      mygid <- mygid[nix]
      mycoef <- mycoef[nix]
      gid1 <- gid1[nix]
      gid2 <- gid2[nix]
    }
    probe.candp <- myprobe[mycoef >= 0]
    probe.candn <- myprobe[mycoef < 0]
    gf <- length(myprobe)
    gf <- c(mapped = gf, total = nrow(x))
    if (verbose) {
      message(sprintf("probe candidates: %i/%i", gf[1], 
                      gf[2]))
    }
    nprobe <- c(probe.candp, probe.candn)
    myw <- c(p = length(probe.candp)/length(nprobe), n = length(probe.candn)/length(nprobe))
    res <- rep(0, nrow(data))
    if (signed) {
      if (length(probe.candp) > 0) {
        res <- myw["p"] * (apply(X = data[, probe.candp, 
                                          drop = FALSE], MARGIN = 1, FUN = sum, na.rm = TRUE)/apply(X = data[, 
                                                                                                             probe.candp, drop = FALSE], MARGIN = 1, FUN = function(x) {
                                                                                                               return(sum(!is.na(x)))
                                                                                                             }))
      }
      if (length(probe.candn) > 0) {
        res <- res - myw["n"] * (apply(X = data[, probe.candn, 
                                                drop = FALSE], MARGIN = 1, FUN = sum, na.rm = TRUE)/apply(X = data[, 
                                                                                                                   probe.candn, drop = FALSE], MARGIN = 1, FUN = function(x) {
                                                                                                                     return(sum(!is.na(x)))
                                                                                                                   }))
      }
    }
    else {
      if (length(probe.candp) > 0) {
        res <- myw["p"] * (apply(X = data[, probe.candp, 
                                          drop = FALSE], MARGIN = 1, FUN = function(x, 
                                                                                    y) {
                                            nix <- is.na(x)
                                            return(sum(x * y, na.rm = TRUE)/sum(y[!nix]))
                                          }, y = abs(mycoef[probe.candp])))
      }
      if (length(probe.candn) > 0) {
        res <- res - myw["n"] * (apply(X = data[, probe.candn, 
                                                drop = FALSE], MARGIN = 1, FUN = function(x, 
                                                                                          y) {
                                                  nix <- is.na(x)
                                                  return(sum(x * y, na.rm = TRUE)/sum(y[!nix]))
                                                }, y = abs(mycoef[probe.candn])))
      }
    }
    return(list(score = res, mapping = gf, probe = cbind(probe = names(gid1), 
                                                         EntrezGene.ID = gid1, new.probe = names(gid2))))
  }
  
  
  
  
  
  
  size = 0
  
  no_grade_1_3 = vector()
  
  i = 1
  c = 1
  k = 1
  
  comparison = numeric()
  ggi_results = list()
  ggi_rankclass = list()
  grades_ggi_genefu = list()
  grades_ggi_rankclass = list()
  probe_ggi = numeric()
  nb_patients = numeric()
  for(i in 1:length(groups)){
    
    
    
    if(i ==17){
      next}
    
    exp = exprs(alldatasets[[i]])
    
    
    data = exp
    
    
    
    
    # Filter healthy people
    
    pheno = pData(alldatasets[[i]])
    
    pos = which(pheno$sample_type%in%"tumor")
    
    pheno = pheno[pos,]
    data = data[,pos]
    grades = pheno$grade
    
    
    
    # Mean normalization of the probes
    rm=rowMeans(data)
    
    data=data-rm
    
    
    
    nb_patients[i] = nrow(pheno)
    # Select platform for the dataset
    
    plat = as.vector(final_summary$biomart[which(final_summary$dataset%in%groups[i])])
    
    plat = plat[which(plat%in%total_annotation$platform)]
    
    # Extra platform information based on maximum number of probes mapped to the annotation provided by Genefu
    plat_added = as.vector(plat_new[which(plat_new$datasets_notplatform%in%groups[i]),2])
    
    plat = c(plat,plat_added)
    plat = plat[!is.na(plat)]
    
    
    if(length(plat)>=1){
      
      # Extract the annotation file from Nick?s approach (Total annotation is the file that gathers all the platforms that I could find)
      NCBI70Affy = total_annotation[which(total_annotation$platform%in%plat),]
      
      
      safe = NCBI70Affy[,c(1:3)]
      
      names(safe) = c("probe","entrez","symbol")
      
      NCBI70Affy = NCBI70[,c(1,4,5)]
      names(NCBI70Affy) =  c("probe","entrez","symbol")
      # We have NICKS. Now we will also combine it with the extra information we may get from merging genefu annot
      
      
      
      var = featureData(alldatasets[[i]])
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      annot = split_annotation_genefu(annot)
      
      annot = annot[,c(1,3)]
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      
      
      names(annot) = c("probe","entrez")
      
      annot$entrez = as.numeric(as.vector(annot$entrez))
      NCBI70Affy$entrez = as.numeric(as.vector(NCBI70Affy$entrez))
      
      annot = annot[which(!is.na(annot$entrez)),]
      NCBI70Affy = NCBI70Affy[which(!is.na(NCBI70Affy$entrez)),]
      
      
      
      
      dat = left_join(NCBI70Affy,annot, by = "entrez")
      
      
      res1 = dat[,c(1:3)]
      res2 = dat[,c(4,2,3)]
      
      names(res1) = c("probe","entrez","symbol")
      names(res2) = c("probe","entrez","symbol")
      
      
      res = rbind(safe, res1,res2)
      
      res = res[complete.cases(res),]
      
      res = res[!duplicated(res),]
      
      
      annot = res
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      
      
      annot = annot[which(annot$entrez%in%NCBI70$EntrezGene.ID),]
      
      annot = annot[,c(1,2)]
      names(annot) = c("probe", "EntrezGene.ID")
      
      annot = annot[complete.cases(annot),]
      
      annot = annot[!duplicated(annot),]
      
      rownames(annot) = annot[,1]
      
      
      # Processing data accordingly so we can run the ggi function of genefu
      
      data = data[which(rownames(data)%in%annot$probe),]
      annot = annot[which(annot$probe%in%rownames(data)),]
      
      
      # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
      # and I will order the new files accordinly.
      
      
      new_sig_ggi = merge(sig.ggi,annot,by.x = "EntrezGene.ID",by.y = "EntrezGene.ID",all.x = T)
      new_sig_ggi = new_sig_ggi[match(rownames(data), new_sig_ggi$probe.y),]
      annot = new_sig_ggi[,c(2,1,10)]
      names(annot) = c("probe","EntrezGene.ID","probe_dataset")
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      
      
      data = data[which(rownames(data)%in%annot$probe_dataset),]
      annot = annot[which(annot$probe_dataset%in%rownames(data)),]
      
      rownames(data) = as.vector(annot$probe)
      
      
      
      if(length(which(!is.na(grades)))==0){
        
        
        
        ggi_results[[i]] = NA
        
        
        
        prueba = data
        
        rowdesc <- rownames(prueba)
        
        sub = annot[which(annot$probe%in%rowdesc),c(1,2)]
        
        prueba = prueba[which(rownames(prueba)%in%as.vector(sub$probe)),]
        rowdesc <- rownames(prueba)
        
        grades_ggi_genefu[[i]] = grades
        
        
        tryCatch(output <- rankclass(prueba, annotation="entrezgene"), error = function(e) { output <<- NULL})
        
        grades_ggi_rankclass[[i]] = as.numeric(as.vector(output))
        
        if(length(output)<1){
          
          ggi_rankclass[[i]] = NA
          
        }else{
          
          rownames(prueba) = as.vector(annot$probe)
          rs1 <- ggi_modified(data=t(prueba), annot=annot, do.mapping=FALSE,hg = output, verbose = FALSE)
          
          if(('1' %!in% output) | ('3' %!in% output)){
            
            no_grade_1_3 = c(no_grade_1_3,groups[i])
          }
          
          ggi_rankclass[[i]] = rs1
          
        }
        
        
        
        
      }else{
        
        
        if(('1' %!in% grades) | ('3' %!in% grades)){
          
          no_grade_1_3 = c(no_grade_1_3,groups[i])
        }
        
        
        comparison[k] = groups[i]
        
        k = k + 1
        
        
        rs_genefu <- ggi_modified(data=t(data), annot=annot, do.mapping=FALSE,hg = grades, verbose = FALSE)
        grades_ggi_genefu[[i]] = as.numeric(as.vector(grades))
        ggi_results[[i]] = rs_genefu
        
        
        # rankclass
        prueba = data
        
        rowdesc <- rownames(prueba)
        
        sub = annot[which(annot$probe%in%rowdesc),c(1,2)]
        
        
        
        tryCatch(output <- rankclass(prueba,annotation="entrezgene"), error = function(e) { output <<- NULL})
        
        grades_ggi_rankclass[[i]] = as.numeric(as.vector(output))
        
        
        # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
        # and I will order the new files accordinly.
        
        rs1 <- ggi_modified(data=t(data), annot=annot, do.mapping=FALSE,hg = output, verbose = FALSE)
        
        ggi_rankclass[[i]] = rs1
      }
      
      
      
      
      
    }else{
      
      
      # We coul not find any platform information for these datasets so we need to rely on the annotation from Genefu
      
      var = featureData(alldatasets[[i]])
      
      
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      # It can happen that some of the annotations are merge with ///. We split them:
      
      annot = split_annotation_genefu(annot)
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      annot[,4] = trimws(annot[,4],"both")
      
      annot$EntrezGene.ID = as.numeric(as.vector(annot$EntrezGene.ID))
      annot = annot[which(!is.na(annot$EntrezGene.ID)),]
      
      # We map this annotation to the original mammaprint annotation file so we can link our genes/probes to Original IDs
      merge1 = merge(NCBI70,annot, by.x="EntrezGene.ID", by.y="EntrezGene.ID", all.x=T)
      
      first = merge1[,c(2,1)]
      names(first) = c("probe","EntrezGene.ID")
      
      second = merge1[,c(10,1)]
      names(second) = c("probe","EntrezGene.ID")
      
      annot = rbind(first,second)
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
      
      
      if(groups[i]%in%names(gpl_datasets)){
        
        pos = which(names(gpl_datasets)%in%groups[i])
        
        annot_extra = gpl_datasets[[pos]]
        
        names(annot_extra) = c("probe","EntrezGene.ID")
        
        annot = rbind(annot,annot_extra)
        
        annot = annot[complete.cases(annot),]
        annot = annot[!duplicated(annot),]
      }
      rownames(annot) = annot[,1]
      
      
      
      # Processing data accordingly so we can run the ggi function of genefu
      
      
      data = data[which(rownames(data)%in%annot$probe),]
      annot = annot[which(annot$probe%in%rownames(data)),]
      
      
      # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
      # and I will order the new files accordinly.
      
      
      new_sig_ggi = merge(sig.ggi,annot,by.x = "EntrezGene.ID",by.y = "EntrezGene.ID",all.x = T)
      new_sig_ggi = new_sig_ggi[match(rownames(data), new_sig_ggi$probe.y),]
      annot = new_sig_ggi[,c(2,1,10)]
      names(annot) = c("probe","EntrezGene.ID","probe_dataset")
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
      
      
      data = data[which(rownames(data)%in%annot$probe_dataset),]
      annot = annot[which(annot$probe_dataset%in%rownames(data)),]
      
      rownames(data) = as.vector(annot$probe)
      
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      
      
      
      # DATASET 17 has no good mappings
      if(nrow(annot)<1){
        
        
        objeto = list(result = NA,ratio_mapped = 0,total_mapped = 0)
        
        ggi_rankclass[[i]] = NA
      }else{
        
        
        if(length(which(!is.na(grades)))==0){
          
          
          
          
          ggi_results[[i]] = NA
          prueba = data
          
          rowdesc <- rownames(prueba)
          
          sub = annot[which(annot$probe%in%rowdesc),c(1,2)]
          
          
          grades_ggi_genefu[[i]] = grades
          
          
          tryCatch(output <- rankclass(prueba, annotation="entrezgene"), error = function(e) { output <<- NULL})
          
          grades_ggi_rankclass[[i]] = as.numeric(as.vector(output))
          
          if(length(output)<1){
            
            ggi_rankclass[[i]] = NA
            
          }else{
            
            
            
            rownames(prueba) = as.vector(annot$probe)
            
            rs1 <- ggi_modified(data=t(prueba), annot=annot, do.mapping=FALSE,hg = output, verbose = FALSE)
            
            ggi_rankclass[[i]] = rs1
            
          }
          
          
          
          
        }else{
          
          
          if(('1' %!in% grades) | ('3' %!in% grades)){
            
            no_grade_1_3 = c(no_grade_1_3,groups[i])
          }
          
          
          comparison[k] = groups[i]
          
          k = k + 1
          
          
          
          rownames(data) = as.vector(annot$probe)
          
          rs_genefu <- ggi_modified(data=t(data), annot=annot, do.mapping=FALSE,hg = grades, verbose = FALSE)
          grades_ggi_genefu[[i]] = as.numeric(as.vector(grades))
          ggi_results[[i]] = rs_genefu
          
          
          # rankclass
          prueba = data
          
          rowdesc <- rownames(prueba)
          
          sub = annot[which(annot$probe%in%rowdesc),c(1,2)]
          
          
          
          tryCatch(output <- rankclass(prueba, annotation="EntrezGene.ID"), error = function(e) { output <<- NULL})
          
          grades_ggi_rankclass[[i]] = as.numeric(as.vector(output))
          rs1 <- ggi_modified(data=t(data), annot=annot, do.mapping=FALSE,hg = output, verbose = FALSE)
          
          ggi_rankclass[[i]] = rs1
        }
        
        
      }
      
      
    }
    
    
    
    print(i)
  }
  
  
  
  
  
  
  
  
  # Patients that do not have grade information.
  positions = sapply(grades_ggi_genefu, function(x) which(is.na(x)==TRUE))
  
  
  # How many do we lose if we dont complete with rankclass grades
  
  missing = sapply(positions, function(x) length(x))
  
  missed_patients = as.data.frame(cbind(groups,missing,nb_patients))
  
  missed_patients[,2] = as.numeric(as.vector(missed_patients[,2]))
  missed_patients[,3] = as.numeric(as.vector(missed_patients[,3]))
  
  missed_patients[,4] = missed_patients[,2]/missed_patients[,3]
  
  
  
  # We will substitute it with the grade information we took from rankclass
  
  a = grades_ggi_genefu
  
  b = grades_ggi_rankclass
  
  c = positions
  
  
  full_grades = Map(function(a,b,c) {a[c]<-b[c]; a},a,b,c)
  
  
  
  
  
  ########################################################################################################################################
  # Now we repeat the GGI analysis. This time we have all the grades for the missing patients
  ########################################################################################################################################
  
  no_grade_1_3 = vector()
  
  i = 1
  c = 1
  k = 1
  
  comparison = numeric()
  ggi_results = list()
  grades_ggi_genefu = list()
  grades_ggi_rankclass = list()
  genes_ggi = numeric()
  for(i in 1:length(groups)){
    
    
    
    if(i ==17){
      next}
    
    exp = exprs(alldatasets[[i]])
    
    
    data = exp
    
    # Filter healthy people
    
    pheno = pData(alldatasets[[i]])
    
    pos = which(pheno$sample_type%in%"tumor")
    
    pheno = pheno[pos,]
    data = data[,pos]
    grades = full_grades[[i]]
    
    
    
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
      NCBI70Affy = total_annotation[which(total_annotation$platform%in%plat),]
      
      
      safe = NCBI70Affy[,c(1:3)]
      
      names(safe) = c("probe","entrez","symbol")
      
      NCBI70Affy = NCBI70[,c(1,4,5)]
      names(NCBI70Affy) =  c("probe","entrez","symbol")
      # We have NICKS. Now we will also combine it with the extra information we may get from merging genefu annot
      
      
      
      var = featureData(alldatasets[[i]])
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      annot = split_annotation_genefu(annot)
      
      annot = annot[,c(1,3)]
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      
      
      names(annot) = c("probe","entrez")
      
      annot$entrez = as.numeric(as.vector(annot$entrez))
      NCBI70Affy$entrez = as.numeric(as.vector(NCBI70Affy$entrez))
      
      annot = annot[which(!is.na(annot$entrez)),]
      NCBI70Affy = NCBI70Affy[which(!is.na(NCBI70Affy$entrez)),]
      
      
      
      
      dat = left_join(NCBI70Affy,annot, by = "entrez")
      
      
      res1 = dat[,c(1:3)]
      res2 = dat[,c(4,2,3)]
      
      names(res1) = c("probe","entrez","symbol")
      names(res2) = c("probe","entrez","symbol")
      
      
      res = rbind(safe, res1,res2)
      
      res = res[complete.cases(res),]
      
      res = res[!duplicated(res),]
      
      
      annot = res
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      
      
      annot = annot[which(annot$entrez%in%NCBI70$EntrezGene.ID),]
      
      annot = annot[,c(1,2)]
      names(annot) = c("probe", "EntrezGene.ID")
      
      annot = annot[complete.cases(annot),]
      
      annot = annot[!duplicated(annot),]
      
      rownames(annot) = annot[,1]
      
      
      # Processing data accordingly so we can run the ggi function of genefu
      
      data = data[which(rownames(data)%in%annot$probe),]
      annot = annot[which(annot$probe%in%rownames(data)),]
      
      
      # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
      # and I will order the new files accordinly.
      
      
      new_sig_ggi = merge(sig.ggi,annot,by.x = "EntrezGene.ID",by.y = "EntrezGene.ID",all.x = T)
      new_sig_ggi = new_sig_ggi[match(rownames(data), new_sig_ggi$probe.y),]
      annot = new_sig_ggi[,c(2,1,10)]
      names(annot) = c("probe","EntrezGene.ID","probe_dataset")
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      
      
      data = data[which(rownames(data)%in%annot$probe_dataset),]
      annot = annot[which(annot$probe_dataset%in%rownames(data)),]
      
      rownames(data) = as.vector(annot$probe)
      
      
      
      
      
      if(length(which(!is.na(grades)))==0){
        
        
        
        ggi_results[[i]] = NA
        prueba = data
        
        rowdesc <- rownames(prueba)
        
        sub = annot[which(annot$probe%in%rowdesc),c(1,2)]
        
        prueba = prueba[which(rownames(prueba)%in%as.vector(sub$probe)),]
        rowdesc <- rownames(prueba)
        
        grades_ggi_genefu[[i]] = grades
        
        
        tryCatch(output <- rankclass(prueba, rowdesc, annotation="entrezgene"), error = function(e) { output <<- NULL})
        
        grades_ggi_rankclass[[i]] = as.numeric(as.vector(output))
        
        if(length(output)<1){
          
          ggi_rankclass[[i]] = NA
          
        }else{
          
          rs1 <- ggi_modified(data=t(prueba), annot=annot, do.mapping=FALSE,hg = output, verbose = FALSE)
          
          if(('1' %!in% output) | ('3' %!in% output)){
            
            no_grade_1_3 = c(no_grade_1_3,groups[i])
          }
          
          ggi_rankclass[[i]] = rs1
          
        }
        
        
        
        
      }else{
        
        
        if(('1' %!in% grades) | ('3' %!in% grades)){
          
          no_grade_1_3 = c(no_grade_1_3,groups[i])
        }
        
        
        comparison[k] = groups[i]
        
        k = k + 1
        
        data = data[which(rownames(data)%in%annot$probe),]
        annot = annot[which(annot$probe%in%rownames(data)),]
        
        
        # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
        # and I will order the new files accordinly.
        
        
        new_sig_ggi = merge(sig.ggi,annot,by.x = "EntrezGene.ID",by.y = "EntrezGene.ID",all.x = T)
        new_sig_ggi = new_sig_ggi[match(rownames(data), new_sig_ggi$probe.y),]
        annot = new_sig_ggi[,c(2,1,10)]
        names(annot) = c("probe","EntrezGene.ID","probe_dataset")
        
        annot = annot[complete.cases(annot),]
        annot = annot[!duplicated(annot),]
        
        
        data = data[which(rownames(data)%in%annot$probe_dataset),]
        annot = annot[which(annot$probe_dataset%in%rownames(data)),]
        
        
        
        ok = unique( rownames(data))
        
        
        rs_genefu <- ggi_modified(data=t(data), annot=annot, do.mapping=FALSE,hg = grades, verbose = FALSE)
        grades_ggi_genefu[[i]] = as.numeric(as.vector(grades))
        ggi_results[[i]] = rs_genefu
        
        
        
      }
      
      
      
    }else{
      
      
      # We could not find any platform information for these datasets so we need to rely on the annotation from Genefu
      
      var = featureData(alldatasets[[i]])
      
      
      annot = pData(var)
      names(annot) = c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
      
      
      # It can happen that some of the annotations are merge with ///. We split them:
      
      annot = split_annotation_genefu(annot)
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      annot[,3] = trimws(annot[,3],"both")
      annot[,4] = trimws(annot[,4],"both")
      
      annot$EntrezGene.ID = as.numeric(as.vector(annot$EntrezGene.ID))
      annot = annot[which(!is.na(annot$EntrezGene.ID)),]
      
      # We map this annotation to the original mammaprint annotation file so we can link our genes/probes to Original IDs
      merge1 = merge(NCBI70,annot, by.x="EntrezGene.ID", by.y="EntrezGene.ID", all.x=T)
      
      first = merge1[,c(2,1)]
      names(first) = c("probe","EntrezGene.ID")
      
      second = merge1[,c(10,1)]
      names(second) = c("probe","EntrezGene.ID")
      
      annot = rbind(first,second)
      
      annot = annot[complete.cases(annot),]
      annot = annot[!duplicated(annot),]
      
      
      if(groups[i]%in%names(gpl_datasets)){
        
        pos = which(names(gpl_datasets)%in%groups[i])
        
        annot_extra = gpl_datasets[[pos]]
        
        names(annot_extra) = c("probe","EntrezGene.ID")
        
        annot = rbind(annot,annot_extra)
        
        annot = annot[complete.cases(annot),]
        annot = annot[!duplicated(annot),]
      }
      rownames(annot) = annot[,1]
      
      annot[,1] = trimws(annot[,1],"both")
      annot[,2] = trimws(annot[,2],"both")
      
      
      # DATASET 17 has no good mappings
      if(nrow(annot)<1){
        
        
        objeto = list(result = NA,ratio_mapped = 0,total_mapped = 0)
        
        ggi_rankclass[[i]] = NA
      }else{
        
        
        if(length(which(!is.na(grades)))==0){
          
          
          
          
          ggi_results[[i]] = NA
          prueba = data
          
          rowdesc <- rownames(prueba)
          
          sub = annot[which(annot$probe%in%rowdesc),c(1,2)]
          
          prueba = prueba[which(rownames(prueba)%in%sub$probe),]
          rowdesc <- rownames(prueba)
          
          grades_ggi_genefu[[i]] = grades
          
          
          tryCatch(output <- rankclass(prueba, rowdesc, annotation="entrezgene"), error = function(e) { output <<- NULL})
          
          grades_ggi_rankclass[[i]] = as.numeric(as.vector(output))
          
          if(length(output)<1){
            
            ggi_rankclass[[i]] = NA
            
          }else{
            
            
            
            # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
            # and I will order the new files accordinly.
            
            
            new_sig_ggi = merge(sig.ggi,annot,by.x = "EntrezGene.ID",by.y = "EntrezGene.ID",all.x = T)
            new_sig_ggi = new_sig_ggi[match(rownames(data), new_sig_ggi$probe.y),]
            annot = new_sig_ggi[,c(2,1,10)]
            names(annot) = c("probe","EntrezGene.ID","probe_dataset")
            
            annot[,1] = trimws(annot[,1],"both")
            annot[,2] = trimws(annot[,2],"both")
            annot[,3] = trimws(annot[,3],"both")
            
            annot = annot[complete.cases(annot),]
            annot = annot[!duplicated(annot),]
            
            
            rownames(prueba) = as.vector(annot$probe)
            
            rs1 <- ggi_modified(data=t(prueba), annot=annot, do.mapping=FALSE,hg = output, verbose = FALSE)
            
            
            
            ggi_rankclass[[i]] = rs1
            
          }
          
          
          
          
        }else{
          
          
          if(('1' %!in% grades) | ('3' %!in% grades)){
            
            no_grade_1_3 = c(no_grade_1_3,groups[i])
          }
          
          
          comparison[k] = groups[i]
          
          k = k + 1
          
          data = data[which(rownames(data)%in%annot$probe),]
          
          annot = annot[which(annot$probe%in%rownames(data)),]
          
          
          # Genefu?s ggi requires that if we DONT use their mapping system WE MUST have the probes as in sig.ggi. Therefore I will map it
          # and I will order the new files accordinly.
          
          
          new_sig_ggi = merge(sig.ggi,annot,by.x = "EntrezGene.ID",by.y = "EntrezGene.ID",all.x = T)
          new_sig_ggi = new_sig_ggi[match(rownames(data), new_sig_ggi$probe.y),]
          annot = new_sig_ggi[,c(2,1,10)]
          names(annot) = c("probe","EntrezGene.ID","probe_dataset")
          
          annot = annot[complete.cases(annot),]
          annot = annot[!duplicated(annot),]
          
          
          data = data[which(rownames(data)%in%annot$probe_dataset),]
          annot = annot[which(annot$probe_dataset%in%rownames(data)),]
          
          rownames(data) = as.vector(annot$probe)
          
          rs_genefu <- ggi_modified(data=t(data), annot=annot, do.mapping=FALSE,hg = grades, verbose = FALSE)
          grades_ggi_genefu[[i]] = as.numeric(as.vector(grades))
          ggi_results[[i]] = rs_genefu
          
        }
        
        
      }
      
      
    }
    
    genes_ggi[i] = length(unique(annot$EntrezGene.ID))/length(unique(NCBI70$EntrezGene.ID))
    
    
    print(i)
  }
  
  return(ggi_results)
  
}


ggi_results = ggi_function()



save(ggi_results,file = "data_bitbucket/final_results/ggi_results")





