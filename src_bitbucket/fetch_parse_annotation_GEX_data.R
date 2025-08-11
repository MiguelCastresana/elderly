## =========================
## 0) Setup & paths
## =========================
# install.packages("BiocManager")
# BiocManager::install(c(
#   "AnnotationDbi","hgu133a.db","hgu133b.db","hgu95av2.db","hgu133plus2.db",
#   "u133x3p.db","hgug4110b.db","hgug4112a.db","hugene10stprobeset.db",
#   "illuminaHumanv2.db","illuminaHumanv3.db",
#   "genefu","MetaGxBreast"
# ))

library(AnnotationDbi)
library(dplyr)
library(stringr)
library(genefu)
library(MetaGxBreast)
library(ggplot2)
library(tidyr)

`%!in%` <- function(x,y)!('%in%'(x,y))

# Base directory
base_dir <- "~/data_bitbucket/"

paths <- list(
  annot_dir   = file.path(base_dir, "annotation_package_ALL_genes"),
  gpl_dir     = file.path(base_dir, "GPL"),
  merged_dir  = file.path(base_dir, "merged_data"),
  summary     = file.path(base_dir, "SUMMARY_datasets"),
  total_tsv   = file.path(base_dir, "total_annotation_all_genes.tsv"),
  gpl_cache   = file.path(base_dir, "gpl_annotation_missing_datasets_all_genes")
)


dir.create(paths$annot_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(paths$merged_dir, showWarnings = FALSE, recursive = TRUE)

## =========================
## 1) Build per-platform annotations (probe -> {symbol, entrez})
## =========================
build_platform_annot <- function(db_pkg, platform_name, out_dir) {
  k <- tryCatch(keys(db_pkg, keytype = "PROBEID"), error = function(e) character(0))
  if (!length(k)) return(NULL)
  
  cols <- c("SYMBOL","ENTREZID","ENSEMBL","GENENAME","CHR","MAP","CYTOBAND")
  cols <- cols[cols %in% columns(db_pkg)]
  
  ann <- AnnotationDbi::select(db_pkg, keys = k, columns = cols, keytype = "PROBEID")
  if (!nrow(ann)) return(NULL)
  
  names(ann)[names(ann)=="PROBEID"]  <- "probe"
  names(ann)[names(ann)=="SYMBOL"]   <- "symbol"
  names(ann)[names(ann)=="ENTREZID"] <- "entrez"
  ann <- ann[, intersect(c("probe","symbol","entrez"), names(ann)), drop = FALSE]
  ann <- unique(ann)
  if ("entrez" %in% names(ann)) ann$entrez <- suppressWarnings(as.character(ann$entrez))
  ann$platform <- platform_name
  
  saveRDS(ann, file = file.path(out_dir, paste0(platform_name, ".rds")))
  ann
}

safeload <- function(pkg) { suppressPackageStartupMessages(require(pkg, character.only = TRUE)) }

# Load DBs you need & build
db_specs <- list(
  list(pkg="hgu133a.db",            name="affy_hg_u133a"),
  list(pkg="hgu133b.db",            name="affy_hg_u133b"),
  list(pkg="hgu95av2.db",           name="affy_hg_u95av2"),
  list(pkg="hgu133plus2.db",        name="affy_hg_u133_plus_2"),
  list(pkg="u133x3p.db",            name="affy_u133_x3p"),
  list(pkg="hgug4110b.db",          name="agilent_human_1A"),
  list(pkg="hgug4112a.db",          name="agilent_human_genome_4x44k"),
  list(pkg="hugene10stprobeset.db", name="affy_hugene_1_0_st_v1"),
  list(pkg="illuminaHumanv2.db",    name="illumina_humanwg_6_v2"),
  list(pkg="illuminaHumanv3.db",    name="illumina_humanht_12_v3")
)

per_platform <- list()
for (sp in db_specs) {
  if (safeload(sp$pkg)) {
    db_obj <- get(sp$pkg)
    per_platform[[sp$name]] <- build_platform_annot(db_obj, sp$name, paths$annot_dir)
  }
}

# Merge all per-platform annotations and write TSV
ann_files <- list.files(paths$annot_dir, pattern="\\.rds$", full.names = TRUE)
all_ann <- do.call(rbind, lapply(ann_files, readRDS))
all_ann <- all_ann %>% filter(!is.na(symbol), symbol != "") %>% distinct()
write.table(all_ann, file = paths$total_tsv, sep="\t", quote=FALSE, row.names=FALSE)

## =========================
## 2) Helper split functions for weird annotations (kept from your code)
## =========================
split_annotation_genefu <- function(annot){
  mm <- apply(annot,1,paste,collapse=" ")
  pos <- grep("///", mm)
  if(length(pos)>0){
    add <- lapply(seq_along(pos), function(k){
      i <- pos[k]
      separar <- lapply(1:ncol(annot), function(j){
        out <- str_split(annot[i,j], "///")
        clean <- lapply(out, function(x) x[x!=""])
        if(length(clean[[1]])<1) clean[[1]] <- NA
        clean
      })
      l <- sapply(separar, function(cl) sapply(cl, length))
      l <- as.numeric(l)
      l_ordered <- l[order(-l)]
      data.frame(
        first_col  = rep(as.vector(annot[i,1]), (max(l)*l_ordered[2])/l[1]),
        second_col = rep(unlist(separar[[2]]),(max(l)*l_ordered[2])/l[2]),
        third_col  = rep(unlist(separar[[3]]),(max(l)*l_ordered[2])/l[3]),
        fourth_col = rep(annot[i,4],(max(l)*l_ordered[2])/l[4]),
        check.names = FALSE
      )
    })
    add <- do.call(rbind.data.frame, add)
    names(add) <- names(annot)
    add <- add[complete.cases(add),]
    add <- add[!(is.na(add[,2])|add[,2]==""),]
    add <- add[!(is.na(add[,3])|add[,3]==""),]
    annot <- annot[-pos,]
    annot <- rbind(annot, add)
  }
  colnames(annot)[1] <- "probe"
  annot
}

split_annotation_gpl <- function(annot){
  mm <- apply(annot,1,paste,collapse=" ")
  pos <- grep("///", mm)
  if(length(pos)>0){
    add <- lapply(seq_along(pos), function(k){
      i <- pos[k]
      separar <- lapply(1:ncol(annot), function(j) str_split(annot[i,j], "///"))
      l <- sapply(separar, function(s) length(s[[1]]))
      data.frame(first_col=rep(as.vector(annot[i,1]),max(l)),
                 second_col=unlist(separar[[2]]),
                 check.names = FALSE)
    })
    add <- do.call(rbind.data.frame, add)
    names(add) <- names(annot)
    add <- add[complete.cases(add),]
    add <- add[!(is.na(add[,2])|add[,2]==""),]
    annot <- annot[-pos,]
    annot <- rbind(annot, add)
  }
  colnames(annot)[1] <- "probe"
  annot
}

split_annotation_gpl_comma <- function(annot){
  mm <- apply(annot,1,paste,collapse=" ")
  pos <- grep(",", mm)
  if(length(pos)>0){
    add <- lapply(seq_along(pos), function(k){
      i <- pos[k]
      separar <- lapply(1:ncol(annot), function(j) str_split(annot[i,j], ","))
      l <- sapply(separar, function(s) length(s[[1]]))
      data.frame(first_col=rep(as.vector(annot[i,1]),max(l)),
                 second_col=unlist(separar[[2]]),
                 check.names = FALSE)
    })
    add <- do.call(rbind.data.frame, add)
    names(add) <- names(annot)
    add <- add[complete.cases(add),]
    add <- add[!(is.na(add[,2])|add[,2]==""),]
    annot <- annot[-pos,]
    annot <- rbind(annot, add)
  }
  colnames(annot)[1] <- "probe"
  annot
}

## =========================
## 3) Read summary & total annotation
## =========================
final_summary <- read.delim(paths$summary, stringsAsFactors = FALSE)
total_annotation <- read.delim(paths$total_tsv, stringsAsFactors = FALSE)

## =========================
## 4) Build GPL rescue annotations (only for odd GPLs)
## =========================
# REQUIRE: `alldatasets` & `groups` set up (MetaGxBreast objects)
# Build `alldatasets` from MetaGxBreast (your original approach):
# These come from your environment:
# - esetsAndDups$esets
# - extra$esets
groups <- c(names(esetsAndDups$esets), names(extra$esets))
alldatasets <- list()
c <- 1
for (i in seq_along(groups)) {
  if (i < 38) {
    alldatasets <- c(alldatasets, esetsAndDups$esets[[i]])
  } else {
    alldatasets <- c(alldatasets, extra$esets[[c]]); c <- c + 1
  }
}

gpl_files <- list.files(paths$gpl_dir, full.names = TRUE)
gpl_datasets <- vector("list", length(gpl_files))
names(gpl_datasets) <- sub("\\.txt.*","", basename(gpl_files))

for (i in seq_along(gpl_files)) {
  annot <- read.delim(gpl_files[i], stringsAsFactors = FALSE)
  g <- names(gpl_datasets)[i]
  pos <- which(groups %in% g)
  if (!length(pos)) { gpl_datasets[[i]] <- NA; next }
  exp <- exprs(alldatasets[[pos]])
  var <- featureData(alldatasets[[pos]])
  annot_genefu <- pData(var)
  names(annot_genefu) <- c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
  
  if (ncol(annot) >= 3 && colnames(annot)[2] == "Gene.title") {
    annot <- annot[,c(1,3)] %>% filter(complete.cases(.)) %>% filter(!is.na(.[,2]) & .[,2]!="")
    annot <- split_annotation_gpl(annot)
    annot <- annot %>% distinct()
    gpl_datasets[[i]] <- annot
    
  } else if (ncol(annot) >= 3 && colnames(annot)[2] == "Comment.OLIGO_ID.") {
    annot <- annot[,c(2,3)] %>% filter(complete.cases(.)) %>% filter(!is.na(.[,2]) & .[,2]!="")
    annot <- split_annotation_gpl(annot)
    names(annot) <- c("probe","symbol")
    annot <- annot %>% distinct()
    gpl_datasets[[i]] <- annot
    
  } else if ("Search_key" %in% colnames(annot)) {
    annot[,1] <- paste0("probe_", annot$ID)
    annot <- annot[annot$ID %in% rownames(exp), c(1,7)]
    names(annot) <- c("probe","symbol")
    gpl_datasets[[i]] <- annot
    
  } else if ("PLATE" %in% colnames(annot)) {
    annot <- annot[which(annot[,1] %in% rownames(exp)), c(1,7)]
    annot <- split_annotation_gpl_comma(annot)
    gpl_datasets[[i]] <- annot
    
  } else {
    gpl_datasets[[i]] <- NA
  }
}

save(gpl_datasets, file = paths$gpl_cache)

## =========================
## 5) Infer missing platforms by majority mapping of probes
## =========================
datasets_notplatform <- final_summary$dataset[which(final_summary$biomart %!in% total_annotation$platform)]
positions <- match(datasets_notplatform, groups)

plat_new <- character(length(datasets_notplatform))
for (i in seq_along(datasets_notplatform)) {
  if (is.na(positions[i])) { plat_new[i] <- NA; next }
  var <- if (positions[i] >= 38) featureData(extra$esets[[positions[i]-37]]) else featureData(esetsAndDups$esets[[positions[i]]])
  annot <- pData(var)
  names(annot) <- c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
  annot <- split_annotation_genefu(annot)
  dat <- total_annotation[ total_annotation$probe %in% annot$probe, ]
  if (nrow(dat)) {
    tab <- table(dat$platform)
    plat_new[i] <- names(tab)[which.max(tab)]
  } else {
    plat_new[i] <- NA
  }
}
plat_new <- data.frame(datasets_notplatform, plat_new, stringsAsFactors = FALSE)

## =========================
## 6) Collapse probes -> genes per dataset and merge all
## =========================
data_clean <- vector("list", length(groups))

for (i in seq_along(groups)) {
  if (i == 17) { data_clean[[i]] <- NA; next }  # your original exclusion
  eset <- alldatasets[[i]]
  exp  <- exprs(eset)
  
  # keep tumor only
  pheno <- pData(eset)
  pos   <- which(pheno$sample_type %in% "tumor")
  pheno <- pheno[pos, , drop=FALSE]
  data  <- exp[, pos, drop=FALSE]
  
  # mean-centering per probe
  data <- sweep(data, 1, rowMeans(data, na.rm = TRUE), FUN = "-")
  
  # preferred platform(s)
  plat <- as.vector(final_summary$biomart[final_summary$dataset %in% groups[i]])
  plat <- plat[plat %in% total_annotation$platform]
  
  # add inferred platform if needed
  plat_added <- as.vector(plat_new$plat_new[plat_new$datasets_notplatform %in% groups[i]])
  plat <- unique(na.omit(c(plat, plat_added)))
  
  meann <- function(v) mean(v, na.rm = TRUE)
  
  if (length(plat) >= 1) {
    # Use our master mapping for those platforms
    NCBI70 <- total_annotation %>% filter(platform %in% plat) %>% select(probe, symbol, entrez) %>% distinct()
    NCBI70 <- NCBI70[NCBI70$probe %in% rownames(data), , drop=FALSE]
    
    # add genefu annotation too, then union
    var <- featureData(eset)
    annot <- pData(var)
    names(annot) <- c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
    annot <- split_annotation_genefu(annot)[, c(1,2,3)]
    names(annot) <- c("probe","symbol","entrez")
    
    res <- bind_rows(NCBI70, annot) %>% distinct()
    res <- res[res$probe %in% rownames(data), , drop=FALSE]
    res$probe  <- str_trim(res$probe)
    res$symbol <- str_trim(res$symbol)
    res <- res %>% filter(complete.cases(.)) %>% distinct()
    
    data2 <- data[rownames(data) %in% res$probe, , drop=FALSE]
    con   <- merge(data2, res, by.x = "row.names", by.y = "probe", all.x = TRUE)
    final <- aggregate(con, list(con$symbol), meann)
    values <- final[, -c(2, ncol(final), (ncol(final)-1))]
    rownames(values) <- values[,"Group.1"]
    values <- values[, -1, drop=FALSE]
    
  } else {
    # fallback: genefu annotation + GPL rescue if available
    var <- featureData(eset)
    annot <- pData(var)
    names(annot) <- c("probe","NCBI.gene.symbol","EntrezGene.ID","best_probe")
    annot <- split_annotation_genefu(annot)
    annot$EntrezGene.ID <- suppressWarnings(as.numeric(as.character(annot$EntrezGene.ID)))
    annot <- annot[!is.na(annot$EntrezGene.ID), c(1,2)]
    names(annot) <- c("probe","symbol")
    
    if (groups[i] %in% names(gpl_datasets)) {
      extra_ann <- gpl_datasets[[groups[i]]]
      if (!is.null(extra_ann) && !is.na(extra_ann)[1]) {
        names(extra_ann) <- c("probe","symbol")
        annot <- bind_rows(annot, extra_ann) %>% filter(complete.cases(.)) %>% distinct()
      }
    }
    annot <- annot[annot$probe %in% rownames(exp), , drop=FALSE]
    data2 <- data[rownames(data) %in% annot$probe, , drop=FALSE]
    con   <- merge(data2, annot, by.x="row.names", by.y="probe", all.x=TRUE)
    final <- aggregate(con, list(con$symbol), meann)
    values <- final[, -c(2, ncol(final))]
    rownames(values) <- values[,"Group.1"]
    values <- values[, -1, drop=FALSE]
  }
  
  pheno <- pheno[, c("sample_name","er"), drop=FALSE]
  data_clean[[i]] <- list(values=values, annot=res, pheno=pheno)
  cat("Processed:", i, groups[i], "\n")
}

save(data_clean, file = file.path(paths$merged_dir, "dataclean.rds"))


