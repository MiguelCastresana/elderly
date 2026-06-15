


## CCS signature


# percentage of genes covered

load("C:/Users/migcas/OneDrive - Karolinska Institutet/Desktop/Postdoc/data_final_final/final_summary_ALL")

datasets = final_summary$dataset

c = 1
percentage_genes = numeric()
i = 1
for(i in 1:length(datasets)){
  
  # pos = which(datasets%in%groups[i])
  # if(i==17){next}
  metabric = data_clean[[i]][[1]]
  
  
  percentage_genes[c] = length(ccs_genes$hgnc_symbol[which(ccs_genes$hgnc_symbol%in%rownames(metabric))])/nrow(ccs_genes)
  
  c = c + 1
  
  print(i)
}

dat = as.data.frame(cbind(groups,percentage_genes))

colnames(dat)[2] = "mapped_cell_cycle_genes"

resumen = dat


# pam50


files = list.files(paste("data_bitbucket/PAM50_ROR_P/PAM50_results_review/",sep=""))

mapped_pam50 = numeric()
i = 1
for(i in 1:length(files)){
  
  a = read.delim(paste("data_bitbucket/PAM50_ROR_P/PAM50_results_review/",files[i],"/MedCenForPAM50_Genes.txt",sep=""))
  mapped_pam50[i] = nrow(a)/50
  
}

dat = as.data.frame(cbind(files,mapped_pam50))

colnames(dat)[1] = "groups"


resumen = left_join(resumen,dat,by = "groups")









# oncotype

mapped_oncotype_miguel = sapply(oncotype_results, function(x) x[[3]])

run_oncotype_miguel <- ifelse(mapped_oncotype_miguel == 16, 1, 0)


dat = as.data.frame(cbind(groups,run_oncotype_miguel,mapped_oncotype_miguel))

dat = dat[,1:2]

resumen = left_join(resumen,dat,by = "groups")


# mammaprint


safe = sapply(mammaprint_results, function(x) x[[2]])

ratio_mapped_mammaprint_miguel = sapply(mammaprint_results, `[[`, 2)

mapped_mammaprint_miguel = sapply(mammaprint_results, `[[`, 3)

dat = as.data.frame(cbind(groups,ratio_mapped_mammaprint_miguel,mapped_mammaprint_miguel))

dat = dat[,1:2]


resumen = left_join(resumen,dat,by = "groups")


# ggi

load("data_bitbucket/ggi/ggi_genes")

dat = as.data.frame(cbind(groups,genes_ggi))

resumen = left_join(resumen,dat,by = "groups")



resumen = resumen[which(resumen$groups%in%unique(jajajaja$dataset)),]


resumen[, -1] <- lapply(resumen[, -1], as.numeric)



# Compute the mean of all columns except the first one

mean_values <- sapply(resumen[, -1], mean, na.rm = TRUE)

# Compute the median of all columns except the first one
median_values <- sapply(resumen[, -1], median, na.rm = TRUE)

