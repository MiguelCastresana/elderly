

cellcycle_function = function(){
  
  
 
  
  
  
  metabric = data_clean[[38]][[1]]
  # write.table(metabric,"data_bitbucket/Cell_cycle/metabric_expression_set",quote = FALSE, row.names = FALSE, sep = "\t")
  
  # We compute the Metabric CCS cutoff values
  
  metabric_ccs = metabric[which(rownames(metabric)%in%ccs_genes$hgnc_symbol),]
  
  metabric_ccs_values = as.numeric(as.vector(colSums(metabric_ccs,na.rm = T)))
  
  metabric_ccs_values = (metabric_ccs_values - min(metabric_ccs_values)) / (max(metabric_ccs_values) - min(metabric_ccs_values))
  
  metabric_ccs_values = as.data.frame(metabric_ccs_values)
  
  
  
  tertiles = quantile(metabric_ccs_values$metabric_ccs_values, c(.33333, .66666, 1))
  
  cutoff_1 = tertiles[1] # 0.3720803
  
  cutoff_2 = tertiles[2] # 0.5175723
  
  
  
  # Now we run cell cycle in each dataset using metabric cutoffs
  
  
  i = 1
  risk_list = list()
  for(i in 1:length(alldatasets)){
    
    
    if(i==17){next}
    metabric = data_clean[[i]][[1]]
    
    
    ccs_genes$hgnc_symbol[which(ccs_genes$hgnc_symbol%in%rownames(metabric))]
    
    
    
    metabric_ccs = metabric[which(rownames(metabric)%in%ccs_genes$hgnc_symbol),]
    
    metabric_ccs_values = as.numeric(as.vector(colSums(metabric_ccs,na.rm = T)))
    
    
    
    metabric_ccs_values = (metabric_ccs_values - min(metabric_ccs_values)) / (max(metabric_ccs_values) - min(metabric_ccs_values))
    
    pos1_1 = which(metabric_ccs_values<=cutoff_1)
    pos1_2 = which(metabric_ccs_values>cutoff_1 & metabric_ccs_values<=cutoff_2)
    
    pos1_3 = which(metabric_ccs_values>cutoff_2)
    
    
    risk = rep("low",length(metabric_ccs_values))
    
    risk[pos1_1] = "Low"
    risk[pos1_2] = "Intermediate"
    risk[pos1_3] = "High"
    
    
    a = as.data.frame(cbind(risk,colnames(data_clean[[i]][[1]]),rep(groups[i],length(risk))))
    
    risk_list[[i]] = a
    
    print(i)
  }
  
  
  cellcycle = do.call(rbind.data.frame, risk_list)
  
  names(cellcycle) = c("cell_cycle_risk","sample_name","dataset")

  
  
  return(cellcycle)

}


cell_cycle_results = cellcycle_function()



save(cell_cycle_results,file = "data_bitbucket/final_results/cell_cycle_results")
