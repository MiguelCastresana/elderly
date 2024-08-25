

# LOAD 55-65 years old patients

sub55 = read.delim("C:/Users/migcas/OneDrive - Karolinska Institutet/Desktop/Postdoc/data_bitbucket/final_results/all_results_55_65", header = T)


# We select ER positive and LN negative
sub55 = sub55[which(sub55$er%in%"positive" & sub55$N==0),]


# LOAD >70 years old patients
sub70 = read.delim("C:/Users/migcas/OneDrive - Karolinska Institutet/Desktop/Postdoc/data_bitbucket/final_results/all_results", header = T)

# We select ER positive and LN negative
sub70 = sub70[which(sub70$er%in%"positive" & sub70$N==0),]


# GGI
vec1 = c(sub55$ggi_risk,sub70$ggi_risk)

vec2 = c(rep("0",length(sub55$ggi_risk)),rep("1",length(sub70$ggi_risk)))
vec2 = as.numeric(as.vector(vec2))


chisq.test(table(vec1,vec2))


# 70-GENE
vec1 = c(sub55$mammaprint.prognosis_binary,sub70$mammaprint.prognosis_binary)


vec2 = c(rep("55-65",length(sub55$mammaprint.prognosis_binary)),rep(">70",length(sub70$mammaprint.prognosis_binary)))

chisq.test(table(vec1,vec2))


# ONCOTYPE
vec1 = c(sub55$oncotype_class,sub70$oncotype_class)

vec2 = c(rep("55-65",length(sub55$oncotype_class)),rep(">70",length(sub70$oncotype_class)))



chisq.test(table(vec1,vec2))



# RORP
vec1 = c(sub55$ROR.P.Group..Subtype...Proliferation.,sub70$ROR.P.Group..Subtype...Proliferation.)

vec2 = c(rep("55-65",length(sub55$ROR.P.Group..Subtype...Proliferation.)),rep(">70",length(sub70$ROR.P.Group..Subtype...Proliferation.)))



chisq.test(table(vec1,vec2))


# CELL CYCLE

vec1 = c(sub55$cell_cycle_risk,sub70$cell_cycle_risk)

vec2 = c(rep("55-65",length(sub55$cell_cycle_risk)),rep(">70",length(sub70$cell_cycle_risk)))



chisq.test(table(vec1,vec2))



# PAM50

vec1 = c(sub55$Call,sub70$Call)

vec2 = c(rep("55-65",length(sub55$Call)),rep(">70",length(sub70$Call)))






chisq.test(table(vec1,vec2))

