


# sub = read.delim("C:/Users/Miguel/OneDrive - Karolinska Institutet/Desktop/Postdoc/data_bitbucket/final_results/all_results", header = T)


sub = results_final

sub$RS <- ifelse(sub$RS == FALSE, 0, sub$RS)

sub = sub %>%
  mutate(pgr=replace(pgr, pgr%in%"+", "positive")) %>%
  as.data.frame()


sub = sub %>%
  mutate(pgr=replace(pgr, pgr%in%"-", "negative")) %>%
  as.data.frame()


sub = sub %>%
  mutate(her2=replace(her2, her2%in%"+", "positive")) %>%
  as.data.frame()


sub = sub %>%
  mutate(her2=replace(her2, her2%in%"-", "negative")) %>%
  as.data.frame()


sub$cell_cycle_risk <- factor(sub$cell_cycle_risk, levels = c("Low", "Intermediate","High"))


trial = sub

Type<-as.factor(trial$unified_survival_status)
trial$unified_survival_status = Type

trial$unified_survival_status = as.vector(unclass(Type))


Type<-as.factor(trial$grade)
trial$grade = Type

trial$unified_survival_days = as.numeric(as.vector(trial$unified_survival_days))


# Censored data at 10 years


trial = mutate(trial,
               unified_survival_days = ifelse(unified_survival_days >= 10*365.25, 10*365.25, unified_survival_days),
               unified_survival_status = ifelse(unified_survival_days >= 10*365.25, 1, unified_survival_status))



trial <- within(trial, ROR.P.Group..Subtype...Proliferation.[ROR.P.Group..Subtype...Proliferation. == 'low'] <- 'Low')

trial <- within(trial, ROR.P.Group..Subtype...Proliferation.[ROR.P.Group..Subtype...Proliferation. == 'med'] <- 'Intermediate')

trial <- within(trial, ROR.P.Group..Subtype...Proliferation.[ROR.P.Group..Subtype...Proliferation. == 'high'] <- 'High')

trial$ROR.P.Group..Subtype...Proliferation. <- factor(trial$ROR.P.Group..Subtype...Proliferation., levels = c("Low", "Intermediate","High"))

trial <- within(trial, oncotype_class[oncotype_class == 'low'] <- 'Low')

trial <- within(trial, oncotype_class[oncotype_class == 'intermediate'] <- 'Intermediate')

trial <- within(trial, oncotype_class[oncotype_class == 'high'] <- 'High')

trial$oncotype_class<-as.factor(trial$oncotype_class)
trial$oncotype_class <- factor(trial$oncotype_class, levels = c("Low", "Intermediate","High"))

trial$treatment = as.factor(trial$treatment)


# trial = trial[which(trial$er%in%"positive" & trial$N==1),]

##########################             GGI            ########################## 



fit = survfit(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk, data = trial)



ggsurv_ggi = ggsurvplot(fit,   
                        size = 1.2,
                        data = trial,  
                        censor = FALSE,
                        palette =
                          c("blue", "red"),
                        risk.table = TRUE,       # show risk table.
                        pval = TRUE,
                        pval.method=TRUE,# show p-value of log-rank test.
                        pval.size = 7,
                        pval.method.size = 7,
                        conf.int = FALSE,         # show confidence intervals for 
                        xscale = "d_y",
                        break.time.by=365.25*2,
                        # point estimates of survival curves.
                        # survival estimates.
                        ylab = "Survival probability",
                        xlab = "Time in years",   # customize X axis label.
                        ggtheme = theme_classic2(base_size=16), # customize plot and risk table with a theme.
                        risk.table.y.text.col = T, # colour risk table text annotations.
                        risk.table.fontsize = 5,
                        risk.table.y.text = FALSE, # show bars instead of names in text annotations
                        risk.table.col = "strata",
                        legend.labs =
                          c("GG1", "GG3")
)



plotggi = ggsurv_ggi$plot +theme(legend.position = c(0.87, 0.25))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))+theme(legend.title = element_text(size=18),
                                                                                 legend.text = element_text(size=18))+guides(color=guide_legend("Grade"))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))


tableggi = ggsurv_ggi$table + theme(legend.position = "none", axis.text=element_text(size=18), axis.title=element_text(size=18))+
  theme(axis.text.x = element_text(colour="black")) +ylab ("Risk")

ggsurv_ggi_final = ggarrange(plotggi, tableggi, heights = c(1.5, 0.8),widths = c(0.3,0.3),
                             ncol = 1, nrow = 2, align = "v")





ggsave("data_bitbucket/final_results/plots/GGI.pdf",plotggi, width=6, height=5,dpi = 300)

ggsave("data_bitbucket/final_results/plots/GGI_table.pdf",tableggi, width=6, height=2,dpi = 300)



##########################             MAMMAPRINT            ########################## 

a = sub[,c("endpoint","unified_survival_status","unified_survival_days","mammaprint.prognosis_binary")]


fit = survfit(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary, data = trial)


ggsurv_mammaprint = ggsurvplot(fit,   
                        size = 1.2,
                        data = trial,  
                        censor = FALSE,
                        palette =
                          c("blue", "red"),
                        risk.table = TRUE,       # show risk table.
                        pval = TRUE,
                        pval.method=TRUE,
                        pval.size = 7,
                        pval.method.size = 7,# show p-value of log-rank test.
                        conf.int = FALSE,         # show confidence intervals for 
                        xscale = "d_y",
                        break.time.by=365.25*2,
                        # point estimates of survival curves.
                        # survival estimates.
                        ylab = "Survival probability",
                        xlab = "Time in years",   # customize X axis label.
                        ggtheme = theme_classic2(base_size=16),
                        risk.table.fontsize = 5,# customize plot and risk table with a theme.
                        risk.table.y.text.col = T, # colour risk table text annotations.
                        risk.table.y.text = FALSE, # show bars instead of names in text annotations
                        risk.table.col = "strata",
                        legend.labs =
                          c("Low", "High")
)




plotmammaprint = ggsurv_mammaprint$plot +theme(legend.position = c(0.87, 0.25))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))+theme(legend.title = element_text(size=18),
                                                                                 legend.text = element_text(size=18))+guides(color=guide_legend("Risk"))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))


tablemammaprint = ggsurv_mammaprint$table + theme(legend.position = "none", axis.text=element_text(size=18), axis.title=element_text(size=18))+
  theme(axis.text.x = element_text(colour="black")) +ylab ("Risk")

ggsurv_mammaprint = ggarrange(plotmammaprint, tablemammaprint, heights = c(1.5, 0.8),widths = c(0.3,0.3),
                       ncol = 1, nrow = 2, align = "v")



ggsave("data_bitbucket/final_results/plots/MAMMAPRINT.pdf",plotmammaprint, width=6, height=5,dpi = 300)

ggsave("data_bitbucket/final_results/plots/MAMMAPRINT_table.pdf",tablemammaprint, width=6, height=2,dpi = 300)



##########################             ONCOTYPE            ########################## 

a = sub[,c("endpoint","unified_survival_status","unified_survival_days","oncotype_class")]


fit = survfit(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class, data = trial)





ggsurv_oncotype = ggsurvplot(fit,   
                               size = 1.2,
                               data = trial,  
                               censor = FALSE,
                               palette =
                                 c("blue", "red","darkgreen"),
                               risk.table = TRUE,       # show risk table.
                               pval = TRUE,
                               pval.method=TRUE,# show p-value of log-rank test.
                             pval.size = 7,
                             pval.method.size = 7,# show p-value of log-rank test.
                             conf.int = FALSE,         # show confidence intervals for 
                             xscale = "d_y",
                             break.time.by=365.25*2,
                             # point estimates of survival curves.
                             # survival estimates.
                             ylab = "Survival probability",
                             xlab = "Time in years",   # customize X axis label.
                             ggtheme = theme_classic2(base_size=16),
                             risk.table.fontsize = 5,# customize plot and risk table with a theme.
                             risk.table.y.text.col = T, # colour risk table text annotations.
                             risk.table.y.text = FALSE, # show bars instead of names in text annotations
                             risk.table.col = "strata",
                               legend.labs =
                                 c("Low","Intermediate","High")
)


plotoncotype = ggsurv_oncotype$plot +theme(legend.position = c(0.87, 0.25))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))+theme(legend.title = element_text(size=18),
                                                                                 legend.text = element_text(size=18))+guides(color=guide_legend("Risk"))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))



tableoncotype = ggsurv_oncotype$table + theme(legend.position = "none", axis.text=element_text(size=18), axis.title=element_text(size=18))+
  theme(axis.text.x = element_text(colour="black")) +ylab ("Risk")

ggsurv_oncotype = ggarrange(plotoncotype, tableoncotype, heights = c(1.5, 0.8),widths = c(0.3,0.3),
                              ncol = 1, nrow = 2, align = "v")



ggsave("data_bitbucket/final_results/plots/ONCOTYPE_new.pdf",plotoncotype, width=6, height=5,dpi = 300)

ggsave("data_bitbucket/final_results/plots/ONCOTYPE_table_new.pdf",tableoncotype, width=6, height=2,dpi = 300)


##########################            PAM50          ########################## 

a = sub[,c("endpoint","unified_survival_status","unified_survival_days","Call")]

# trial_t = trial[which(trial$Call%!in%"Normal"),]
trial_t = trial
trial_t$Call <- factor(trial_t$Call, levels = c("LumA", "LumB", "Her2", "Basal","Normal"))



fit = survfit(Surv(unified_survival_days, unified_survival_status) ~ Call, data = trial_t)




trial$Call<-as.factor(trial$Call)

trial$Call <- factor(trial$Call, levels = c("LumA","LumB","Her2","Basal","Normal"))


# Calculate the log-rank test to get the p-value
# Calculate Log-rank test
log_rank_test <- survdiff(Surv(unified_survival_days, unified_survival_status) ~ Call, data = trial_t)
p_value <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)

# Format the p-value
formatted_p_value <- if(p_value < 0.001) {
  "p < 0.001"
} else {
  paste("p =", formatC(p_value, format = "f", digits = 3))
}

# Create the Kaplan-Meier plot
ggsurv_pam50 <- ggsurvplot(
  fit, 
  size = 1.2,
  data = trial_t,
  censor = FALSE,
  palette = c("blue", "skyblue2", "hotpink", "red", "green4"),
  risk.table = TRUE,
  pval = FALSE, # Set to FALSE because we're adding it manually
  conf.int = FALSE,
  xscale = "d_y",
  break.time.by = 365.25*2,
  ylab = "Survival probability",
  xlab = "Time in years",
  ggtheme = theme_classic2(base_size = 16),
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 5,
  risk.table.y.text = FALSE,
  risk.table.col = "strata",
  legend.labs = c("LumA", "LumB", "Her2", "Basal", "Normal")
)

# Annotate the plot with "Log-rank" and the formatted p-value
# Adjust the x and y positions as needed for your plot
ggsurv_pam50$plot <- ggsurv_pam50$plot + 
  annotate("text", x = Inf, y = Inf, label = paste("Log-rank", "\n", formatted_p_value), 
           hjust = 1.1, vjust = 4, size = 7.8, color = "black")



plotpam50 = ggsurv_pam50$plot +theme(legend.position = c(0.87, 0.25))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))+theme(legend.title = element_text(size=18),
                                                                                 legend.text = element_text(size=18))+guides(color=guide_legend("Subtype"))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

savepam50 = plotpam50

tablepam50 = ggsurv_pam50$table  + theme(legend.position = "none", axis.text=element_text(size=18), axis.title=element_text(size=18))+
  theme(axis.text.x = element_text(colour="black")) +ylab ("Subtype")

ggsurv_pam50 = ggarrange(plotpam50, tablepam50, heights = c(1.5, 0.8),widths = c(0.3,0.3),
                         ncol = 1, nrow = 2, align = "v")

ggsurv_pam50

ggsave("data_bitbucket/final_results/plots/PAM50.pdf",plotpam50, width=6, height=5,dpi = 300)

ggsave("data_bitbucket/final_results/plots/PAM50_table.pdf",tablepam50, width=6, height=2.1,dpi = 300)

##########################            RORP           ########################## 

a = sub[,c("endpoint","unified_survival_status","unified_survival_days","cell_cycle_risk")]





fit = survfit(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation., data = trial)




ggsurv_rorp = ggsurvplot(fit,   
                             size = 1.2,
                             data = trial,  
                             censor = FALSE,
                             palette =
                               c("blue", "red","darkgreen"),
                             risk.table = TRUE,       # show risk table.
                         pval = TRUE,
                         pval.method=TRUE,# show p-value of log-rank test.
                         pval.size = 7,
                         pval.method.size = 7,
                         conf.int = FALSE,         # show confidence intervals for 
                         xscale = "d_y",
                         break.time.by=365.25*2,
                         # point estimates of survival curves.
                         # survival estimates.
                         ylab = "Survival probability",
                         xlab = "Time in years",   # customize X axis label.
                         ggtheme = theme_classic2(base_size=16), # customize plot and risk table with a theme.
                         risk.table.y.text.col = T, # colour risk table text annotations.
                         risk.table.fontsize = 5,
                         risk.table.y.text = FALSE, # show bars instead of names in text annotations
                         risk.table.col = "strata",
                             legend.labs =
                               c("Low", "Intermediate","High")
)


plotrorp = ggsurv_rorp$plot +theme(legend.position = c(0.87, 0.25))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))+theme(legend.title = element_text(size=18),
                                                                                 legend.text = element_text(size=18))+guides(color=guide_legend("Risk"))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))


tablerorp = ggsurv_rorp$table + theme(legend.position = "none", axis.text=element_text(size=18), axis.title=element_text(size=18))+
  theme(axis.text.x = element_text(colour="black")) +ylab ("Risk")

ggsurv_rorp = ggarrange(plotrorp, tablerorp, heights = c(1.5, 0.8),widths = c(0.3,0.3),
                            ncol = 1, nrow = 2, align = "v")

ggsurv_rorp


ggsave("data_bitbucket/final_results/plots/RORP.pdf",plotrorp, width=6, height=5,dpi = 300)

ggsave("data_bitbucket/final_results/plots/RORP_table.pdf",tablerorp, width=6, height=2,dpi = 300)


##########################            Cell cycle           ########################## 

a = sub[,c("endpoint","unified_survival_status","unified_survival_days","cell_cycle_risk")]




fit = survfit(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk, data = trial)




ggsurv_cellcycle = ggsurvplot(fit,   
                             size = 1.2,
                             data = trial,  
                             censor = FALSE,
                             palette =
                               c("blue", "red","darkgreen"),
                             risk.table = TRUE,       # show risk table.
                             pval = TRUE,
                             pval.method=TRUE,# show p-value of log-rank test.
                             pval.size = 7,
                             pval.method.size = 7,
                             conf.int = FALSE,         # show confidence intervals for 
                             xscale = "d_y",
                             break.time.by=365.25*2,
                             # point estimates of survival curves.
                             # survival estimates.
                             ylab = "Survival probability",
                             xlab = "Time in years",   # customize X axis label.
                             ggtheme = theme_classic2(base_size=16), # customize plot and risk table with a theme.
                             risk.table.y.text.col = T, # colour risk table text annotations.
                             risk.table.fontsize = 5,
                             risk.table.y.text = FALSE, # show bars instead of names in text annotations
                             risk.table.col = "strata",
                             legend.labs =
                               c("Low","Intermediate","High")
)


plotcellcycle = ggsurv_cellcycle$plot +theme(legend.position = c(0.87, 0.25))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18))+theme(legend.title = element_text(size=18),
                                                                                 legend.text = element_text(size=18))+guides(color=guide_legend("Risk"))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))


tablecellcycle = ggsurv_cellcycle$table + theme(legend.position = "none", axis.text=element_text(size=18), axis.title=element_text(size=18))+
  theme(axis.text.x = element_text(colour="black")) +ylab ("Risk")

ggsurv_cellcycle = ggarrange(plotcellcycle, tablecellcycle, heights = c(1.5, 0.8),widths = c(0.3,0.3),
                            ncol = 1, nrow = 2, align = "v")

ggsurv_cellcycle




ggsave("data_bitbucket/final_results/plots/CELL_CYCLE.pdf",plotcellcycle, width=6, height=5,dpi = 300)

ggsave("data_bitbucket/final_results/plots/CELL_CYCLE_table.pdf",tablecellcycle, width=6, height=2,dpi = 300)





##################   COX REGRESSION    ##################################



cox_trial = trial[,c("sample_name","er","pgr","her2","N","T","grade","treatment","tumor_size","unified_survival_days", "unified_survival_status","Call","oncotype_class",
                     "mammaprint.prognosis_binary","ggi_risk","ROR.P.Group..Subtype...Proliferation.","cell_cycle_risk","chemo","hormono","chemo_hormono")]



cox_trial$her2 <- factor(cox_trial$her2, levels = c("negative", "positive"))

cox_trial$pgr <- factor(cox_trial$pgr, levels = c("positive", "negative"))

cox_trial$ggi_risk = as.integer(as.logical(cox_trial$ggi_risk))

cox_trial$ggi_risk <- factor(cox_trial$ggi_risk, levels = c("0", "1"))


cox_trial$mammaprint.prognosis_binary = as.integer(as.logical(cox_trial$mammaprint.prognosis_binary))

cox_trial$mammaprint.prognosis_binary <- factor(cox_trial$mammaprint.prognosis_binary, levels = c("0", "1"))

 cox_trial$Call <- factor(cox_trial$Call, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))



cox_trial$er <- factor(cox_trial$er, levels = c("positive", "negative"))

cox_trial_t = cox_trial

cox_trial_all = cox_trial_t

#######    ALL  PATIENTS #################################



# PAM50
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ Call +  grade + er + N + tumor_size + hormono , data = cox_trial_t))





# RORP
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + er + N + grade + tumor_size + hormono , data = cox_trial_t))




# ONCOTYPE DX
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class + er + grade + N + tumor_size + hormono, data = cox_trial_t))




# MAMMAPRINT
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary + er + grade + N + tumor_size + hormono , data = cox_trial_t))


# GGI
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk + er + grade + N + tumor_size+  hormono, data = cox_trial_t))



# CELL CYCLE
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk + er + grade + N + tumor_size+  hormono , data = cox_trial_t))




#######    ER+/LN+ PATIENTS     #################################


cox_trial = cox_trial_all
cox_trial_erpos_lnpos = cox_trial[which(cox_trial$er=="positive" & cox_trial$N>=1),]


# PAM50
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ Call + grade + tumor_size+ hormono , data = cox_trial_erpos_lnpos))


# ROR-P
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + grade + tumor_size+  hormono, data = cox_trial_erpos_lnpos))




# ONCOTYPE DX
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class +  grade + tumor_size+  hormono , data = cox_trial_erpos_lnpos))




# MAMMAPRINT
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary +  grade + tumor_size+ hormono , data = cox_trial_erpos_lnpos))


# GGI
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk +  grade + tumor_size+ hormono , data = cox_trial_erpos_lnpos))



# CELL CYCLE
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk +  grade + tumor_size+ hormono, data = cox_trial_erpos_lnpos))









#######    ER+/LN-      #################################

cox_trial = cox_trial_all

cox_trial_erpos_lnneg = cox_trial[which(cox_trial$er=="positive" & cox_trial$N==0),]



# PAM50
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ Call + grade + tumor_size+ hormono , data = cox_trial_erpos_lnneg))


# ROR-P
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + grade + tumor_size+  hormono, data = cox_trial_erpos_lnneg))




# ONCOTYPE DX
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class +  grade + tumor_size+  hormono , data = cox_trial_erpos_lnneg))




# MAMMAPRINT
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary +  grade + tumor_size+ hormono , data = cox_trial_erpos_lnneg))


# GGI
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk +  grade + tumor_size+ hormono , data = cox_trial_erpos_lnneg))



# CELL CYCLE
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk +  grade + tumor_size+ hormono, data = cox_trial_erpos_lnneg))





#######    ER+/LN-/HER2-      #################################

cox_trial = cox_trial_all

cox_trial_erpos_lnneg = cox_trial[which(cox_trial$er=="positive" & cox_trial$N==0 & cox_trial$her2%in%"negative"),]



# PAM50
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ Call + grade + tumor_size+ hormono , data = cox_trial_erpos_lnneg))


# ROR-P
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + grade + tumor_size+  hormono, data = cox_trial_erpos_lnneg))




# ONCOTYPE DX
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class +  grade + tumor_size+  hormono , data = cox_trial_erpos_lnneg))




# MAMMAPRINT
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary +  grade + tumor_size+ hormono , data = cox_trial_erpos_lnneg))


# GGI
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk +  grade + tumor_size+ hormono , data = cox_trial_erpos_lnneg))



# CELL CYCLE
summary(coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk +  grade + tumor_size+ hormono, data = cox_trial_erpos_lnneg))





# Compare models:


library("epiDisplay")

################ ALL PATIENTS ######################

cox_trial_t = cox_trial_all

# PAM50
model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size + hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ Call +  er + N +grade + tumor_size + hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance




# GGI
model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono , data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance

# OncotypeDX

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance

# cell cycle

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance


# mammaprint

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance



# rorp

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance






################ ER+/LN- ######################


cox_trial_erpos_lnneg = cox_trial[which(cox_trial$er=="positive" & cox_trial$N==0),]

cox_trial_t = cox_trial_erpos_lnneg

# PAM50
model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size + hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ Call +  er + N +grade + tumor_size + hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance




# GGI
model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono , data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance

# OncotypeDX

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance

# cell cycle

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance


# mammaprint

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance



# rorp

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance





################ ER+/LN+ ######################

cox_trial_erpos_lnpos = cox_trial[which(cox_trial$er=="positive" & cox_trial$N>0),]

cox_trial_t = cox_trial_erpos_lnpos

# PAM50
model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size + hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ Call +  er + N +grade + tumor_size + hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance




# GGI
model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono , data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ ggi_risk +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance

# OncotypeDX

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ oncotype_class +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance

# cell cycle

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance


# mammaprint

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance



# rorp

model0 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ grade + er + N + tumor_size+ hormono, data = cox_trial_t)

model1 <- coxph(Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. +  grade + er + N + tumor_size+ hormono, data = cox_trial_t)
lrtest (model0, model1)
model1$concordance
