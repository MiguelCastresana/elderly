## =========================
## Libraries
## =========================
library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(ggpubr)

## =========================
## Inputs
## =========================
sub <- results_final              
out_dir <- "data_bitbucket/final_results/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## Utilities
## =========================
recode_low_med_high <- function(x) {
  x <- dplyr::recode(x,
                     "low"="Low","Low"="Low",
                     "med"="Intermediate","Medium"="Intermediate","Intermediate"="Intermediate",
                     "high"="High","High"="High",
                     .default = as.character(x)
  )
  factor(x, levels = c("Low","Intermediate","High"))
}

## =========================
## Clean / harmonize data
## =========================
trial <- sub %>%
  mutate(
    RS = as.integer(!!RS),
    pgr  = dplyr::recode(pgr,  `+`="positive", `-`="negative", .default = pgr),
    her2 = dplyr::recode(her2, `+`="positive", `-`="negative", .default = her2),
    unified_survival_days   = as.numeric(unified_survival_days),
    unified_survival_status = as.integer(as.factor(unified_survival_status)),
  ) %>%
  # Censor at 10y
  mutate(
    unified_survival_days   = pmin(unified_survival_days, 10 * 365.25),
    unified_survival_status = ifelse(unified_survival_days > 10 * 365.25, 1L, unified_survival_status)
  ) %>%
  # Signature harmonization
  mutate(
    ROR.P.Group..Subtype...Proliferation. = recode_low_med_high(ROR.P.Group..Subtype...Proliferation.),
    oncotype_class       = recode_low_med_high(oncotype_class),
    cell_cycle_risk      = recode_low_med_high(cell_cycle_risk),
    grade                = as.factor(grade),
    treatment            = as.factor(treatment),
    Call = factor(Call, levels = c("LumA","LumB","Her2","Basal","Normal"))
  )

## =========================
## KM plot helper
## =========================
km_plot <- function(data, time, status, group, palette, legend_labs, legend_title,
                    file_stub, show_pval = TRUE, manual_pval = FALSE) {
  
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time, status, group))
  fit <- survfit(f, data = data)
  
  pval_arg <- FALSE
  pval_method_arg <- FALSE
  ann_layer <- NULL
  
  if (manual_pval) {
    lr <- survdiff(f, data = data)
    p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
    ptxt <- if (p < 0.001) "p < 0.001" else paste0("p = ", formatC(p, digits = 3, format = "f"))
    ann_layer <- annotate("text", x = Inf, y = Inf, label = paste("Log-rank\n", ptxt),
                          hjust = 1.1, vjust = 4, size = 7.8, color = "black")
  } else if (show_pval) {
    pval_arg <- TRUE
    pval_method_arg <- TRUE
  }
  
  g <- ggsurvplot(
    fit,
    data = data,
    size = 1.2,
    censor = FALSE,
    palette = palette,
    risk.table = TRUE,
    pval = pval_arg,
    pval.method = pval_method_arg,
    pval.size = 7,
    pval.method.size = 7,
    conf.int = FALSE,
    xscale = "d_y",
    break.time.by = 365.25 * 2,
    ylab = "Survival probability",
    xlab = "Time in years",
    ggtheme = theme_classic2(base_size = 16),
    risk.table.y.text.col = TRUE,
    risk.table.fontsize = 5,
    risk.table.y.text = FALSE,
    risk.table.col = "strata",
    legend.labs = legend_labs
  )
  
  plot_top <- g$plot +
    theme(legend.position = c(0.87, 0.25),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) +
    guides(color = guide_legend(legend_title))
  
  if (!is.null(ann_layer)) plot_top <- plot_top + ann_layer
  
  plot_tab <- g$table +
    theme(legend.position = "none",
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(colour = "black")) +
    ylab(legend_title)
  
  arranged <- ggarrange(plot_top, plot_tab, heights = c(1.5, 0.8), widths = c(0.3, 0.3),
                        ncol = 1, nrow = 2, align = "v")
  
  ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), plot_top, width = 6, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(file_stub, "_table.pdf")), plot_tab, width = 6, height = 2, dpi = 300)
  
  invisible(arranged)
}

## =========================
## KM plots
## =========================
# GGI (binary; assumes 0/1 or FALSE/TRUE in trial$ggi_risk)
trial <- trial %>% mutate(ggi_risk = factor(as.integer(as.logical(ggi_risk)), levels = c(0,1), labels = c("GG1","GG3")))
km_plot(trial, "unified_survival_days", "unified_survival_status", "ggi_risk",
        palette = c("blue","red"),
        legend_labs = c("GG1","GG3"),
        legend_title = "Grade",
        file_stub = "GGI")

# MammaPrint (binary)
trial <- trial %>% mutate(mammaprint.prognosis_binary = factor(as.integer(as.logical(mammaprint.prognosis_binary)),
                                                               levels = c(0,1), labels = c("Low","High")))
km_plot(trial, "unified_survival_days", "unified_survival_status", "mammaprint.prognosis_binary",
        palette = c("blue","red"),
        legend_labs = c("Low","High"),
        legend_title = "Risk",
        file_stub = "MAMMAPRINT")

# Oncotype (Low/Intermediate/High)
km_plot(trial, "unified_survival_days", "unified_survival_status", "oncotype_class",
        palette = c("blue","red","darkgreen"),
        legend_labs = c("Low","Intermediate","High"),
        legend_title = "Risk",
        file_stub = "ONCOTYPE_new")

# PAM50 (with manual p-value annotation)
km_plot(trial, "unified_survival_days", "unified_survival_status", "Call",
        palette = c("blue","skyblue2","hotpink","red","green4"),
        legend_labs = c("LumA","LumB","Her2","Basal","Normal"),
        legend_title = "Subtype",
        file_stub = "PAM50",
        show_pval = FALSE, manual_pval = TRUE)

# ROR-P (Low/Intermediate/High)
km_plot(trial, "unified_survival_days", "unified_survival_status", "ROR.P.Group..Subtype...Proliferation.",
        palette = c("blue","red","darkgreen"),
        legend_labs = c("Low","Intermediate","High"),
        legend_title = "Risk",
        file_stub = "RORP")

# Cell cycle (Low/Intermediate/High)
km_plot(trial, "unified_survival_days", "unified_survival_status", "cell_cycle_risk",
        palette = c("blue","red","darkgreen"),
        legend_labs = c("Low","Intermediate","High"),
        legend_title = "Risk",
        file_stub = "CELL_CYCLE")

## =========================
## Cox setup
## =========================
cox_trial <- trial %>%
  transmute(
    sample_name, er = factor(er, levels = c("positive","negative")),
    pgr = factor(pgr, levels = c("positive","negative")),
    her2 = factor(her2, levels = c("negative","positive")),
    N, T, grade, treatment, tumor_size,
    unified_survival_days, unified_survival_status,
    Call = factor(Call, levels = c("LumA","LumB","Her2","Basal","Normal")),
    oncotype_class,
    mammaprint.prognosis_binary = factor(as.integer(as.logical(mammaprint.prognosis_binary)), levels = c(0,1)),
    ggi_risk = factor(as.integer(as.logical(ggi_risk)), levels = c(0,1)),
    ROR.P.Group..Subtype...Proliferation.,
    cell_cycle_risk,
    chemo, hormono, chemo_hormono
  )

cox_trial_all <- cox_trial

## =========================
## Cox helpers
## =========================
fit_and_print <- function(formula, data) {
  fm <- as.formula(formula)
  s <- summary(coxph(fm, data = data))
  print(s)
  invisible(s)
}

compare_add_marker <- function(base_terms, marker, data) {
  base_f <- as.formula(paste0("Surv(unified_survival_days, unified_survival_status) ~ ", base_terms))
  add_f  <- as.formula(paste0("Surv(unified_survival_days, unified_survival_status) ~ ", marker, " + ", base_terms))
  m0 <- coxph(base_f, data = data)
  m1 <- coxph(add_f,  data = data)
  print(epiDisplay::lrtest(m0, m1))
  print(m1$concordance)
  invisible(list(m0=m0, m1=m1))
}

## =========================
## Cox: All patients
## =========================
dat <- cox_trial_all
base <- "grade + er + N + tumor_size + hormono"

fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ Call +  grade + er + N + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + er + N + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ oncotype_class + er + grade + N + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary + er + grade + N + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ggi_risk + er + grade + N + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk + er + grade + N + tumor_size + hormono", dat)

compare_add_marker(base, "Call", dat)
compare_add_marker(base, "ggi_risk", dat)
compare_add_marker(base, "oncotype_class", dat)
compare_add_marker(base, "cell_cycle_risk", dat)
compare_add_marker(base, "mammaprint.prognosis_binary", dat)
compare_add_marker(base, "ROR.P.Group..Subtype...Proliferation.", dat)

## =========================
## Cox: ER+/LN-
## =========================
dat <- cox_trial_all %>% filter(er == "positive", N == 0)
base <- "grade + er + N + tumor_size + hormono"

fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ Call + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ oncotype_class + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ggi_risk + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk + grade + tumor_size + hormono", dat)

compare_add_marker(base, "Call", dat)
compare_add_marker(base, "ggi_risk", dat)
compare_add_marker(base, "oncotype_class", dat)
compare_add_marker(base, "cell_cycle_risk", dat)
compare_add_marker(base, "mammaprint.prognosis_binary", dat)
compare_add_marker(base, "ROR.P.Group..Subtype...Proliferation.", dat)

## =========================
## Cox: ER+/LN+
## =========================
dat <- cox_trial_all %>% filter(er == "positive", N > 0)
base <- "grade + er + N + tumor_size + hormono"

fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ Call + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ oncotype_class + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ggi_risk + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk + grade + tumor_size + hormono", dat)

compare_add_marker(base, "Call", dat)
compare_add_marker(base, "ggi_risk", dat)
compare_add_marker(base, "oncotype_class", dat)
compare_add_marker(base, "cell_cycle_risk", dat)
compare_add_marker(base, "mammaprint.prognosis_binary", dat)
compare_add_marker(base, "ROR.P.Group..Subtype...Proliferation.", dat)

## =========================
## Cox: ER+/LN-/HER2-
## =========================
dat <- cox_trial_all %>% filter(er == "positive", N == 0, her2 == "negative")

fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ Call + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ROR.P.Group..Subtype...Proliferation. + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ oncotype_class + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ mammaprint.prognosis_binary + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ ggi_risk + grade + tumor_size + hormono", dat)
fit_and_print("Surv(unified_survival_days, unified_survival_status) ~ cell_cycle_risk + grade + tumor_size + hormono", dat)
