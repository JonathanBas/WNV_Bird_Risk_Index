# Code to reproduce:
# Mapping the Bird Risk Index for West Nile virus in Europe and its relationship with disease occurrence in humans
# J. Bastard, R. Metras, B. Durand

###################################################################
## PART 5 - Relationship between the BRI and WNV cases (Model 2) ##
###################################################################

rm(list=ls())
library(dplyr)
library(terra)
library(regions)
library(parallel)
library(readxl)
library(ggplot2)
library(ggspatial)
library(ggpubr)

## Data on area of NUTS2 and NUTS3 regions in Europe

arnuts3 <- as.data.frame(read.table(file="./data/reg_area_NUTS3.csv", sep=",", header=T) %>%
                           dplyr::select(geo, OBS_VALUE) %>%
                           filter(grepl("BE|BG|CZ|DK|DE|EE|IE|EL|ES|FR|HR|IT|CY|LV|LT|LU|HU|MT|NL|AT|PL|PT|RO|SI|SK|FI|SE", geo)) %>%
                           rename(area_nuts3=OBS_VALUE, code_nuts3=geo) %>%
                           mutate(code_nuts2 = substr(code_nuts3, 1, nchar(code_nuts3)-1),
                                  code_cntr = substr(code_nuts3, 1, 2)))

arnuts2 <- as.data.frame(read.table(file="./data/reg_area_NUTS2.csv", sep=",", header=T) %>%
                           dplyr::select(geo, OBS_VALUE) %>%
                           filter(grepl("BE|BG|CZ|DK|DE|EE|IE|EL|ES|FR|HR|IT|CY|LV|LT|LU|HU|MT|NL|AT|PL|PT|RO|SI|SK|FI|SE", geo)) %>%
                           rename(area_nuts2=OBS_VALUE, code_nuts2=geo))

cor_nuts_23 <- merge(x=arnuts3, y=arnuts2, by="code_nuts2")[,c("code_cntr","code_nuts2","code_nuts3")]

## Data: BRI and associated variables
PIOU_moyen <- vect("./data/PIOU_moyen_par_nuts3.shp") ; names(PIOU_moyen)[7] <- "PIOU_moyen"
spec_rich <- vect("./data/spec_richness_par_nuts3.shp") ; names(spec_rich)[7] <- "spec_rich"
propor_pix_50 <- vect("./data/propor_pix_50_par_nuts3.shp") ; names(propor_pix_50)[7] <- "propor_pix_50"
var_map <- PIOU_moyen %>%
  terra::merge(y = spec_rich, by = c("FID","LEVL_CODE","NUTS_ID","CNTR_CODE","NUTS_NAME","agg_n")) %>%
  terra::merge(y = propor_pix_50, by = c("FID","LEVL_CODE","NUTS_ID","CNTR_CODE","NUTS_NAME","agg_n"))

## Population data
popnuts <- read.table(file="./data/demo_r_pjangrp3_linear.csv", sep=",", header=T) %>%
  filter(sex=="T", age=="TOTAL", TIME_PERIOD==2016) %>%
  rename(popul = OBS_VALUE)

## Aggregates cases data at the NUTS2 level for BE,DE,MT,NL
cases <- as.data.frame(read.table(file = "./data/clean_db_cas_human.csv", sep=";", header=T) %>%
                         filter(! is.na(code_nuts3)) %>%
                         merge(y=cor_nuts_23, by="code_nuts3", all.x=T) %>%
                         mutate(code_utilis = code_nuts3))
cases$code_utilis[cases$code_cntr %in% c("BE","DE","MT","NL")] <- cases$code_nuts2[cases$code_cntr %in% c("BE","DE","MT","NL")]
cases <- as.data.frame(cases %>%
                         group_by(Country,code_cntr,Year,code_utilis) %>%
                         summarise(Date_first_case = min(Date_first_case),
                                   Total_cases = sum(Total_cases)))

## Aggregates cases in all years
cases_aggr <- as.data.frame(cases %>%
                              group_by(code_utilis) %>%
                              summarise(any_case = any(Total_cases >0),
                                        nb_cases = sum(Total_cases),
                                        nb_years_with_cases = sum(Total_cases >0)))

## Adds case data to the raster
nuts_cases <- terra::merge(x=var_map, y=cases_aggr, by.x="FID", by.y="code_utilis", all.x=T)
nuts_cases <- terra::merge(x=nuts_cases, y=popnuts[,c("geo","popul")], by.x="FID", by.y="geo", all.x=T)
nuts_cases$any_case[is.na(nuts_cases$any_case)] <- 0
nuts_cases$nb_cases[is.na(nuts_cases$nb_cases)] <- 0
nuts_cases$nb_years_with_cases[is.na(nuts_cases$nb_years_with_cases)] <- 0

## Country shapefiles
nuts1 <- vect("./data/CNTR_RG_10M_2016_3035.shp")
nuts1 <- nuts1[nuts1$CNTR_ID %in% c("BE","BG","CZ","DK","DE","EE","IE","EL","ES","FR","HR","IT","CY","LV","LT","LU","HU","MT","NL","AT","PL","PT","RO","SI","SK","FI","SE","LI","NO","CH","BA","ME","MD","MK","AL","RS","TR","UA","XK","BY","UK","RU","MA","DZ","TN","LY","EG","IL","PS","LB","SY","JO","IS","AD","GE","GL","SJ","FO")]
nuts1_sf <- sf::st_as_sf(nuts1)

## Associates some countries together

nuts_cases$CNTR_CODE[nuts_cases$CNTR_CODE == "CY"] <- "EL"
nuts_cases$CNTR_CODE[nuts_cases$CNTR_CODE == "LU"] <- "BE"
nuts_cases$CNTR_CODE[nuts_cases$CNTR_CODE == "MT"] <- "IT"

## Prepares models

library(INLA) # INLA version 23.04.24
library(spdep)
library(parallel)
library(spatialreg)
library(ape)
library(stringr)
library(gstat)
library(bizicount)

nuts_cases_sf <- sf::st_as_sf(nuts_cases)

## Computes adjacency matrices
adj_mat <- poly2nb(nuts_cases_sf)
rowstand_adj <- nb2mat(adj_mat, style = "W", zero.policy = T) ; rowstand_adj[rowstand_adj == 0] = 0.001

## Functions to compute posterior predictions and PIT histograms
options(mc.cores = 7)
source("./def_fun_predict_PIT.R")

## Run models

run_modBYM <- function(outcom_var, mod_predic, likel, incl_bym_rand_eff=T){
  mod_formula <- as.formula(paste(outcom_var, paste(mod_predic, collapse=" + "), sep=" ~ "))
  if(incl_bym_rand_eff){
    mod_formula <- update(mod_formula, . ~. +
                            f(NUTS_ID, model = "bym2", graph = rowstand_adj,
                              values = as.factor(unique(nuts_cases_sf$NUTS_ID))) +
                            f(CNTR_CODE, model = "iid") )
  }
  mod.bym <- inla(mod_formula,
                  data = as.data.frame(nuts_cases_sf),
                  family = likel, Ntrials=8,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE, config = TRUE),
                  control.predictor = list(compute = TRUE))
  
  save(mod.bym, file = paste("./res/mod", outcom_var, likel, paste(mod_predic,collapse="-"), "spat", incl_bym_rand_eff, ".Rdata"))
}

## Model 2 - Number of years with cases
run_modBYM("nb_years_with_cases", c("PIOU_moyen","log(popul)"), "binomial") # Main analysis
run_modBYM("nb_years_with_cases", c("1"), "binomial", incl_bym_rand_eff=T) # Null model WITH spatial effect
run_modBYM("nb_years_with_cases", c("1"), "binomial", incl_bym_rand_eff=F) # Null model WITHOUT spatial effect

run_modBYM("nb_years_with_cases", c("log(popul)"), "binomial") # Sensitivity analysis
run_modBYM("nb_years_with_cases", c("spec_rich","log(popul)"), "binomial") # Sensitivity analysis
run_modBYM("nb_years_with_cases", c("propor_pix_50","log(popul)"), "binomial") # Sensitivity analysis

## Model 2bis - Cumulated number of cases
run_modBYM("nb_cases", c("PIOU_moyen","offset(log(popul))"), "nbinomial") # Main analysis
run_modBYM("nb_cases", c("1"), "nbinomial", incl_bym_rand_eff=T) # Null model WITH spatial effect
run_modBYM("nb_cases", c("1"), "nbinomial", incl_bym_rand_eff=F) # Null model WITHOUT spatial effect

## Model predictions and deviance residuals

lklhd_2A <- "binomial"
load(file = paste0("./res/mod nb_years_with_cases ",lklhd_2A," PIOU_moyen-log(popul) spat TRUE .Rdata")) ; mod_2A<-mod.bym ; rm(mod.bym)
nuts_cases_sf$predic_numb_years <- predic_mod(mod_2A)
nuts_cases_sf$resid_deviance_num_yr1 <- sign(nuts_cases_sf$nb_years_with_cases - nuts_cases_sf$predic_numb_years) * sqrt(2* (log(choose(8,nuts_cases_sf$nb_years_with_cases) * (nuts_cases_sf$nb_years_with_cases/8)^nuts_cases_sf$nb_years_with_cases * (1-(nuts_cases_sf$nb_years_with_cases/8))^(8-nuts_cases_sf$nb_years_with_cases)) - log(choose(8,nuts_cases_sf$nb_years_with_cases) * (nuts_cases_sf$predic_numb_years/8)^nuts_cases_sf$nb_years_with_cases * (1-(nuts_cases_sf$predic_numb_years/8))^(8-nuts_cases_sf$nb_years_with_cases)))) # Definition by [Feng et al 2020, BMC Medical Research Methodology]

lklhd_2B <- "nbinomial"
load(file = paste0("./res/mod nb_cases ",lklhd_2B," PIOU_moyen-offset(log(popul)) spat TRUE .Rdata")) ; mod_2B<-mod.bym ; rm(mod.bym)
nuts_cases_sf$predic_numb_cas <- predic_mod(mod_2B)
param_k <- mod_2B$summary.hyperpar["size for the nbinomial observations (1/overdispersion)", "mean"]
nuts_cases_sf$resid_deviance_num_case <- sign(nuts_cases_sf$nb_cases - nuts_cases_sf$predic_numb_cas) * sqrt(2 * (log((nuts_cases_sf$nb_cases/nuts_cases_sf$predic_numb_cas) ^nuts_cases_sf$nb_cases) - log(((nuts_cases_sf$nb_cases + param_k) /(nuts_cases_sf$predic_numb_cas + param_k)) ^(nuts_cases_sf$nb_cases + param_k))))

## Summarize results

display_mod_BYM <- function(outcom_var, mod_predic, likel, incl_bym_rand_eff=T, plot_fit=F, common_scale=T){
  load(file = paste("./res/mod", outcom_var, likel, paste(mod_predic,collapse="-"), "spat", incl_bym_rand_eff, ".Rdata"))
  mod_formula <- paste(outcom_var, paste(mod_predic, collapse=" + "), sep=" ~ ")
  if(incl_bym_rand_eff){
    mod_formula <- paste0(mod_formula, " + random effects")
  }
  
  tcoef <- mod.bym$summary.fixed#[-1,]
  if(nrow(tcoef)>0){
    par_names <- rownames(tcoef)
    par_val <- paste0(round(tcoef$`0.5quant`,2)," [",round(tcoef$`0.025quant`,2),"; ",round(tcoef$`0.975quant`,2),"]")
  }else{
    par_names = par_val = "-"
  }
  resul_mod <- data.frame(model = c(mod_formula, rep("", max(0, nrow(tcoef)-1))),
                          likelihood = c(likel, rep("", max(0, nrow(tcoef)-1))),
                          DIC = c(round(mod.bym$dic$dic, 1), rep("", max(0, nrow(tcoef)-1))),
                          parameter = par_names,
                          estimate = par_val)
  
  plotmap <- NULL
  if(plot_fit){
    nuts_cases_sf$BYM <- tmarg(mod.bym$marginals.fitted.values)
    title_fig <- ifelse(outcom_var == "nb_years_with_cases",
                        "number\nof years with\ncases",
                        ifelse(outcom_var == "nb_cases",
                               "number\nof cases",
                               ifelse(outcom_var == "any_case", "occurence\nof cases", NA)))
    plotmap <- map_plot(c(outcom_var,"BYM"),
                        c(paste0("Observed ", title_fig, ":"), paste0("Predicted ", title_fig, " (BYM):")),
                        sf_db=nuts_cases_sf, comm_scale=common_scale)
  }
  return(resul_mod)
}

an_sens_mod_2A <- rbind(display_mod_BYM("nb_years_with_cases", c("PIOU_moyen","log(popul)"), "binomial"),
                        display_mod_BYM("nb_years_with_cases", c("1"), "binomial", incl_bym_rand_eff=T),
                        display_mod_BYM("nb_years_with_cases", c("1"), "binomial", incl_bym_rand_eff=F),
                        display_mod_BYM("nb_years_with_cases", c("log(popul)"), "binomial"),
                        display_mod_BYM("nb_years_with_cases", c("spec_rich","log(popul)"), "binomial"),
                        display_mod_BYM("nb_years_with_cases", c("propor_pix_50","log(popul)"), "binomial"))
write.table(an_sens_mod_2A, file = "./res/an_sens_mod_2.csv", sep=";", row.names=F)

## PIT histograms

library(tscount)
graphics.off()
png(filename=paste0("./res/PIT_histo_models_fulldata_2A",lklhd_2A,"_2B",lklhd_2B,".png"),pointsize=6,res=300,width = 8, height = 4.2, units = "cm")
par(mfrow=c(1,2), mgp=c(1.5,0.3,0), mar=c(3,2.8,1.2,0.5))
comp_pit(vec_obs = nuts_cases_sf$nb_years_with_cases, modres = mod_2A, tit_fig="Model 2")
comp_pit(vec_obs = nuts_cases_sf$nb_cases, modres = mod_2B, tit_fig="Model 2bis", y_lab=NA)
dev.off()
graphics.off()

## Cross-validation

library(tscount)
alread_comput <- T
fixed_predic <- "PIOU_moyen"
form_2A <- as.formula(paste0("nb_years_with_cases ~ ",fixed_predic," + log(popul)"))
likel_2A <- "binomial"
form_2B <- as.formula(paste0("nb_cases ~ ",fixed_predic," + offset(log(popul))"))
likel_2B <- "nbinomial"
spat_eff <- "bym2"

nuts_cases_df2 = as.data.frame(nuts_cases_sf)

nuts_cases_df2$categ <- paste(nuts_cases_df2$CNTR_CODE,
                              nuts_cases_df2$any_case,
                              sep="_")
if(alread_comput){
  nuts_cases_df2$grp <- read.table(file="./res/NUTS_group_affiliation.csv", sep=";", header=T)$grp
}else{
  nuts_cases_df2$grp <- NA
  all_NUTS_same_grp <- T
  while(any(all_NUTS_same_grp)){
    all_NUTS_same_grp <- c()
    for(cat_i in unique(nuts_cases_df2$categ)){
      n_in_categ <- sum(nuts_cases_df2$categ == cat_i)
      which_group <- sample(x=1:5, size=n_in_categ, replace=T)
      nuts_cases_df2$grp[nuts_cases_df2$categ == cat_i] <- which_group
      all_NUTS_same_grp <- c(all_NUTS_same_grp, ((length(unique(which_group)) == 1) & n_in_categ!=1) | (length(unique(which_group))<5 & n_in_categ>=8))
    }
  }
  write.table(nuts_cases_df2[,c("NUTS_ID","CNTR_CODE","grp")], file="./res/NUTS_group_affiliation.csv", row.names=F, sep=";")
}

## Cross-validation of Model 2:

resCV_2A <- data.frame()
graphics.off()
png(filename=paste0("./res/PIT_histo_model_2A_",likel_2A,".png"),pointsize=6,res=300,width = 13, height = 6, units = "cm")
layout(mat = matrix(1:10, nrow=2, ncol=5, byrow=F))
par(mgp=c(1.5,0.3,0), mar=c(3,2.8,1.2,0.5))
for (fold_k in 1:5){
  print(paste("Running model 2A - Fold", fold_k))
  
  ### Observations in the testing dataset
  which_test <- (nuts_cases_df2$grp == fold_k)
  
  ### In the training dataset, outcomes is NA
  train_dat <- nuts_cases_df2
  train_dat[which_test, c("nb_years_with_cases", "nb_cases")] <- NA
  
  ## Fit the models on the training data
  if(! alread_comput){
    modCV_null_2A <- inla(nb_years_with_cases ~ 1 +
                            f(NUTS_ID, model = spat_eff, graph = rowstand_adj, values = as.factor(unique(nuts_cases_sf$NUTS_ID))) +
                            f(CNTR_CODE, model = "iid"),
                          data = train_dat,
                          family = likel_2A, Ntrials=8,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE, config = TRUE),
                          control.predictor = list(compute = TRUE, link = 1))
    save(modCV_null_2A, file = paste0("./res/modCV_2A ",likel_2A," nb_years_with_cases null_withSpatial ",fold_k,".Rdata"))
    
    modCV_2A <- inla(update(form_2A, . ~. +
                              f(NUTS_ID, model = spat_eff, graph = rowstand_adj, values = as.factor(unique(nuts_cases_sf$NUTS_ID))) +
                              f(CNTR_CODE, model = "iid")),
                     data = train_dat,
                     family = likel_2A, Ntrials=8,
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE, config = TRUE),
                     control.predictor = list(compute = TRUE, link = 1))
    save(modCV_2A, file = paste0("./res/modCV_2A ",likel_2A," nb_years_with_cases PIOU_moyen-log(popul) ",fold_k,".Rdata"))
    
  }else{
    load(file = paste0("./res/modCV_2A ",likel_2A," nb_years_with_cases PIOU_moyen-log(popul) ",fold_k,".Rdata"))
    load(file = paste0("./res/modCV_2A ",likel_2A," nb_years_with_cases null_withSpatial ",fold_k,".Rdata")) ; modCV_null_withSpatial_2A <- modCV_null_2A ; rm(modCV_null_2A)
    
    modCV_null_noSpatial_2A <- inla(as.formula(paste(as.character(form_2A)[2], 1, sep=" ~ ")),
                                    data = train_dat,
                                    family = likel_2A, Ntrials=8,
                                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE, config = TRUE),
                                    control.predictor = list(compute = TRUE, link = 1))
    
    ### Posterior coefficient estimates
    estim_interc <- paste0(round(modCV_2A$summary.fixed["(Intercept)","0.5quant"],2),"\n[",round(modCV_2A$summary.fixed["(Intercept)","0.025quant"],2),"; ",round(modCV_2A$summary.fixed["(Intercept)","0.975quant"],2),"]")
    estim_PIOU_moy <- paste0(round(modCV_2A$summary.fixed["PIOU_moyen","0.5quant"],2),"\n[",round(modCV_2A$summary.fixed["PIOU_moyen","0.025quant"],2),"; ",round(modCV_2A$summary.fixed["PIOU_moyen","0.975quant"],2),"]")
    estim_log_popul <- paste0(round(modCV_2A$summary.fixed["log(popul)","0.5quant"],2),"\n[",round(modCV_2A$summary.fixed["log(popul)","0.025quant"],2),"; ",round(modCV_2A$summary.fixed["log(popul)","0.975quant"],2),"]")
    
    precis_BYM2 <- paste0(round(modCV_2A$summary.hyperpar["Precision for NUTS_ID","0.5quant"],2),"\n[",round(modCV_2A$summary.hyperpar["Precision for NUTS_ID","0.025quant"],2),"; ",round(modCV_2A$summary.hyperpar["Precision for NUTS_ID","0.975quant"],2),"]")
    phi_BYM2 <- paste0(round(modCV_2A$summary.hyperpar["Phi for NUTS_ID","0.5quant"],2),"\n[",round(modCV_2A$summary.hyperpar["Phi for NUTS_ID","0.025quant"],2),"; ",round(modCV_2A$summary.hyperpar["Phi for NUTS_ID","0.975quant"],2),"]")
    precis_countr_eff <- paste0(round(modCV_2A$summary.hyperpar["Precision for CNTR_CODE","0.5quant"],3),"\n[",round(modCV_2A$summary.hyperpar["Precision for CNTR_CODE","0.025quant"],3),"; ",round(modCV_2A$summary.hyperpar["Precision for CNTR_CODE","0.975quant"],3),"]")
    
    ### Check GoF
    predic_2A <- predic_mod(modCV_2A)
    
    val_R2_train_noSpatial <- round(R2_mod(modCV_2A, modCV_null_noSpatial_2A, vec_subset = !which_test, dat_vect=nuts_cases_df2$nb_years_with_cases), 2)
    val_R2_test_noSpatial <- round(R2_mod(modCV_2A, modCV_null_noSpatial_2A, vec_subset = which_test, dat_vect=nuts_cases_df2$nb_years_with_cases), 2)
    
    val_R2_train_withSpatial <- round(R2_mod(modCV_2A, modCV_null_withSpatial_2A, vec_subset = !which_test, dat_vect=nuts_cases_df2$nb_years_with_cases), 2)
    val_R2_test_withSpatial <- round(R2_mod(modCV_2A, modCV_null_withSpatial_2A, vec_subset = which_test, dat_vect=nuts_cases_df2$nb_years_with_cases), 2)
    
    deltaDIC_nullmod_noSpatial <- round(modCV_2A$dic$dic - modCV_null_noSpatial_2A$dic$dic, 1)
    deltaDIC_nullmod_withSpatial <- round(modCV_2A$dic$dic - modCV_null_withSpatial_2A$dic$dic, 1)
    
    resCV_2A <- rbind(resCV_2A,
                      data.frame(Numb = fold_k,
                                 Intercept = estim_interc,
                                 Mean_RI = estim_PIOU_moy,
                                 log_popul = estim_log_popul,
                                 precis_BYM2 = precis_BYM2,
                                 phi_BYM2 = phi_BYM2,
                                 precis_countr_eff = precis_countr_eff,
                                 R2_test_noSpatial = val_R2_test_noSpatial,
                                 R2_test_withSpatial = val_R2_test_withSpatial,
                                 deltaDIC_nullmod_noSpatial = deltaDIC_nullmod_noSpatial,
                                 deltaDIC_nullmod_withSpatial = deltaDIC_nullmod_withSpatial,
                                 R2_train_noSpatial = val_R2_train_noSpatial,
                                 R2_train_withSpatial = val_R2_train_withSpatial,
                                 log_score_train = round(-mean(log(modCV_2A$cpo$cpo[! which_test])), 2),
                                 nullmod_precis_BYM2 = round(modCV_null_withSpatial_2A$summary.hyperpar["Precision for NUTS_ID","mean"], 3),
                                 nullmod_phi_BYM2 = round(modCV_null_withSpatial_2A$summary.hyperpar["Phi for NUTS_ID","mean"], 3),
                                 nullmod_precis_countr_eff = round(modCV_null_withSpatial_2A$summary.hyperpar["Precision for CNTR_CODE","mean"], 3)))
    
    p_gof = ggplot() +
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      geom_point(aes(x=nuts_cases_df2$nb_years_with_cases, y=predic_2A, col=which_test)) +
      xlab("Observed") + ylab("Predicted")
    
    comp_pit(vec_obs = nuts_cases_df2$nb_years_with_cases, modres = modCV_2A, which_id = which(! which_test), tit_fig=paste0("Fold ",fold_k," - Train"), x_lab=NA, y_lab=ifelse(fold_k==1,"Density",NA))
    comp_pit(vec_obs = nuts_cases_df2$nb_years_with_cases, modres = modCV_2A, which_id = which(which_test), tit_fig=paste0("Fold ",fold_k," - Test"), y_lab=ifelse(fold_k==1,"Density",NA))
    
  }
}
dev.off()
write.table(resCV_2A, file=paste0("./res/summary_cross-validation_model_2A_",likel_2A,".csv"), sep=";", row.names=F)
graphics.off()

## Cross-validation of Model 2bis:

resCV_2B <- data.frame()
graphics.off()
png(filename=paste0("./res/PIT_histo_model_2B_",likel_2B,".png"),pointsize=6,res=300,width = 13, height = 6, units = "cm")
layout(mat = matrix(1:10, nrow=2, ncol=5, byrow=F))
par(mgp=c(1.5,0.3,0), mar=c(3,2.8,1.2,0.5))
for (fold_k in 1:5){
  print(paste("Running model 2B - Fold", fold_k))
  
  ### Observations in the testing dataset
  which_test <- (nuts_cases_df2$grp == fold_k)
  
  ### In the training dataset, outcomes is NA
  train_dat <- nuts_cases_df2
  train_dat[which_test, c("nb_years_with_cases", "nb_cases")] <- NA
  
  test_dat <- nuts_cases_df2[which_test,]
  
  ## Fit the models on the training data
  if(! alread_comput){
    
    modCV_null_2B <- inla(nb_cases ~ 1 + f(NUTS_ID, model = spat_eff, graph = rowstand_adj, values = as.factor(unique(nuts_cases_sf$NUTS_ID))) + f(CNTR_CODE, model = "iid"),
                          data = train_dat,
                          family = likel_2B,
                          control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE, config = TRUE),
                          control.predictor = list(compute = TRUE, link = 1))
    save(modCV_null_2B, file = paste0("./res/modCV_2B ",likel_2B," nb_cases null_withSpatial ",fold_k,".Rdata"))
    
    modCV_2B <- inla(update(form_2B, . ~. + f(NUTS_ID, model = spat_eff, graph = rowstand_adj, values = as.factor(unique(nuts_cases_sf$NUTS_ID))) + f(CNTR_CODE, model = "iid") ),
                     data = train_dat,
                     family = likel_2B,
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE),
                     control.predictor = list(compute = TRUE, link = 1))
    save(modCV_2B, file = paste0("./res/modCV_2B ",likel_2B," nb_cases PIOU_moyen-offset(log(popul)) ",fold_k,".Rdata"))
    
  }else{
    load(file = paste0("./res/modCV_2B ",likel_2B," nb_cases PIOU_moyen-offset(log(popul)) ",fold_k,".Rdata"))
    load(file = paste0("./res/modCV_2B ",likel_2B," nb_cases null_withSpatial ",fold_k,".Rdata")) ; modCV_null_withSpatial_2B <- modCV_null_2B ; rm(modCV_null_2B)
    
    modCV_null_noSpatial_2B <- inla(as.formula(paste(as.character(form_2B)[2], 1, sep=" ~ ")),
                                    data = train_dat,
                                    family = likel_2B,
                                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE),
                                    control.predictor = list(compute = TRUE, link = 1))
    
    ### Posterior coefficient estimates
    estim_interc <- paste0(round(modCV_2B$summary.fixed["(Intercept)","0.5quant"],2)," [",round(modCV_2B$summary.fixed["(Intercept)","0.025quant"],2),"; ",round(modCV_2B$summary.fixed["(Intercept)","0.975quant"],2),"]")
    estim_PIOU_moy <- paste0(round(modCV_2B$summary.fixed["PIOU_moyen","0.5quant"],2)," [",round(modCV_2B$summary.fixed["PIOU_moyen","0.025quant"],2),"; ",round(modCV_2B$summary.fixed["PIOU_moyen","0.975quant"],2),"]")
    
    precis_BYM2 <- paste0(round(modCV_2B$summary.hyperpar["Precision for NUTS_ID","0.5quant"],2),"\n[",round(modCV_2B$summary.hyperpar["Precision for NUTS_ID","0.025quant"],2),"; ",round(modCV_2B$summary.hyperpar["Precision for NUTS_ID","0.975quant"],2),"]")
    phi_BYM2 <- paste0(round(modCV_2B$summary.hyperpar["Phi for NUTS_ID","0.5quant"],2),"\n[",round(modCV_2B$summary.hyperpar["Phi for NUTS_ID","0.025quant"],2),"; ",round(modCV_2B$summary.hyperpar["Phi for NUTS_ID","0.975quant"],2),"]")
    precis_countr_eff <- paste0(round(modCV_2B$summary.hyperpar["Precision for CNTR_CODE","0.5quant"],3),"\n[",round(modCV_2B$summary.hyperpar["Precision for CNTR_CODE","0.025quant"],3),"; ",round(modCV_2B$summary.hyperpar["Precision for CNTR_CODE","0.975quant"],3),"]")
    
    comp_pit(vec_obs = nuts_cases_df2$nb_cases, modres = modCV_2B, which_id = which(! which_test), tit_fig=paste0("Fold ",fold_k," - Train"), x_lab=NA, y_lab=ifelse(fold_k==1,"Density",NA))
    comp_pit(vec_obs = nuts_cases_df2$nb_cases, modres = modCV_2B, which_id = which(which_test), tit_fig=paste0("Fold ",fold_k," - Test"), y_lab=ifelse(fold_k==1,"Density",NA))
  }
}
dev.off()
write.table(resCV_2B, file=paste0("./res/summary_cross-validation_model_2B_",likel_2B,".csv"), sep=";", row.names=F)
graphics.off()


