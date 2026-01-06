# Code to reproduce:
# Mapping the Bird Risk Index for West Nile virus in Europe and its relationship with disease occurrence in humans
# J. Bastard, R. Metras, B. Durand

#############################################################################
## PART 1 - Prediction of WNV seroprevalence in 150 bird species (Model 1) ##
#############################################################################

rm(list=ls())
library(tidyr)
library(dplyr)
library(dvmisc)
library(pracma)
library(ggplot2)
library(lme4)
library(parameters)
library(ggpubr)

## Prepares data

lif_histo <- read.table("./data/Supplementary_Data_1.csv", sep=",", header=T, check.names=F) # Data on bird species' life history traits
sero_data1 <- read.table("./data/Supplementary_Data_2.csv", sep=",", header=T, check.names=F) # Data from WNV seroprevalence studies

sero_data2 <- as.data.frame(sero_data1
                            %>% uncount(Sampled)
                            %>% group_by(ref, latin)
                            %>% summarise(spec_pres = c(base::rep(1, unique(Pos)), base::rep(0, length(Pos) - unique(Pos))))
                            %>% merge(y=lif_histo, by.x="latin", by.y="Species", all.x=T))

sero_data2$`Migratory status` <- relevel(as.factor(sero_data2$`Migratory status`), ref = "No")
sero_data2$`Nest height` <- relevel(as.factor(sero_data2$`Nest height`), ref = "Intermediate")
sero_data2$`Use of urban/suburban habitats` <- relevel(as.factor(sero_data2$`Use of urban/suburban habitats`), ref = "No")
sero_data2$`Nocturnal gregariousness` <- relevel(as.factor(sero_data2$`Nocturnal gregariousness`), ref = "No")
sero_data2$`Exposure of nestlings` <- relevel(as.factor(sero_data2$`Exposure of nestlings`), ref = "Precocial")
sero_data2$`Breeding sociality` <- relevel(as.factor(sero_data2$`Breeding sociality`), ref = "No")

## Runs Model 1

modpred <- glmer(spec_pres ~ `Migratory status` + mass_scaled + `Nest height` + `Use of urban/suburban habitats` + `Nocturnal gregariousness` + `Exposure of nestlings` + `Breeding sociality` + (1|ref),
                 data = sero_data2, family="binomial", control=glmerControl(optimizer ='bobyqa'))

tab_coef = data.frame(param = rownames(summary(modpred)$coefficients),
                      estim_CI = paste0(round(exp(modpred@beta),2)," [",round(exp(ci(modpred)$CI_low),2),"; ",round(exp(ci(modpred)$CI_high),2),"]"),
                      p_val = round(summary(modpred)$coefficients[,4],3),
                      p_val_raw = summary(modpred)$coefficients[,4])
print(tab_coef)

pred = sigmoid(predict(modpred, newdata = lif_histo, re.form=NA))
names(pred) = lif_histo$Species
write.table(data.frame(spec=names(pred), pred=pred), file="./res/pred_full_data.csv", row.names=F, sep=";")

## Runs bootstrap for WNV seroprevalence predictions (Model 1)

n_repeat = 100
pred_boots = matrix(NA, length(lif_histo$Species), n_repeat)
for (rep_i in 1:n_repeat){
  cat("\r", paste0(rep_i, "/", n_repeat))
  modpred <- glmer(spec_pres ~ `Migratory status` + mass_scaled + `Nest height` + `Use of urban/suburban habitats` + `Nocturnal gregariousness` + `Exposure of nestlings` + `Breeding sociality` + (1|ref),
                   data = sero_data2[sample(1:nrow(sero_data2), size=round(1 * nrow(sero_data2)), replace=T),],
                   family="binomial", control=glmerControl(optimizer ='bobyqa'))
  
  pred_boots[,rep_i] = sigmoid(predict(modpred, newdata = lif_histo, re.form=NA))
}
rownames(pred_boots) = lif_histo$Species

write.table(pred_boots, file="./res/pred_boots_bobyqa.csv", row.names=T, sep=";")

## Makes Supplementary Figure S3

plotsero <- data.frame(spec = lif_histo$Species,
                      pred_fulldata = sigmoid(predict(modpred, newdata = lif_histo, re.form=NA)),
                      ci_lwr = apply(X=as.matrix(pred_boots), MARGIN=1, FUN=function(i){quantile(as.numeric(i), probs=0.025)}),
                      ci_upr = apply(X=as.matrix(pred_boots), MARGIN=1, FUN=function(i){quantile(as.numeric(i), probs=0.975)}))
plotsero <- plotsero[order(plotsero$spec),]
plotsero$colum <- ifelse((1:nrow(plotsero)) < nrow(plotsero)/2, "A", "B")

p_pred_serop <- function(colu){
  p <- ggplot(data = filter(plotsero, colum==colu), aes(y=spec)) +
    geom_pointrange(aes(x=pred_fulldata, xmin=ci_lwr, xmax=ci_upr)) +
    xlab("Seroprevalence") +
    xlim(0, max(plotsero$ci_upr)) +
    theme_bw() +
    theme(axis.text.y = element_text(vjust=0.3, hjust=1),
          legend.box.background = element_rect(color="black"))
  if(colu=="A"){
    p = p + ylab("Species") + theme(legend.position=c(0.76, 0.94))
  }else{
    p = p + ylab(NULL) + theme(legend.position="none")
  }
  return(p)
}

p_predic_sero <- ggarrange(plotlist = list(p_pred_serop("A"), p_pred_serop("B")), ncol=2)

png(filename="./res/Figure_S3.png",pointsize=6,res=300,width = 30, height = 22, units = "cm")
p_predic_sero
dev.off()


