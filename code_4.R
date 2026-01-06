# Code to reproduce:
# Mapping the Bird Risk Index for West Nile virus in Europe and its relationship with disease occurrence in humans
# J. Bastard, R. Metras, B. Durand

########################################################
## PART 4 - Aggregation of BRI maps at the NUTS level ##
########################################################

rm(list=ls())
library(dplyr)
library(terra)
library(ggplot2)
library(tidyterra)
library(ggplot2)
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

## For what countries should we use NUTS2 instead of NUTS3 ?

sum_ar <- as.data.frame(arnuts3 %>%
                          merge(y=arnuts2, by="code_nuts2") %>%
                          group_by(code_cntr) %>%
                          summarise(n_nuts2 = length(unique(area_nuts2)),
                                    mean_area_nuts2 = mean(unique(area_nuts2), na.rm=T),
                                    n_nuts3 = length(area_nuts3),
                                    mean_area_nuts3 = mean(unique(area_nuts3), na.rm=T)))
mean_nuts3_ar_cntr <- mean(sum_ar$mean_area_nuts3)

diff_mean_ar <- data.frame(country = sum_ar$code_cntr,
                           diff_nuts2 = abs(sum_ar$mean_area_nuts2 - mean_nuts3_ar_cntr),
                           diff_nuts3 = abs(sum_ar$mean_area_nuts3 - mean_nuts3_ar_cntr))

sum_ar$opt_nuts_lvl <- 1+apply(X=diff_mean_ar[,c(2,3)], MARGIN=1, FUN=which.min)
print(sum_ar[,c("code_cntr", "opt_nuts_lvl")])

## Therefore we use NUTS2 instead of NUTS3 for BE, DE, MT and NL (for CY and LU, NUTS2=NUTS3)

## NUTS shapefiles

nuts1 <- vect("./data/CNTR_RG_10M_2016_3035.shp")
nuts3a <- vect("./data/NUTS_RG_10M_2016_3035_LEVL_3.shp")
nuts3b <- nuts3a[!nuts3a$CNTR_CODE %in% c("AL","CH","IS","LI","ME","MK","NO","RS","TR","UK") & !nuts3a$NUTS_ID %in% c(c("FRY10","FRY20","FRY30","FRY40","FRY50","PT200","PT300","ES703","ES704","ES705","ES706","ES707","ES708","ES709"))]

nuts3b_bis <- terra::merge(x=nuts3b, y=cor_nuts_23, by.x="FID", by.y="code_nuts3", all.x=T)
nuts3b_bis$nuts_bis_id <- nuts3b_bis$FID
nuts3b_bis$nuts_bis_id[nuts3b_bis$CNTR_CODE %in% c("BE","DE","MT","NL")] <- nuts3b_bis$code_nuts2[nuts3b_bis$CNTR_CODE %in% c("BE","DE","MT","NL")]
nuts3b_bis <- terra::aggregate(nuts3b_bis, by="nuts_bis_id", fun=base::min)
names(nuts3b_bis)[2] <- "LEVL_CODE"
nuts3b_bis$LEVL_CODE[nuts3b_bis$CNTR_CODE %in% c("BE","DE","MT","NL")] <- 2
nuts3b_bis$FID = nuts3b_bis$NUTS_ID = nuts3b_bis$nuts_bis_id
plot(nuts3b_bis)
same.crs(nuts3b, nuts3b_bis) # must be TRUE

## Aggregates BRI maps at the NUTS level
## PIOU_moyen = mean Bird Risk Index
## propor_pix_50 = P50

aggreg_nuts3 <- function(tif_file_input, shp_file_output){
  rast <- terra::rast(paste0("./data/", tif_file_input))

  if(! same.crs(rast, nuts3b_bis)){
    rast <- project(x=rast, y=nuts3b_bis)
  }
  same.crs(rast, nuts3b_bis)
  
  aggr_pix_val2 <- data.frame(FID=as.character(NULL), val=as.numeric(NULL))
  for(nut_i in 1:length(nuts3b_bis$FID)){
    cat(paste0(nut_i, "/", length(nuts3b_bis$FID), "\r"))
    
    reg_i <- nuts3b_bis[nuts3b_bis$FID == nuts3b_bis$FID[nut_i]]
    rast_i <- terra::extract(x=rast, y=reg_i)
    mean_i <- mean(rast_i[,2], na.rm=T)
    
    aggr_pix_val2 <- rbind(aggr_pix_val2,
                           data.frame(FID=nuts3b_bis$FID[nut_i],
                                      val=mean_i))
  }
  
  nuts3c <- terra::merge(x=nuts3b_bis, y=aggr_pix_val2, by="FID")
  nuts3c <- terra::subset(x=nuts3c, subset=rep(T,nrow(nuts3c)), select=c("FID","LEVL_CODE","NUTS_ID","CNTR_CODE","NUTS_NAME","agg_n","val"))
  terra::plot(nuts3c, "val", type="continuous")
  writeVector(nuts3c, filename=paste0("./data/", shp_file_output, "_par_nuts3.shp"), overwrite=TRUE)
}

aggreg_nuts3("propor_0.5_pix_1000_bootstraps.tif", "propor_pix_50")
aggreg_nuts3("mean_1000_bootstraps.tif", "PIOU_moyen")
aggreg_nuts3("richness.tif", "spec_richness")
