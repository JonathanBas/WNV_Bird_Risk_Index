# Code to reproduce:
# Mapping the Bird Risk Index for West Nile virus in Europe and its relationship with disease occurrence in humans
# J. Bastard, R. Metras, B. Durand

############################################################################
## PART 3 - Computes maps of BRI's mean, coefficient of variation and P50 ##
############################################################################

rm(list=ls())
library(dplyr)
library(terra)
library(raster)
library(readxl)
q_val <- 0.5

## NUTS shapefiles

nuts3a <- vect("./data/NUTS_RG_10M_2016_3035_LEVL_3.shp")
nuts3b <- nuts3a[!nuts3a$CNTR_CODE %in% c("AL","CH","IS","LI","ME","MK","NO","RS","TR","UK") & !nuts3a$NUTS_ID %in% c(c("FRY10","FRY20","FRY30","FRY40","FRY50","PT200","PT300","ES703","ES704","ES705","ES706","ES707","ES708","ES709"))]

## List of 1,000 BRI rasters

setwd("./data/bootmaps/")
list_rast <- list.files(pattern = ".tif", recursive=T)
n_rast <- length(list_rast)

rast_1_full <- terra::rast(list_rast[1])
nuts3c <- project(x=nuts3b, y=rast_1_full)

## Computes rasters of BRI mean, coefficient of variation and P50

mean_val_full = mean_valsquared_full = propor_val_sup_crop = 0

for(i in 1:n_rast){
  cat(paste0(i, "/", n_rast, "\r"))
  rast_i_full <- terra::rast(list_rast[i])
  
  if(same.crs(nuts3c, rast_i_full)){ # Should always be TRUE but just a check
    rast_i_crop <- terra::crop(x=rast_i_full, y=nuts3c, mask=T, overwrite=T)
    
    val_i_full <- terra::values(rast_i_full)
    valsquared_i_full <- val_i_full ^2
    val_i_crop <- terra::values(rast_i_crop)

    mean_val_full <- mean_val_full + (val_i_full /n_rast)
    mean_valsquared_full <- mean_valsquared_full + (valsquared_i_full /n_rast)
    
    quant_i_crop <- stats::quantile(val_i_crop, probs=q_val, na.rm=T)
    val_sup_i_crop <- as.numeric(val_i_crop >= quant_i_crop)
    propor_val_sup_crop <- propor_val_sup_crop + (val_sup_i_crop /n_rast)
  }
}

SD_val_full <- sqrt(mean_valsquared_full - (mean_val_full ^2))
CV_val_full <- SD_val_full /mean_val_full

rast_mean = rast_CV = rast_i_full
values(rast_mean) <- mean_val_full
values(rast_CV) <- CV_val_full

rast_prop_pix = rast_i_crop
values(rast_prop_pix) <- propor_val_sup_crop

writeRaster(rast_mean, paste0("./mean_",n_rast,"_bootstraps.tif"), overwrite=T)
writeRaster(rast_CV, paste0("./CV_",n_rast,"_bootstraps.tif"), overwrite=T)
writeRaster(rast_prop_pix, paste0("./propor_",q_val,"_pix_",n_rast,"_bootstraps.tif"), overwrite=T)

## Computes lower and upper bounds of the BRI 95% confidence interval from a subset

if(n_rast >=200){
  
  list_samp_rast <- c(list.files(path="./boot_1_100", pattern=".tif", recursive=T, full.names=T),
                      list.files(path="./boot_101_200", pattern=".tif", recursive=T, full.names=T))
  
  rast_samp <- rast(list_samp_rast)
  rast_CI_low <- terra::app(rast_samp, fun = function(x){quantile(x, probs=0.025, na.rm=T)})
  rast_CI_high <- terra::app(rast_samp, fun = function(x){quantile(x, probs=0.975, na.rm=T)})
  
  writeRaster(rast_CI_low, paste0("./CI_low_SampleOf200_bootstraps.tif"), overwrite=T)
  writeRaster(rast_CI_high, paste0("./CI_high_SampleOf200_bootstraps.tif"), overwrite=T)
}


