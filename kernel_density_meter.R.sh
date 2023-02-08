#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/kernel_density_meter.R.sh.%j.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/kernel_density_meter.R.sh.%j.err
#SBATCH --job-name=kernel_density_meter.R.sh
#SBATCH --mem=40G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/kernel_density_meter.R.sh

#Goal is to make kernel density map from coordinates, with dimensions below
# class       : SpatRaster 
# dimensions  : 13200, 15600, 1  (nrow, ncol, nlyr)
# resolution  : 0.0008333333, 0.0008333333  (x, y)
# extent      : 2, 15, 4, 15  (xmin, xmax, ymin, ymax)
# coord. ref. : lon/lat WGS 84 (EPSG:4326) 
# source      : elevation.tif 
# name        : elevation 

Rscript --vanilla  -e '

setwd("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/pa_kernel")

#https://www.samuelbosch.com/2014/02/creating-kernel-density-estimate-map-in.html

library(KernSmooth)
library(raster)
library(rgdal)

records <- read.table("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header.txt" , sep=" " ,  header = TRUE)
coordinates_wgs84 <- subset(records , pa ==  1 )[,1:2]
coordinates(coordinates_wgs84) <- ~x+y
proj4string(coordinates_wgs84)=CRS("+init=epsg:4326")

# create a raster projected in meter and get the extend in the Equidistant projection
# nigeria https://epsg.io/26392
# +proj=tmerc +lat_0=4 +lon_0=8.5 +k=0.99975 +x_0=670553.98 +y_0=0 +a=6378249.145 +rf=293.465 +towgs84=-92,-93,122,0,0,0,0 +units=m +no_defs +type=crs

coordinates_meter <- spTransform(coordinates_wgs84,CRS("+proj=tmerc +lat_0=4 +lon_0=8.5 +k=0.99975 +x_0=670553.98 +y_0=0 +a=6378249.145 +rf=293.465 +towgs84=-92,-93,122,0,0,0,0 +units=m +no_defs +type=crs"))

reaster_test = raster(nrows=13200, ncols=15600, xmn=2, xmx=15, ymn=4, ymx=15, vals=1, crs="+init=epsg:4326" )

projected_raster <- projectRaster(reaster_test, crs = "+proj=tmerc +lat_0=4 +lon_0=8.5 +k=0.99975 +x_0=670553.98 +y_0=0 +a=6378249.145 +rf=293.465 +towgs84=-92,-93,122,0,0,0,0 +units=m +no_defs +type=crs" , res=100)
ext= extent(projected_raster)

projected_raster

# compute the 2D binned kernel density estimate

est <- bkde2D(coordinates_meter@coords,
              bandwidth=c(200000,200000),  # 200 000 is in meter  so 200 km  
              gridsize=c(12276, 14474),
              range.x=list(c(ext[1] ,ext[2] ),c(ext[3],ext[4])))

# create raster
est.raster = raster(list(x=est$x1,y=est$x2, z=est$fhat))
projection(est.raster) <- "+proj=tmerc +lat_0=4 +lon_0=8.5 +k=0.99975 +x_0=670553.98 +y_0=0 +a=6378249.145 +rf=293.465 +towgs84=-92,-93,122,0,0,0,0 +units=m +no_defs +type=crs"

est.raster_wgs84 <- projectRaster(est.raster, crs = "+init=epsg:4326"  , res=0.008333333333)

est.raster_wgs84

est.raster_wgs84_crop = crop ( est.raster_wgs84 , extent ( 2, 15, 4, 15  )  )

est.raster_wgs84_crop

# the  est.raster_wgs84_crop as extent slithly different that the other raster , I would export to tif and adjust with gdal_edit.py 

raster::writeRaster (est.raster_wgs84_crop, "kernel.tif", gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="FLT4S" , NAflag=-9999)

save.image(paste0("kernel.RData"))

'

