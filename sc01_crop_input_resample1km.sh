#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 10:00:00
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc01_crop_input_resample1km.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc01_crop_input_resample1km.sh.%J.err
#SBATCH --job-name=sc01_crop_input.sh
#SBATCH --mem=40G

module load StdEnv
source ~/bin/gdal3
source ~/bin/pktools
             
##### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc01_crop_input_resample1km.sh


ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO

echo geomorpho90m 
for var in  $ONCHO/input/geomorpho90m/*.tif  ;  do 
gdalwarp  -r average   -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333  -srcnodata -9999  -dstnodata -9999 -overwrite $var  $ONCHO/input/geomorpho90m_1km/$(basename  $var )
done

echo hydrography90m
for var in   $ONCHO/input/hydrography90m/*.tif   ;  do
gdalwarp  -r max    -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333  -srcnodata -2147483648   -dstnodata -2147483648  -overwrite $var  $ONCHO/input/hydrography90m_1km/$(basename  $var )
done

echo  chelsa

for var in  $ONCHO/input/chelsa/CHELSA_bio?.tif $ONCHO/input/chelsa/CHELSA_bio??.tif  ; do 
cp  $var  $ONCHO/input/chelsa_1km/$(basename  $var .tif)_r.tif 
done 

echo   soiltemp

for var in  $( ls $ONCHO/input/soiltemp/*.tif | grep -v -e _acc -e _msk  -e _r.tif  )  ; do
cp  $var  $ONCHO/input/soiltemp_1km/$(basename  $var .tif)_r.tif	
done

echo soilgrids 

for var in $(ls $ONCHO/input/soilgrids/*.tif  | grep -v -e _acc -e _msk  -e _r.tif )  ; do
gdalwarp  -r average   -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333  -srcnodata -9999999   -dstnodata -9999999  -overwrite $var  $ONCHO/input/soilgrids_1km/$(basename  $var .tif)_r.tif 
done 

for var in $(ls $ONCHO/input/soilgrids/*.tif  | grep  _acc  )  ; do
gdalwarp -r max -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333  -srcnodata -9999999   -dstnodata -9999999  -overwrite $var  $ONCHO/input/soilgrids_1km/$(basename  $var )
done

echo vegetation index 

gdalwarp -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -srcnodata -19  -dstnodata -19 -overwrite $ONCHO/input/sentinel/ndti.tif $ONCHO/input/sentinel_1km/ndti.tif  

gdalwarp -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -srcnodata -19  -dstnodata -19 -overwrite $ONCHO/input/sentinel/ndvi.tif $ONCHO/input/sentinel_1km/ndvi.tif

gdalwarp -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -srcnodata -19  -dstnodata -19 -overwrite $ONCHO/input/sentinel/ndwi.tif $ONCHO/input/sentinel_1km/ndwi.tif

echo livestock

for var in  $ONCHO/input/livestock/LS_??.tif  ; do
cp $var $ONCHO/input/livestock_1km/$(basename  $var )
done

echo  population 

echo "GHSpop_90m"  
gdalwarp -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -srcnodata -200  -dstnodata -200 -overwrite $ONCHO/input/population/GHSpop_90m.tif $ONCHO/input/population_1km/GHSpop_90m.tif  

echo water 
echo "occurence" 
gdalwarp -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -overwrite $ONCHO/input/water/occurence.tif $ONCHO/input/water_1km/occurence.tif

echo "occurence_proximity" 
gdalwarp -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -overwrite $ONCHO/input/water/occurence_proximity.tif $ONCHO/input/water_1km/occurence_proximity.tif

echo "flow_mean"
cp $ONCHO/input/water/flow_mean.tif $ONCHO/input/water_1km/flow_mean.tif

echo   landcover 

echo "LC2021" 
gdalwarp -r mode -co COMPRESS=DEFLATE -co ZLEVEL=9 -tr 0.00833333333333 0.00833333333333 -overwrite $ONCHO/input/landcover/LC2021.tif  $ONCHO/input/landcover_1km/LC2021.tif 




exit 

        
