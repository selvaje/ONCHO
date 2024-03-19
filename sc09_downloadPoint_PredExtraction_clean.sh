#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc09_downloadPoint_PredExtraction_clean.sh.%j.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc09_downloadPoint_PredExtraction_clean.sh.%j.err
#SBATCH --job-name=sc09_downloadPoint_PredExtraction_clean.sh 
#SBATCH --mem=40G

###### sbatch  /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc09_downloadPoint_PredExtraction_clean.sh 

####  it allows to extract environmental/physical information at point level using gdallocationinfo 

export ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
cd $ONCHO/vector

source ~/bin/gdal3
source ~/bin/pktools

tr '\r' '\n' < $ONCHO/vector/abridged_consolidated_file.csv  | awk '!/^[[:space:]]*$/' > $ONCHO/vector/abridged_consolidated_file_unix.csv

echo "x y pa" > $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.txt

STATUS=adultpresent  ##### change the columns number   larvaepresent=$20 adultpresent=$21 flypresent=$22

awk -F "," '{ if (NR>1) { print $8 , $7 , $21 }}' $ONCHO/vector/abridged_consolidated_file_unix.csv | awk '{ if (NF==3) print }' >> $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.txt
awk  '{ if (NR>1) print $1 , $2 }'       $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.txt   >  $ONCHO/vector/abridged_consolidated_x_y_$STATUS.txt
awk  '{ if (NR>1) print $1 , $2 , $3 }'  $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.txt   >  $ONCHO/vector/abridged_consolidated_x_y_pa_noheader_$STATUS.txt

rm -fr  $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.gpkg
pkascii2ogr   -n "PA"  -a_srs epsg:4326 -f GPKG  -i $ONCHO/vector/abridged_consolidated_x_y_pa_noheader_$STATUS.txt -o $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.gpkg
gdal_rasterize  -co COMPRESS=DEFLATE -co ZLEVEL=9  -ot Byte -a_nodata 255   -a "PA" -tr 0.000833333333333 0.000833333333333  -te 2 4 15 15  -a_srs EPSG:4326 $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.gpkg  $ONCHO/vector/abridged_consolidated_x_y_pa_$STATUS.tif

#### This section allows to create a raster ID layer. It is used to classify the points that fall in the same pixel. When the raster ID layer is created then canbe commented. 

# echo "ncols        15600"                  >  $ONCHO/input/id/raster_ID.asc 
# echo "nrows        13200"                  >> $ONCHO/input/id/raster_ID.asc 
# echo "xllcorner    2"                      >> $ONCHO/input/id/raster_ID.asc 
# echo "yllcorner    4"                      >> $ONCHO/input/id/raster_ID.asc 
# echo "cellsize     0.000833333333333"      >> $ONCHO/input/id/raster_ID.asc 

# awk ' BEGIN {  
# for (row=1 ; row<=13200 ; row++)  { 
#      for (col=1 ; col<=15600 ; col++) { 
#          printf ("%i " ,  col+(row-1)*15600  ) } ; printf ("\n")  }}' >> $ONCHO/input/id/raster_ID.asc 
 
# # transform the created arcinfo ascii grid in a tif.
 
# gdal_translate   -a_srs EPSG:4326 -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/input/id/raster_ID.asc  $ONCHO/input/id/raster_ID.tif 
# rm   $ONCHO/input/id/raster_ID.asc

#### extract raster ID information and calculate the average of point info having the same ID. If the value is 0 then is absence, if the value is > 0 then is presence. 

gdallocationinfo -geoloc -wgs84 -valonly  $ONCHO/input/id/raster_ID.tif  < $ONCHO/vector/abridged_consolidated_x_y_$STATUS.txt  > $ONCHO/vector/abridged_consolidated_x_y_ID_$STATUS.txt 
paste -d " "  $ONCHO/vector/abridged_consolidated_x_y_pa_noheader_$STATUS.txt   $ONCHO/vector/abridged_consolidated_x_y_ID_$STATUS.txt | sort -k 4,4   > $ONCHO/vector/abridged_consolidated_x_y_pa_ID_$STATUS.txt 
bash   /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/average.sh $ONCHO/vector/abridged_consolidated_x_y_pa_ID_$STATUS.txt  $ONCHO/vector/abridged_consolidated_x_y_pa_ID_uniq_$STATUS.txt  <<EOF 
n
4
8
EOF

awk '{if($3==0) {PA=0} else {PA=1}  ;   print $1 , $2 , PA  }'   $ONCHO/vector/abridged_consolidated_x_y_pa_ID_uniq_$STATUS.txt > $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_noheader_$STATUS.txt

rm $ONCHO/vector/abridged_consolidated_x_y_ID_$STATUS.txt $ONCHO/vector/abridged_consolidated_x_y_pa_ID_uniq_$STATUS.txt

awk  '{ print $1 , $2 }'      $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_noheader_$STATUS.txt  >  $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt

echo "x y pa" >  $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_header_$STATUS.txt
awk  '{ print $1 , $2 , $3 }' $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_noheader_$STATUS.txt  >>  $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_header_$STATUS.txt

echo "x y" >  $ONCHO/vector/abridged_consolidated_x_y_uniq_header_$STATUS.txt
awk  '{ print $1 , $2}' $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_noheader_$STATUS.txt  >>  $ONCHO/vector/abridged_consolidated_x_y_uniq_header_$STATUS.txt


#### crate unique point vector file 

pkascii2ogr -n "PA" -a_srs epsg:4326 -f GPKG -i $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_noheader_$STATUS.txt -o $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_$STATUS.gpkg
 
echo  start the predictors  extraction  using gdallocationinfo

rm -f $ONCHO/vector/pred_*.txt 

echo geomorpho90m hydrography90m 

for var in geomorpho90m hydrography90m ; do 

gdalbuildvrt -separate  -overwrite $ONCHO/input/$var/all_tif.vrt $(ls $ONCHO/input/$var/*.tif | grep -v -e _acc -e _msk)
BB=$(ls $ONCHO/input/$var/*.tif | grep -v -e _acc -e _msk  | wc -l ) # count the numbers of tif that are stack in the vrt. This is usefull for transfer from vertical to orizontal the gdallocationinfo information

for var1 in $( ls $ONCHO/input/$var/*.tif  | grep -v -e _acc -e _msk )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}.txt
echo "" >> $ONCHO/vector/pred_${var}.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/$var/all_tif.vrt < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}.txt

done   

echo  chelsa soilgrids soiltemp 

for var in chelsa soilgrids soiltemp ; do 

gdalbuildvrt -separate  -overwrite $ONCHO/input/$var/all_tif.vrt $(ls $ONCHO/input/$var/*_r.tif | grep -v -e _acc -e _msk)
BB=$(ls $ONCHO/input/$var/*_r.tif | grep -v -e _acc -e _msk  | wc -l )

for var1 in $( ls $ONCHO/input/$var/*_r.tif  | grep -v -e _acc -e _msk )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}.txt
echo "" >> $ONCHO/vector/pred_${var}.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/$var/all_tif.vrt < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt | awk -v BB=$BB  'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}.txt

done #### close the for var in chelsa soilgrids soiltemp 

var=soilgrids
gdalbuildvrt -separate  -overwrite $ONCHO/input/$var/all_tif_acc.vrt $(ls $ONCHO/input/$var/*.tif | grep -e _acc | grep -v _r.tif  )
BB=$(ls $ONCHO/input/$var/*.tif | grep -e _acc | grep -v _r.tif | wc -l )

for var1 in $( ls $ONCHO/input/$var/*.tif  | grep -e _acc | grep -v _r.tif )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}_acc.txt
echo "" >> $ONCHO/vector/pred_${var}_acc.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/$var/all_tif_acc.vrt < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}_acc.txt

echo vegetation index 

echo "ndti" > $ONCHO/vector/pred_ndti.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndti.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_ndti.txt
echo "ndvi" > $ONCHO/vector/pred_ndvi.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndvi.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_ndvi.txt
echo "ndwi" > $ONCHO/vector/pred_ndwi.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndwi.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_ndwi.txt

echo livestock 

gdalbuildvrt -separate  -overwrite $ONCHO/input/livestock/all_tif.vrt  $ONCHO/input/livestock/LS_??.tif 
BB=$(ls  $ONCHO/input/livestock/LS_??.tif  | wc -l )

for var in  $ONCHO/input/livestock/LS_??.tif; do echo -n $(basename $var .tif)" " ; done > $ONCHO/vector/pred_LS.txt
echo "" >> $ONCHO/vector/pred_LS.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/livestock/all_tif.vrt < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_LS.txt


echo  population 

echo "GHSpop_90m" > $ONCHO/vector/pred_pop.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/population/GHSpop_90m.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_pop.txt

echo water 
echo "occurence" > $ONCHO/vector/pred_occ.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/water/occurence.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_occ.txt

echo "occurence_proximity" > $ONCHO/vector/pred_occ_prox.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/water/occurence_proximity.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_occ_prox.txt

echo "flow_mean" > $ONCHO/vector/pred_flow.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/water/flow_mean.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_flow.txt 

echo   landcover 

echo "LC2021" > $ONCHO/vector/pred_landcover.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/landcover/LC2021.tif  < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_landcover.txt

echo ecorigion 

# removed due to
# Error: Missing data in columns: ER2017.
# echo "ER2017" > $ONCHO/vector/pred_ecoreg.txt  
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/ecoregions/ER2017.tif   < $ONCHO/vector/abridged_consolidated_x_y_uniq_noheader_$STATUS.txt  >> $ONCHO/vector/pred_ecoreg.txt


echo merging predictors 

paste -d " " $ONCHO/vector/abridged_consolidated_x_y_pa_uniq_header_$STATUS.txt $ONCHO/vector/pred_*.txt |  sed  's,  , ,g' > $ONCHO/vector/x_y_pa_predictors_$STATUS.txt
awk '{ print $1="", $2="", $0  }' $ONCHO/vector/x_y_pa_predictors_$STATUS.txt |  sed  's,    ,,g'  > $ONCHO/vector/x_y_pa_predictors4R_$STATUS.txt


