#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc16_prob_masking.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc16_prob_masking.sh.%J.err
#SBATCH --job-name=sc16_prob_masking.sh
#SBATCH --mem=10G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc16_prob_masking.sh 

export ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
source ~/bin/gdal3
source ~/bin/pktools

GDAL_CACHEMAX=30000

cd $ONCHO/

# pksetmask  -co  COMPRESS=LZW  -co ZLEVEL=9  -m $ONCHO/input/openstreetmap/gis_osm_roads_free_1_buf.tif  -msknodata 0 -nodata -9999 -i    $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif -o  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf.tif

# gdal_translate -of XYZ $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf.tif $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf.xyz

# grep -v  "\-9999" $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf.xyz > $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf_clean.xyz

# sort -k 3,3 -g   $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf_clean.xyz  >  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf_clean_s.xyz


pkgetmask -ot Byte    -co  COMPRESS=LZW  -co ZLEVEL=9 -min 0.6 -max 1 -data 1 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o $ONCHO/prediction_all/prediction_seedP_all_msk_mean_60.tif 

pkgetmask -ot Byte    -co  COMPRESS=LZW  -co ZLEVEL=9 -min 0.55 -max 0.6 -data 1 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o $ONCHO/prediction_all/prediction_seedP_all_msk_mean_55.tif 

pkgetmask -ot Byte    -co  COMPRESS=LZW  -co ZLEVEL=9 -min 0.50 -max 0.55 -data 1 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o $ONCHO/prediction_all/prediction_seedP_all_msk_mean_50.tif 

rm -fr $ONCHO/prediction_all/rand_point??.gpkg
pkextractogr -rand 100000  -r point  -srcnodata 0 -i  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_60.tif    -f  "GPKG" -o $ONCHO/prediction_all/rand_point60.gpkg

pkextractogr -rand 10000   -r point  -srcnodata 0 -i  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_55.tif    -f  "GPKG" -o $ONCHO/prediction_all/rand_point55.gpkg

pkextractogr -rand 10000   -r point  -srcnodata 0 -i  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_50.tif    -f  "GPKG" -o $ONCHO/prediction_all/rand_point50.gpkg

exit 
### on my michine 
rm -f sqlite rand_pointP??.gpkg
ogr2ogr -sql "SELECT ST_PointOnSurface(geom), * FROM points" -dialect sqlite rand_pointP50.gpkg  rand_point50.gpkg
ogr2ogr -sql "SELECT ST_PointOnSurface(geom), * FROM points" -dialect sqlite rand_pointP55.gpkg  rand_point55.gpkg
ogr2ogr -sql "SELECT ST_PointOnSurface(geom), * FROM points" -dialect sqlite rand_pointP60.gpkg  rand_point60.gpkg

rm -f sqlite rand_pointP??.txt 
ogrinfo -al rand_pointP50.gpkg | grep POINT  | awk '{ gsub("\\(" , ""); gsub("\\)" , ""); print $2 , $3  }' > rand_pointP50.txt
ogrinfo -al rand_pointP55.gpkg | grep POINT  | awk '{ gsub("\\(" , ""); gsub("\\)" , ""); print $2 , $3  }' > rand_pointP55.txt
ogrinfo -al rand_pointP60.gpkg | grep POINT  | awk '{ gsub("\\(" , ""); gsub("\\)" , ""); print $2 , $3  }' > rand_pointP60.txt

rm -f rand_pointP??_prob.*
pkextractogr -f "GPKG"  -i prediction_seedP_all_msk_mean.tif -s rand_pointP50.gpkg  -o rand_pointP50_prob.gpkg
pkextractogr -f "GPKG"  -i prediction_seedP_all_msk_mean.tif -s rand_pointP55.gpkg  -o rand_pointP55_prob.gpkg
pkextractogr -f "GPKG"  -i prediction_seedP_all_msk_mean.tif -s rand_pointP60.gpkg  -o rand_pointP60_prob.gpkg

rm -f rand_pointP??_prob.csv
paste -d "," <(awk '{ print $1 "," $2 }' rand_pointP50.txt )  <(awk '{ print $1 "," $2 }' rand_pointP50.txt ) <(gdallocationinfo -geoloc -valonly prediction_seedP_all_msk_mean.tif < rand_pointP50.txt) > rand_pointP50_prob.csv
paste -d "," <(awk '{ print $1 "," $2 }' rand_pointP55.txt )  <(awk '{ print $1 "," $2 }' rand_pointP55.txt ) <(gdallocationinfo -geoloc -valonly prediction_seedP_all_msk_mean.tif < rand_pointP55.txt) > rand_pointP55_prob.csv
paste -d "," <(awk '{ print $1 "," $2 }' rand_pointP60.txt )  <(awk '{ print $1 "," $2 }' rand_pointP60.txt ) <(gdallocationinfo -geoloc -valonly prediction_seedP_all_msk_mean.tif < rand_pointP60.txt) > rand_pointP60_prob.csv

rm -f rand_pointP??_prob.gpkg
pkascii2ogr -f "GPKG" -fs "," -n "Lon"  -n "Lat" -n "prob_mean"  -i   rand_pointP60_prob.csv  -o rand_pointP60_prob.gpkg 
pkascii2ogr -f "GPKG" -fs "," -n "Lon"  -n "Lat" -n "prob_mean"  -i   rand_pointP55_prob.csv  -o rand_pointP55_prob.gpkg 
pkascii2ogr -f "GPKG" -fs "," -n "Lon"  -n "Lat" -n "prob_mean"  -i   rand_pointP50_prob.csv  -o rand_pointP50_prob.gpkg 



paste -d "," <( awk '{if (NR>1) print $1 "," $2 }' NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.txt)  <(awk '{ if (NR>1)  print $1 "," $2 "," $3 "," $4 "," $5  }' NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.txt) > NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.csv 

rm NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.gpkg
pkascii2ogr -f "GPKG" -fs ","  -n "Lon"  -n "Lat" -n "PA" -n "prob_mean" -n "prob_stdev"  -i  NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.csv -o  NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.gpkg    



exit 



rm -fr /tmp/mylocation 
apptainer exec  ~/bin/grass8.sif bash <<EOF
grass -f --text -c  -c $ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf60.tif  -e /tmp/mylocation
grass /tmp/mylocation/PERMANENT --exec r.external  input=$ONCHO/prediction_all/prediction_seedP_all_msk_mean_mskbuf60.tif output=mask  --overwrite
grass /tmp/mylocation/PERMANENT --exec g.region raster=mask  
grass /tmp/mylocation/PERMANENT --exec r.random -z -s input=mask  npoints=10 vector=rand_point  --o
grass /tmp/mylocation/PERMANENT --exec v.out.ogr format=GPKG input=rand_point output=$ONCHO/prediction_all/rand_point.gpkg --o
EOF
