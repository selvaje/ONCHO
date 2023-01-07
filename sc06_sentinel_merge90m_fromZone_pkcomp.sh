#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc06_sentinel_merge90m_fromZone_pkcomp.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc06_sentinel_merge90m_fromZone_pkcomp.sh.%J.err
#SBATCH --job-name=sc06_sentinel_merge90m_fromZone_pkcomp.sh 
#SBATCH --mem=40G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc06_sentinel_merge90m_fromZone_pkcomp.sh
source ~/bin/gdal3
source ~/bin/pktools 


ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_zone
OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_merge
RAM=/dev/shm

# pkcomposite  -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co ZLEVEL=9   -ulx 2 -uly 15 -lrx 15 -lry 4 -co COMPRESS=DEFLATE -co ZLEVEL=9 -cr mean \
# $( for file in  $IN/???.tif ; do echo "-i" $file ; done )    -o $OUT/tile_all_90m_fromZone.tif 

# gdal_edit.py -a_nodata 0 $OUT/tile_all_90m_fromZone.tif  
# gdalinfomm $OUT/tile_all_90m_fromZone.tif  

pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9   -f min -i $OUT/tile_all_90m_fromZone.tif -o  $OUT/tile_all_90m_fromZone_min_tmp.tif 
gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9  $OUT/tile_all_90m_fromZone_min_tmp.tif  $OUT/tile_all_90m_fromZone_min.tif
rm $OUT/tile_all_90m_fromZone_min_tmp.tif 
