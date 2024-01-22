#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc11_merge_iteration.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc11_merge_iteration.sh.%J.err
#SBATCH --job-name=sc11_merge_iteration.sh 
#SBATCH --mem=60G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc11_merge_iteration.sh 

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
source ~/bin/gdal3
source ~/bin/pktools

GDAL_CACHEMAX=30000

cd $ONCHO/
gdalbuildvrt  -separate -overwrite -b 1 $ONCHO/prediction_all/prediction_seedP_all_msk.vrt $ONCHO/prediction_seed?/prediction_seed*P_all_msk.tif $ONCHO/prediction_seed??/prediction_seed*P_all_msk.tif    $ONCHO/prediction_seed???/prediction_seed*P_all_msk.tif   

pkstatprofile -co COMPRESS=LZW -co ZLEVEL=9 -nodata -9999 -f mean -f stdev -i $ONCHO/prediction_all/prediction_seedP_all_msk.vrt -o $ONCHO/prediction_all/prediction_seedP_all_msk_tmp.tif

gdal_translate -b 1 -co COMPRESS=LZW -co ZLEVEL=9 $ONCHO/prediction_all/prediction_seedP_all_msk_tmp.tif $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  
gdal_translate -b 2 -co COMPRESS=LZW -co ZLEVEL=9 $ONCHO/prediction_all/prediction_seedP_all_msk_tmp.tif $ONCHO/prediction_all/prediction_seedP_all_msk_stdev.tif


gdal_translate -tr 0.0083333333333333 0.0083333333333333 -r average  -co COMPRESS=LZW -co ZLEVEL=9  $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif   $ONCHO/prediction_all/prediction_seedP_all_msk_mean_1km.tif
gdal_translate -tr 0.0083333333333333 0.0083333333333333 -r average -co COMPRESS=LZW -co ZLEVEL=9 $ONCHO/prediction_all/prediction_seedP_all_msk_stdev.tif $ONCHO/prediction_all/prediction_seedP_all_msk_stdev_1km.tif


exit 

gdal_translate  -ot Byte    -co  COMPRESS=LZW  -co ZLEVEL=9 -scale $(pkstat -nodata -9999  -min -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif | awk '{ print $2  }')  $(pkstat -max -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif | awk '{ print $2  }') 1 255 $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif $ONCHO/prediction_all/prediction_seedP_all_msk_mean_scale.tif

pkcreatect -min 0 -max 255 > /tmp/color.txt 

pkcreatect -ot Byte -co COMPRESS=DEFLATE -co ZLEVEL=9 -ct /tmp/color.txt -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean_scale.tif -o $ONCHO/prediction_all/prediction_seedP_all_msk_mean_scale_ct.tif

gdal_translate -b 1   -co ZLEVEL=9    -of  MBTiles  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_scale.tif  $ONCHO/prediction_all/prediction_seedP_all_msk_mean.MBTiles

pkgetmask -ot Byte -co  COMPRESS=LZW  -co ZLEVEL=9  -min 0.6 -max 1            -data 60 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_60.tif 
pkgetmask -ot Byte -co  COMPRESS=LZW  -co ZLEVEL=9  -min 0.5 -max 0.599999999  -data 50 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_50.tif 
pkgetmask -ot Byte -co  COMPRESS=LZW  -co ZLEVEL=9  -min 0.4 -max 0.499999999  -data 40 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_40.tif 
pkgetmask -ot Byte -co  COMPRESS=LZW  -co ZLEVEL=9  -min 0.3 -max 0.399999999  -data 30 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_30.tif 
pkgetmask -ot Byte -co  COMPRESS=LZW  -co ZLEVEL=9  -min 0.2 -max 0.299999999  -data 20 -nodata 0 -i $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -o  $ONCHO/prediction_all/prediction_seedP_all_msk_mean_20.tif 




