#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc06_sentinel_histeq90m_grass.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc06_sentinel_histeq90m_grass.sh.%J.err
#SBATCH --job-name=sc06_sentinel_histeq90m_grass.sh
#SBATCH --mem=50G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc06_sentinel_histeq90m_grass.sh

#### for tif in  *.tif ; do gdalinfo -mm  $tif  | grep Computed | awk '{ gsub(/[=,]/," " , $0 ); print int($4) }'  ; done  | sort -g | tail    10800 
####  

source ~/bin/gdal3
source ~/bin/pktools



export ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
export SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
export  OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_tile_histeq
export   IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_tile
export  MERGE=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_merge_histeq

# apptainer exec  ~/bin/grass8.sif bash <<'EOF'

# rm -fr /tmp/mylocation 
# grass -f --text  -c $IN/all_tif.vrt -e /tmp/mylocation 

# # -ulx 2 -uly 15 -lrx 15 -lry 4 
# # grass /tmp/mylocation/PERMANENT --exec  g.region n=15  s=4 e=25  w=2

# #### g.extension extension=i.histo.match
# for file in $IN/*.tif ; do
# grass /tmp/mylocation/PERMANENT --exec  r.external  input=$file   output=$(basename  $file .tif)  --overwrite
# done

# cd /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_tile
# for band  in 1 2 3 4 ; do

# grass /tmp/mylocation/PERMANENT --exec i.histo.match input=$(ls *.tif  | tr "\n" " "  | sed  "s/.tif /.$band,/g" | sed 's/.$//' )  max=11000 --overwrite --verbose

# for file in $IN/*.tif ; do
# filename=$(basename $file .tif )
# grass /tmp/mylocation/PERMANENT --exec r.out.gdal -f --o -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=UInt16 format=GTiff nodata=0 input=${filename}.${band}.match output=$OUT/${filename}_${band}_match.tif

# done 
# done 

# EOF

rm -fr /tmp/mylocation 

for band  in 1 2 3 4 ; do
pkcomposite  -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co ZLEVEL=9 -ulx 2 -uly 15 -lrx 15 -lry 4  -cr median $( for file in  $OUT/5*_${band}_match.tif ; do echo "-i" $file ; done )  -o $MERGE/mosaic${band}_match.tif
done 

gdalbuildvrt -overwrite  -separate $MERGE/mosaic_all.vrt $MERGE/mosaic1_match.tif $MERGE/mosaic2_match.tif $MERGE/mosaic3_match.tif $MERGE/mosaic4_match.tif
gdal_translate   -co COMPRESS=DEFLATE -co ZLEVEL=9 $MERGE/mosaic_all.vrt  $MERGE/mosaic_all.tif
