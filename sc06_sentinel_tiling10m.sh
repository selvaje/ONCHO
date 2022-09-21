#!/bin/bash 
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc06_sentinel_tiling10m.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc06_sentinel_tiling10m.sh.%J.err
#SBATCH --job-name=sc06_sentinel_tiling10m.sh
#SBATCH --mem=40G
#SBATCH --array=1-41

##   --array=1-41

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc06_sentinel_tiling10m.sh

source ~/bin/gdal3
source ~/bin/pktools 

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_10m
OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_10m_tile

#  x 2 15 
#  y 4 15                                                                                                              remove the first row that includes only see area
# for x in $(seq 2 2 14 )  ; do for y in $(seq 4 2 14 ) ; do echo $x $(expr $x + 2 ) $y $(expr $y + 2 ) ; done ; done  | awk '{ if (NR>1) print  }'  >   $ONCHO/vector/tile_list.txt 

tile=$( head -$SLURM_ARRAY_TASK_ID  /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/tile_list.txt | tail -1 )

ulx=$( echo $tile | awk '{ print $1 }' )
uly=$( echo $tile | awk '{ print $4 }' )
lrx=$( echo $tile | awk '{ print $2 }' ) 
lry=$( echo $tile | awk '{ print $3 }' )

pkcomposite -co COMPRESS=DEFLATE -co ZLEVEL=9 -co  BIGTIFF=YES  -srcnodata 0 -dstnodata 0 -ulx $ulx -uly $uly -lrx $lrx -lry $lry -co COMPRESS=DEFLATE -co ZLEVEL=9 -cr mean \
$( for file in  $IN/*.tif ; do echo "-i" $file ; done )    -o $OUT/tile_${ulx}_${uly}_${lrx}_${lry}.tif 

