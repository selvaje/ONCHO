#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc06_sentinel_merge90m_from10m_otb.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc06_sentinel_merge90m_from10m_otb.sh.%J.err
#SBATCH --job-name=sc06_sentinel_merge90m_from10m_otb.sh 
#SBATCH --mem=120G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc06_sentinel_merge90m_from10m_otb.sh 

source ~/bin/gdal3



ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m
OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_from10m_otb
RAM=/dev/shm

# pkcomposite  -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co ZLEVEL=9   -ulx 2 -uly 15 -lrx 15 -lry 4 -co COMPRESS=DEFLATE -co ZLEVEL=9 -cr mean \
# $( for file in  $IN/???.tif ; do echo "-i" $file ; done )    -o $IN/tile_all_90m.tif 


##### oft 
##### cd ~/bin
##### apptainer build  orfeotoolbox.sif   docker://orfeotoolbox/otb

apptainer exec  ~/bin/orfeotoolbox.sif bash <<EOF
export GDAL_CACHEMAX=80000
export OTB_MAX_RAM_HINT=80000

otbcli_Mosaic -interpolator linear  -output.spacingx 0.0008333333333333  -output.spacingy 0.0008333333333333 -nodata 0  -comp.feather large -il  /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_10m/*.tif  -out $OUT/tile_all_90m_otb1.tif
gdal_translate -ot  UInt16 -projwin 2 15 15 4 -a_ullr 2 15 15 4  -co COMPRESS=DEFLATE -co ZLEVEL=9   $OUT/tile_all_90m_otb1.tif   $OUT/tile_all_90m_otb_large_from10m.tif
rm $OUT/tile_all_90m_otb1.tif

otbcli_Mosaic -interpolator linear  -output.spacingx 0.0008333333333333  -output.spacingy 0.0008333333333333 -nodata 0  -comp.feather slim -comp.feather.slim.length 0.008333333333333 -il /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_10m/*.tif  -out $OUT/tile_all_90m_otb2.tif
gdal_translate -ot  UInt16 -projwin 2 15 15 4 -a_ullr 2 15 15 4  -co COMPRESS=DEFLATE -co ZLEVEL=9   $OUT/tile_all_90m_otb2.tif   $OUT/tile_all_90m_otb_slim_from10m.tif
rm $OUT/tile_all_90m_otb2.tif

otbcli_Mosaic -interpolator linear  -output.spacingx 0.0008333333333333  -output.spacingy 0.0008333333333333 -nodata 0 -harmo.method band -comp.feather slim -comp.feather.slim.length  0.008333333333333 -il /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_10m/*.tif -out $OUT/tile_all_90m_otb3.tif
gdal_translate -ot  UInt16 -projwin 2 15 15 4 -a_ullr 2 15 15 4  -co COMPRESS=DEFLATE -co ZLEVEL=9   $OUT/tile_all_90m_otb3.tif   $OUT/tile_all_90m_otb_slim_harm_from10m.tif
rm $OUT/tile_all_90m_otb3.tif
EOF






