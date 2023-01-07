#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc05_sentinel_reproj90m_pertile.sh.%A_%a.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc05_sentinel_reproj90m_pertile.sh.%A_%a.err
#SBATCH --job-name=sc05_sentinel_reproj90m_pertile.sh
#SBATCH --mem=10G
#SBATCH --array=1-71

##   --array=1-6

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc05_sentinel_reproj90m_pertile.sh 

source ~/bin/gdal3

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_orig
OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_tile

# SLURM_ARRAY_TASK_ID=1
file=$(ls $IN/*.tif | head -$SLURM_ARRAY_TASK_ID | tail -1 )
filename=$(echo ${file: -29})

GDAL_CACHEMAX=8000
gdalwarp -overwrite -srcnodata 0 -dstnodata 0  -wm 8000000  -co COMPRESS=DEFLATE -co ZLEVEL=9  -r bilinear -t_srs EPSG:4326  -tr  0.00083333333333333333333333 0.00083333333333333333333333 $file  $OUT/$filename
