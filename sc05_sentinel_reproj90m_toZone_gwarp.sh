#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc05_sentinel_reproj90m_toZone_gwarp.sh.%A_%a.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc05_sentinel_reproj90m_toZone_gwarp.sh.%A_%a.err
#SBATCH --job-name=sc05_sentinel_reproj90m_toZone_gwarp.sh
#SBATCH --mem=50G
#SBATCH --array=1-6

##   --array=1-4

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc05_sentinel_reproj90m_toZone_gwarp.sh
#### perform warping for each UTM zone 

source ~/bin/gdal3

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_orig
OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_zone

# SLURM_ARRAY_TASK_ID=1
zone=$(ls $SENT/shp_tile/percentile*.gpkg   | head -$SLURM_ARRAY_TASK_ID | tail -1 )
zonename=$(basename $file .gpkg)

zonen=$(ogrinfo  -al  -geom=NO $zone  | grep location | awk '{ if (NR==2) print substr( $4, length($4)-28 , 3  )   }')

gdalbuildvrt  -overwrite  $IN/$zonen.vrt   $(ogrinfo  -al  -geom=NO $zone  | grep location | awk '{ if(NR>1)  print "/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_orig/"$4  }')

GDAL_CACHEMAX=40000
gdalwarp -overwrite  -wm 40000000  -co COMPRESS=DEFLATE -co ZLEVEL=9  -r bilinear -t_srs EPSG:4326  -tr  0.00083333333333333333333333 0.00083333333333333333333333 $IN/$zonen.vrt $OUT/$zonen.tif


exit 

# remove below 
   # xmin ymin xmax ymax
   #   2   4    15   15 
echo "NW  2  10     8   15"  > $OUT/zone.txt  
echo "NE  8  10    15   15" >> $OUT/zone.txt  
echo "SW  2  4      8   10" >> $OUT/zone.txt  
echo "SE  8  4     15   10" >> $OUT/zone.txt  

zone=$( awk -v ID=SLURM_ARRAY_TASK_ID '{ if (NR==ID) print $1  }' $OUT/zone.txt  )
xmin=$( awk -v ID=SLURM_ARRAY_TASK_ID '{ if (NR==ID) print $2  }' $OUT/zone.txt  )
ymin=$( awk -v ID=SLURM_ARRAY_TASK_ID '{ if (NR==ID) print $3  }' $OUT/zone.txt  )
xmax=$( awk -v ID=SLURM_ARRAY_TASK_ID '{ if (NR==ID) print $4  }' $OUT/zone.txt  )
ymax=$( awk -v ID=SLURM_ARRAY_TASK_ID '{ if (NR==ID) print $5  }' $OUT/zone.txt  )


# SLURM_ARRAY_TASK_ID=1
GDAL_CACHEMAX=40000
gdalwarp -overwrite  -wm 40000000  -co COMPRESS=DEFLATE -co ZLEVEL=9  -r average -t_srs EPSG:4326 -te $xmin $ymin $xmax $ymax  -tr  0.00083333333333333333333333 0.00083333333333333333333333 $IN/$zonen.vrt $OUT/$zonen.tif
