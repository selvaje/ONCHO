#!/bin/bash
#SBATCH -p transfer
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stderr/sc04_wget_sentinel.sh.%J.out  
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stdout/sc04_wget_sentinel.sh.%J.err
#SBATCH --job-name=sc01_crop_input.sh
#SBATCH --mem=5G


#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc04_wget_sentinel.sh 

source ~/bin/gdal3
source ~/bin/pktools


ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel


### source https://ghsl.jrc.ec.europa.eu/download.php?ds=compositeS2 


cd $SENT/vrt_orig

for tile in 31N 31P 32N 32P 33N 33P ; do 
wget -A vrt  -r -nH --cut-dirs=7 --no-parent --reject="index.html*"  https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_composite_S2_L1C_2017-2018_GLOBE_R2020A/GHS_composite_S2_L1C_2017-2018_GLOBE_R2020A_UTM_10/V1-0/$tile/
rm -f $SENT/shp_zone/$tile.shp
gdaltindex $SENT/shp_zone/$tile.shp $tile.vrt
done 

cd $SENT/tif_orig

for tile in 31N 31P 32N 32P 33N 33P ; do 
wget -A tif  -r -nH --cut-dirs=7 --no-parent --reject="index.html*"   https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_composite_S2_L1C_2017-2018_GLOBE_R2020A/GHS_composite_S2_L1C_2017-2018_GLOBE_R2020A_UTM_10/V1-0/$tile/ 
gdaltindex $SENT/shp_tile/$tile.shp  *.tif
done 

exit 
