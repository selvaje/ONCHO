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
gdaltindex $SENT/shp_zone/$tile.shp ${tile}_UTM.vrt
done 

cd $SENT/tif_orig

for tile in 31N 31P 32N 32P 33N 33P ; do 
wget -A tif  -r -nH --cut-dirs=7 --no-parent --reject="index.html*"   https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_composite_S2_L1C_2017-2018_GLOBE_R2020A/GHS_composite_S2_L1C_2017-2018_GLOBE_R2020A_UTM_10/V1-0/$tile/ 
done 

### the 30N appear after the download 
cd $SENT/tif_orig
for tile in 30N  31N 31P 32N 32P 33N 33P ; do
gdaltindex $SENT/shp_tile/percentile_30_$tile.gpkg $( for file in *.tif ; do echo $file $(gdalinfo $file | grep PROJCRS) ; done | grep S2_percentile_30_UTM | grep $tile | awk '{print $1 }')
gdaltindex $SENT/shp_tile/percentile_$tile.gpkg $( for file in *.tif ; do echo $file $(gdalinfo $file | grep PROJCRS) ; done | grep S2_percentile_UTM | grep $tile | awk '{print $1 }')
done 

exit 
# rename by hand   to be consistent. 

mv percentile_30N.gpkg percentile_31P.gpkg
mv percentile_30_30N.gpkg percentile_31N_new.gpkg
mv percentile_30_32N.gpkg percentile_32N.gpkg
mv percentile_31N.gpkg percentile_32P.gpkg
mv percentile_31N_new.gpkg percentile_31N.gpkg
mv percentile_30_33N.gpkg  percentile_33N.gpkg


