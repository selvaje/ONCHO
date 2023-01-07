#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc07_sentinel_merge90m_index.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc07_sentinel_merge90m_index.sh.%J.err
#SBATCH --job-name=sc07_sentinel_merge90m_index.sh 
#SBATCH --mem=40G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc07_sentinel_merge90m_index.sh

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
SENT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel
 IN=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_merge
OUT=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/sentinel/tif_wgs84_90m_index
RAM=/dev/shm

# conda create -n  gdalcalc_env 

module load miniconda/4.9.2
source /gpfs/gibbs/project/sbsc/ga254/conda_envs/gdalcalc_env/lib/python3.10/venv/scripts/common/activate

# Band 1 Block=256x256 Type=UInt16, ColorInterp=Gray
#   Description = red_p25
#   NoData Value=0
# Band 2 Block=256x256 Type=UInt16, ColorInterp=Undefined
#   Description = green_p25
#   NoData Value=0
# Band 3 Block=256x256 Type=UInt16, ColorInterp=Undefined
#   Description = blue_p25
#   NoData Value=0
# Band 4 Block=256x256 Type=UInt16, ColorInterp=Undefined
#   Description = nir_p25
#   NoData Value=0
# https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel-2/indexdb/

GDAL_CACHEMAX=40000
####  NDVI = (NIR-RED)/(NIR+RED) 
/gpfs/gibbs/project/sbsc/ga254/conda_envs/gdalcalc_env/bin/gdal_calc.py  --type=Float32 -A $IN/tile_all_90m_fromZone.tif --A_band=4 -B $IN/tile_all_90m_fromZone.tif --B_band=1 \
 --co=COMPRESS=DEFLATE --co=ZLEVEL=9   --overwrite  \
--calc="( (B.astype(float32)  - A.astype(float32)) /  ((B.astype(float32)  + A.astype(float32) + 0.01) )    )" --outfile=$OUT/ndvi_tmp.tif 

#### Normalized Difference Water Index (NDWI)  = gren - near / gren + nir 

/gpfs/gibbs/project/sbsc/ga254/conda_envs/gdalcalc_env/bin/gdal_calc.py  --type=Float32 -A $IN/tile_all_90m_fromZone.tif --A_band=2 -B $IN/tile_all_90m_fromZone.tif --B_band=4 \
 --co=COMPRESS=DEFLATE --co=ZLEVEL=9   --overwrite  \
--calc="( (B.astype(float32)  - A.astype(float32)) /  ((B.astype(float32)  + A.astype(float32) + 0.01) )    )" --outfile=$OUT/ndwi_tmp.tif 

####  Normalized Difference Turbidity Index (NDTI) = red - gren  / red + gren 
####  to estimate the turbidity of water (Kuhn et al. 2019; Lacaux et al. 2007; Pahlevan et al. 2019):

/gpfs/gibbs/project/sbsc/ga254/conda_envs/gdalcalc_env/bin/gdal_calc.py  --type=Float32 -A $IN/tile_all_90m_fromZone.tif --A_band=1 -B $IN/tile_all_90m_fromZone.tif --B_band=2 \
 --co=COMPRESS=DEFLATE --co=ZLEVEL=9   --overwrite  \
--calc="( (B.astype(float32)  - A.astype(float32)) /  ((B.astype(float32)  + A.astype(float32) + 0.01) )    )" --outfile=$OUT/ndti_tmp.tif 


conda  deactivate

source ~/bin/gdal3
source ~/bin/pktools
pksetmask   -co COMPRESS=DEFLATE  -co ZLEVEL=9 -m $IN/tile_all_90m_fromZone_min.tif -msknodata 0 -nodata -9 -i $OUT/ndvi_tmp.tif -o $OUT/ndvi.tif  
rm -f $OUT/ndvi_tmp.tif
gdal_edit.py -a_nodata  -9 $OUT/ndvi.tif 
gdalinfo -mm  $OUT/ndvi.tif     | grep Computed | awk '{ gsub(/[=,]/," " , $0 ); print $3 , $4 }' >  $OUT/ndvi.mm
gdal_edit.py -a_nodata -19 $OUT/ndvi.tif # fake the nodata to facilitate the ingestion in R

pksetmask   -co COMPRESS=DEFLATE  -co ZLEVEL=9 -m $IN/tile_all_90m_fromZone_min.tif -msknodata 0 -nodata -9 -i $OUT/ndwi_tmp.tif -o $OUT/ndwi.tif  
rm -f $OUT/ndwi_tmp.tif
gdal_edit.py -a_nodata  -9 $OUT/ndwi.tif 
gdalinfo -mm  $OUT/ndwi.tif     | grep Computed | awk '{ gsub(/[=,]/," " , $0 ); print $3 , $4 }' >  $OUT/ndwi.mm
gdal_edit.py -a_nodata -19 $OUT/ndwi.tif 

pksetmask   -co COMPRESS=DEFLATE  -co ZLEVEL=9 -m $IN/tile_all_90m_fromZone_min.tif -msknodata 0 -nodata -9 -i $OUT/ndti_tmp.tif -o $OUT/ndti.tif  
rm -f $OUT/ndti_tmp.tif
gdal_edit.py -a_nodata  -9 $OUT/ndti.tif 
gdalinfo -mm  $OUT/ndti.tif     | grep Computed | awk '{ gsub(/[=,]/," " , $0 ); print $3 , $4 }' >  $OUT/ndti.mm
gdal_edit.py -a_nodata -19 $OUT/ndti.tif 

