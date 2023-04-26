#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc20_site_selection.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc20_site_selection.sh.%J.err
#SBATCH --job-name=sc20_site_selection.sh
#SBATCH --mem=10G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc20_site_selection.sh

export ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
export OSM=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/openstreetmap
export POP=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/population
export NGA=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/nga
export SITEI=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/site/input
source ~/bin/gdal3
source ~/bin/pktools

GDAL_CACHEMAX=30000

cd $ONCHO/

gdal_rasterize -co COMPRESS=DEFLATE -co ZLEVEL=9 -ot Byte -a_nodata 0 -burn 1 -tr 0.000833333333333333 0.000833333333333333 -te 2 4 15 15 -a_srs EPSG:4326 $SITEI/Niger_Ento_Surveillance_Sites.gpkg   $SITEI/Niger_Ento_Surveillance_Sites.tif 

rm -fr /tmp/mylocation 
apptainer exec  ~/bin/grass8.sif bash <<EOF
grass -f --text -c  -c $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  -e /tmp/mylocation

grass /tmp/mylocation/PERMANENT --exec g.extension extension=r.cell.area

grass /tmp/mylocation/PERMANENT --exec r.external  input=$ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  output=prob   --overwrite

# #### road  criteria ######
grass /tmp/mylocation/PERMANENT --exec r.external  input=$OSM/gis_osm_roads_free_1.tif output=road   --overwrite
grass /tmp/mylocation/PERMANENT --exec r.buffer -z  input=road output=road_buf distances=2000                         # create a buffer: 1 road 2 buffer 
grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=road_buf output=$SITEI/road_buf.tif

##### river  criteria ######
grass /tmp/mylocation/PERMANENT --exec r.external  input=$NGA/NGA_rvrsl_1m_dcw.tif          output=nga_river_dcw   --overwrite
grass /tmp/mylocation/PERMANENT --exec r.external  input=$NGA/NGA_rvrsl_1m_esri.tif         output=nga_river_esri   --overwrite
grass /tmp/mylocation/PERMANENT --exec r.external  input=$OSM/gis_osm_waterways_free_1.tif  output=osm_waterways   --overwrite
grass /tmp/mylocation/PERMANENT --exec r.external  input=$OSM/gis_osm_water_a_free_1.tif      output=osm_water   --overwrite

grass /tmp/mylocation/PERMANENT --exec r.patch  -s input=nga_river_dcw,nga_river_esri,osm_waterways,osm_water   output=river  # merge all the water layers
grass /tmp/mylocation/PERMANENT --exec r.buffer   input=river output=river_buf distances=5000                          # create a buffer: 1 stream 2 buffer 
grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=river output=$SITEI/river.tif
grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=river_buf output=$SITEI/river_buf.tif

###### population criteria , comunity with  ######
grass /tmp/mylocation/PERMANENT --exec r.external  input=$POP/NGA_population_v2_0_gridded.tif   output=pop   --overwrite
grass /tmp/mylocation/PERMANENT --exec r.mapcalc "pop_1 = if (pop > 0, 1 , null() )"                         --overwrite
grass /tmp/mylocation/PERMANENT --exec r.clump input=pop_1  output=pop_clump                                 --overwrite
grass /tmp/mylocation/PERMANENT --exec r.stats.zonal base=pop_clump  cover=pop  output=clump_pop_sum  method=sum --overwrite

grass /tmp/mylocation/PERMANENT --exec r.cell.area output=cellarea  units=m2
grass /tmp/mylocation/PERMANENT --exec r.stats.zonal base=pop_clump  cover=cellarea   output=clump_area_sum   method=sum --overwrite

grass /tmp/mylocation/PERMANENT --exec r.mapcalc "pop_comunity = if ((clump_pop_sum / (clump_area_sum * 1000000)) => 300, (clump_pop_sum / ( clump_area_sum * 1000000)), null())" --overwrite
grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=UInt32 format=GTiff input=pop_comunity output=$SITEI/pop_comunity.tif
grass /tmp/mylocation/PERMANENT --exec r.buffer   input=pop_comunity  output=pop_comunity_buf  distances=5000                       # create a buffer: 1 pop  2 buffer 
grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=UInt32 format=GTiff input=pop_comunity_buf output=$SITEI/pop_comunity_buf.tif

#### ento site criteria > 10 km  ####### 
grass /tmp/mylocation/PERMANENT --exec r.external  input=$SITEI/Niger_Ento_Surveillance_Sites.tif output=ento   --overwrite
grass /tmp/mylocation/PERMANENT --exec r.buffer   input=ento  output=ento_buf  distances=10000                  --overwrite
grass /tmp/mylocation/PERMANENT --exec r.mapcalc "ento_bufinv  = if ( isnull(ento_buf) , 1 , null() )"          --overwrite

# grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=ento_buf     output=$SITEI/ento_buf.tif 
# grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=ento_bufinv  output=$SITEI/ento_bufinv.tif 


#### probability  criteria ####### 
grass /tmp/mylocation/PERMANENT --exec r.mapcalc "prob5  = if (  prob >= 0.5 && prob < 0.6 , 1 , null() )"             --overwrite
grass /tmp/mylocation/PERMANENT --exec r.mapcalc "prob6  = if (  prob >= 0.6 && prob < 0.7 , 1 , null() )"             --overwrite

####### merge criterias  #### 

grass /tmp/mylocation/PERMANENT --exec r.mapcalc "criterias_prob5  = if ((  road_buf + river_buf + pop_comunity_buf + ento_bufinv + prob5  ) > 0 , 1 , null() ) "  --overwrite 
grass /tmp/mylocation/PERMANENT --exec r.mapcalc "criterias_prob6  = if ((  road_buf + river_buf + pop_comunity_buf + ento_bufinv + prob6  ) > 0 , 1 , null() ) "  --overwrite 

grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=criterias_prob5 output=$SITEI/criterias_prob5.tif 
grass /tmp/mylocation/PERMANENT --exec r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=0 type=Byte format=GTiff input=criterias_prob6 output=$SITEI/criterias_prob6.tif

### https://grass.osgeo.org/grass82/manuals/v.random.html 


EOF

exit 




