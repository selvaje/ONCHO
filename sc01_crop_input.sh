#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 10:00:00
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc01_crop_input.sh.%J.out
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc01_crop_input.sh.%J.err
#SBATCH --job-name=sc01_crop_input.sh
#SBATCH --mem=40G

module load StdEnv
source ~/bin/gdal3
source ~/bin/pktools
             
##### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc01_crop_input.sh

MERIT=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT/geomorphometry_90m_wgs84
ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO


#### geomorpho90m detaset
#### instruction for download in https://spatial-ecology.net/docs/build/html/GEODATA/geomorpho90m/geomorpho90m.html
#### use of a global vrt (store on local server)  to crop nigeria. In case of downloading identify the Nigeria tiles

for TOPO in geom aspect aspect-sine cti dev_scale dxx dy eastness pcurv roughness slope tcurv tri aspect-cosine convergence dev_magnitude dx dxy dyy elev-stdev northness rough-magnitude rough-scale spi tpi vrm ; do 
     gdalbuildvrt -overwrite $MERIT/${TOPO}/all_${TOPO}_90M.vrt $MERIT/${TOPO}/${TOPO}_90M_???????.tif 
     gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $MERIT/${TOPO}/all_${TOPO}_90M.vrt $ONCHO/input/geomorpho90m/${TOPO}.tif 
done

#### change - to _ to avoid issues with R

mv $ONCHO/input/geomorpho90m/elev-stdev.tif   $ONCHO/input/geomorpho90m/elev_stdev.tif 
mv $ONCHO/input/geomorpho90m/aspect-sine.tif   $ONCHO/input/geomorpho90m/aspect_sine.tif 
mv $ONCHO/input/geomorpho90m/aspect-cosine.tif   $ONCHO/input/geomorpho90m/aspect_cosine.tif 
mv $ONCHO/input/geomorpho90m/rough-magnitude.tif   $ONCHO/input/geomorpho90m/rough_magnitude.tif 
mv $ONCHO/input/geomorpho90m/rough-scale.tif   $ONCHO/input/geomorpho90m/rough_scale.tif
#### remove cti ans spi better use the one in hydrogrphy90m 
rm -f $ONCHO/input/geomorpho90m/aspect.tif $ONCHO/input/geomorpho90m/geom.tif $ONCHO/input/geomorpho90m/cti.tif $ONCHO/input/geomorpho90m/spi.tif

#### reclass -9999 (water areas) to 0 

pkreclass -c -9999 -r 0 -co COMPRESS=DEFLATE -co ZLEVEL=9  -i   $ONCHO/input/geomorpho90m/aspect_cosine.tif  -o  $ONCHO/input/geomorpho90m/aspect_cosine_tmp.tif 
mv $ONCHO/input/geomorpho90m/aspect_cosine_tmp.tif   $ONCHO/input/geomorpho90m/aspect_cosine.tif 

pkreclass -c -9999 -r 0 -co COMPRESS=DEFLATE -co ZLEVEL=9  -i   $ONCHO/input/geomorpho90m/aspect_sine.tif  -o  $ONCHO/input/geomorpho90m/aspect_sine_tmp.tif 
mv $ONCHO/input/geomorpho90m/aspect_sine_tmp.tif   $ONCHO/input/geomorpho90m/aspect_sine.tif

### download elevation data from the MERIT DEM https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/
### the tiles for nigeria are  n00e005_dem.tif n00e010_dem.tif n05e000_dem.tif n05e005_dem.tif   n05e010_dem.tif n10e000_dem.tif  n10e005_dem.tif   n10e010_dem.tif
### download the corrispondent *.tar.gz

gdal_translate -projwin 2 15 15 4  -co COMPRESS=DEFLATE -co ZLEVEL=9 /gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT/input_tif/all_tif.vrt  $ONCHO/input/geomorpho90m/elevation.tif


#### hydrography90m detaset
#### instruction for download in https://hydrography.org/hydrography90m/hydrography90m_layers & https://hydrography.org/hydrography90m/hydrography90m_download_script

#### use of a global vrt (store on local server)  to crop nigeria. In case of downloading identify the Nigeria tiles

HYDRO=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT_HYDRO/hydrography90m_v.1.0

for HYDROD in   flow.index  r.stream.distance  r.stream.slope  r.watershed  ; do 
    for DIR in $(ls $HYDRO/$HYDROD ) ; do 
     gdal_translate  -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $HYDRO/$HYDROD/$DIR/$(basename $DIR _tiles20d).vrt /tmp/$(basename $DIR _tiles20d).tif 
     pksetmask  -co COMPRESS=DEFLATE -co ZLEVEL=9 -m /tmp/$(basename $DIR _tiles20d).tif -msknodata -9999  -nodata 0 -i  /tmp/$(basename $DIR _tiles20d).tif   -o $ONCHO/input/hydrography90m/$(basename $DIR _tiles20d).tif 
     rm /tmp/$(basename $DIR _tiles20d).tif 
     gdal_edit.py -a_nodata -9999 $ONCHO/input/hydrography90m/$(basename $DIR _tiles20d).tif 
    done 
done 

#### remove hydrography90m variables that do not have a continues information (only value at stream level)

rm -f $ONCHO/input/hydrography90m/channel_*.tif $ONCHO/input/hydrography90m/order_*.tif $ONCHO/input/hydrography90m/{segment.tif,regional_unit.tif,basin.tif,depression.tif,direction.tif,sub_catchment.tif,outlet.tif}

#### reclass no data value 

pkreclass -c -9999999  -r -2147483648  -co COMPRESS=DEFLATE -co ZLEVEL=9  -i   $ONCHO/input/hydrography90m/accumulation.tif  -o  $ONCHO/input/hydrography90m/accumulation_tmp.tif 
mv $ONCHO/input/hydrography90m/accumulation_tmp.tif    $ONCHO/input/hydrography90m/accumulation.tif 

pkreclass -c -9999999  -r -2147483648 -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $ONCHO/input/hydrography90m/slope_curv_max_dw_cel.tif -o $ONCHO/input/hydrography90m/slope_curv_max_dw_cel_tmp.tif
mv $ONCHO/input/hydrography90m/slope_curv_max_dw_cel_tmp.tif  $ONCHO/input/hydrography90m/slope_curv_max_dw_cel.tif

pkreclass -c -9999999  -r -2147483648 -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $ONCHO/input/hydrography90m/slope_curv_min_dw_cel.tif -o $ONCHO/input/hydrography90m/slope_curv_min_dw_cel_tmp.tif
mv $ONCHO/input/hydrography90m/slope_curv_min_dw_cel_tmp.tif  $ONCHO/input/hydrography90m/slope_curv_min_dw_cel.tif

pkreclass -c -9999999  -r -2147483648 -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $ONCHO/input/hydrography90m/slope_grad_dw_cel.tif -o $ONCHO/input/hydrography90m/slope_grad_dw_cel_tmp.tif
mv $ONCHO/input/hydrography90m/slope_grad_dw_cel_tmp.tif  $ONCHO/input/hydrography90m/slope_grad_dw_cel.tif 

#### CHELSA dataset for climatic data
#### instruction for download at https://chelsa-climate.org/downloads/

#### for the annual species distribution I used the Bioclimatic variables   https://chelsa-climate.org/bioclim/
#### in case of monthly or seasonal species distribution consider the monthly precipitation and temperature  
        
CHELSA=/gpfs/gibbs/pi/hydro/hydro/dataproces/CHELSA/climatologies/bio

for CHELSAD  in $CHELSA/CHELSA_bio*_1981-2010_V.2.1.tif  ; do
     filename=$(basename $CHELSAD _1981-2010_V.2.1.tif )
     gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $CHELSAD  $ONCHO/input/chelsa/${filename}.tif
     gdal_translate   -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/input/chelsa/${filename}.tif  $ONCHO/input/chelsa/${filename}_r.tif
done


#### rm the following bio clim that have low variability data (constant value in large areas). 
rm   $ONCHO/input/chelsa/CHELSA_bio9_r.tif  $ONCHO/input/chelsa/CHELSA_bio9.tif 
rm   $ONCHO/input/chelsa/CHELSA_bio19_r.tif  $ONCHO/input/chelsa/CHELSA_bio19.tif
rm   $ONCHO/input/chelsa/CHELSA_bio8_r.tif  $ONCHO/input/chelsa/CHELSA_bio8.tif  
rm   $ONCHO/input/chelsa/CHELSA_bio18_r.tif  $ONCHO/input/chelsa/CHELSA_bio18.tif  

#### Soil temperature https://onlinelibrary.wiley.com/doi/10.1111/gcb.16060 
#### files at https://zenodo.org/record/4558732#.YhZQoHWYWV5

SOILT=/gpfs/gibbs/pi/hydro/hydro/dataproces/SOILTEMP/input

for SOIL  in $( ls  $SOILT/*.tif | grep -v msk )  ; do
filename=$(basename  $SOIL .tif)
gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $SOIL  $ONCHO/input/soiltemp/${filename}.tif
pkgetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -min  -9999999  -max 9999999 -data 1 -nodata 0   -i $ONCHO/input/soiltemp/${filename}.tif -o $ONCHO/input/soiltemp/${filename}_msk01.tif
pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $ONCHO/input/soiltemp/${filename}_msk01.tif -msknodata 0 -nodata -9999 -i $ONCHO/input/soiltemp/${filename}.tif -o $ONCHO/input/soiltemp/${filename}_msk.tif

#### fill the nodata that belongs to stream, taking the nearby pixel value

gdal_fillnodata.py  -co COMPRESS=DEFLATE -co ZLEVEL=9 -nomask  -md 40  -si 1  $ONCHO/input/soiltemp/${filename}_msk.tif $ONCHO/input/soiltemp/${filename}.tif 
gdal_translate   -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ONCHO/input/soiltemp/${filename}.tif  $ONCHO/input/soiltemp/${filename}_r.tif
done

gdalbuildvrt -overwrite -separate $ONCHO/input/soiltemp/all_tif_msk01.vrt   $ONCHO/input/soiltemp/*_msk01.tif
pkstatprofile -ot Byte  -of GTiff  -co COMPRESS=DEFLATE -co ZLEVEL=9 -f sum  -i $ONCHO/input/soiltemp/all_tif_msk01.vrt -o $ONCHO/input/soiltemp/soiltmp_msk_tmp.tif 
pkgetmask -ot Byte   -co COMPRESS=DEFLATE -co ZLEVEL=9  -min   0.5  -max 50  -data 1  -nodata 0   -i $ONCHO/input/soiltemp/soiltmp_msk_tmp.tif   -o $ONCHO/input/soiltemp/soiltmp_msk.tif 
rm $ONCHO/input/soiltemp/soiltmp_msk_tmp.tif $ONCHO/input/soiltemp/*_msk01.tif   $ONCHO/input/soiltemp/all_tif_msk01.vrt 


#### Soil type dataset.  https://soil.copernicus.org/articles/7/217/2021/soil-7-217-2021.html 
#### instruction for download at https://soilgrids.org/ 
#### donwload can be done using  
#### gdal_translate -projwin 1733000   528000   1734000 526000   -co COMPRESS=DEFLATE -co ZLEVEL=9 /vsicurl/https://files.isric.org/soilgrids/latest/data/clay/clay_0-5cm_mean.vrt  test.tif
#### Consider that the original data are in Goode's homolosine projection, so Nigeria's coordinates in meter need to be defined first. 

SOILGRIDS=/gpfs/gibbs/pi/hydro/hydro/dataproces/SOILGRIDS

for SOIL in $SOILGRIDS/*_acc/*_WeigAver.vrt ; do
     filename=$(basename  $SOIL  .vrt)
     gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $SOIL  $ONCHO/input/soilgrids/${filename}_acc.tif
     gdal_edit.py -a_nodata -9999  $ONCHO/input/soilgrids/soilgrids_msk.tif  # to fake the no-data
done

for SOIL in $(ls $SOILGRIDS/*/*_WeigAver.tif |  grep -v  -e _acc -e out_  )  ; do
filename=$(basename  $SOIL  .tif)
# gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $SOIL  $ONCHO/input/soilgrids/${filename}_tmp.tif
# gdal_fillnodata.py  -co COMPRESS=DEFLATE -co ZLEVEL=9 -nomask  -md 40  -si 1  $ONCHO/input/soilgrids/${filename}_tmp.tif  $ONCHO/input/soilgrids/${filename}.tif
rm -f $ONCHO/input/soilgrids/${filename}_tmp.tif
gdal_translate   -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ONCHO/input/soilgrids/${filename}.tif $ONCHO/input/soilgrids/${filename}_r.tif
done

pkgetmask -ot Byte -co COMPRESS=DEFLATE -co ZLEVEL=9 -min 65534 -max 65536 -data 0 -nodata 1 -i $ONCHO/input/soilgrids/AWCtS_WeigAver.tif -o $ONCHO/input/soilgrids/soilgrids_msk.tif
gdal_edit.py -a_nodata 0 $ONCHO/input/soilgrids/soilgrids_msk.tif

##### population density 
POP=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/population

#### population density from grid3 

https://data.grid3.org/maps/GRID3::grid3-nigeria-gridded-population-estimates-version-2-0/about 
cd /gpfs/loomis/project/sbsc/ga254/dataproces/ONCHO/input/population 
wget https://wopr.worldpop.org/download/495
unzip 495


#### populaton density from GHSL  https://ghsl.jrc.ec.europa.eu/download.php?ds=pop 

wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R9_C19.zip
wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C19.zip
wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C20.zip
wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C19.zip

unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C19.zip  
unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C20.zip  
unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R9_C19.zip
unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R9_C20.zip

gdalbuildvrt -overwrite $POP/pop_100m.vrt  $POP/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_*.tif

gdalwarp -overwrite -ot Int16   -te $(getCorners4Gwarp $POP/../geomorpho90m/elevation.tif ) -t_srs EPSG:4326 -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9 $POP/pop_100m.vrt $POP/GHSpop_90m.tif
gdal_edit.py -a_nodata -9999 $POP/GHSpop_90m.tif   

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m /tmp/$(basename $DIR _tiles20d).tif -msknodata -9999  -nodata 0 -i  /tmp/$(basename $DIR _tiles20d).tif   -o $ONCHO/input/hydrography90m/$(basename $DIR _tiles20d).tif 


### glw livestock
### manual downloa from https://dataverse.harvard.edu/dataverse/glw_3


LS=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/livestock
for file in $LS/5_??_2010_Da.tif  ; do 
gdal_translate -ot Int32   -projwin  $(getCorners4Gtranslate  $LS/../geomorpho90m/elevation.tif )  -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9 $file  $LS/$(basename $file .tif )_90m.tif 
pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $LS/../geomorpho90m/elevation.tif -msknodata -9999 -nodata -9999 \
                                            -m $LS/$(basename $file .tif )_90m.tif  -msknodata -0.5 -p "<"  -nodata -9999 \
-i $LS/$(basename $file .tif )_90m.tif -o $LS/$(basename $file .tif )_90m_msk.tif
gdal_edit.py  -a_nodata -32768 $LS/$(basename $file .tif )_90m_msk.tif  # faeka the nodata
done

mv 5_Bf_2010_Da_90m_msk.tif    LS_Bf.tif  # romeve it ... only 0
mv 5_Ch_2010_Da_90m_msk.tif    LS_Ch.tif
mv 5_Ct_2010_Da_90m_msk.tif    LS_Ct.tif
mv 5_Dk_2010_Da_90m_msk.tif    LS_Dk.tif
mv 5_Gt_2010_Da_90m_msk.tif    LS_Gt.tif
mv 5_Ho_2010_Da_90m_msk.tif    LS_Ho.tif
mv 5_Pg_2010_Da_90m_msk.tif    LS_Pg.tif
mv 5_Sh_2010_Da_90m_msk.tif    LS_Sh.tif


#### https://ecoregions.appspot.com/ https://academic.oup.com/bioscience/article/67/6/534/3102935 
## wget https://storage.googleapis.com/teow2016/Ecoregions2017.zip
## unzip Ecoregions2017.zip

gdal_rasterize    -a_nodata 65535 -ot UInt16   -te  $(getCorners4Gwarp   $ER/../geomorpho90m/elevation.tif )  -l "Ecoregions2017"  -a "ECO_ID"   -tr 0.00083333333333333 0.00083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ER/Ecoregions2017.shp  $ER/Ecoregions2017.tif 
mv $ER/Ecoregions2017.tif     $ER/ER2017.tif
gdal_edit.py -unsetnodata     $ER/ER2017.tif
 
##### landcover 
## wget https://lulctimeseries.blob.core.windows.net/lulctimeseriespublic/lc2021/lulc2021.zip

LC=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input_orig/landcover
unzip 

for zone in 33N 33P 31N 31P 32N 32P ; do
gdalwarp -overwrite    -co COMPRESS=DEFLATE -co ZLEVEL=9  -r near   -t_srs EPSG:4326  -tr  0.000083333333333333333333333 0.000083333333333333333333333 $LC/lc2021/${zone}_20210101-20220101.tif $LC/${zone}_wgs84.tif 
done 
gdalbuildvrt -overwrite    $LC/lc2021_wgs84.vrt  $LC/???_wgs84.tif

gdal_translate -projwin $(getCorners4Gtranslate $LC/../geomorpho90m/elevation.tif) -co COMPRESS=DEFLATE -co ZLEVEL=9 $LC/lc2021_wgs84.vrt $LC/../../input/landcover/lc2021_wgs84.tif # run it

gdal_translate    -projwin $(getCorners4Gtranslate $LC/../geomorpho90m/elevation.tif) -tr  0.00083333333333333333333333 0.0008333333\
3333333333333333 -r mode     -co COMPRESS=DEFLATE -co ZLEVEL=9  $LC/lc2021_wgs84.vrt   $LC/../../input/landcover/lc2021_wgs84_r.tif 
mv $LC/../../input/landcover/lc2021_wgs84_r.tif $LC/../../input/landcover/LC2021.tif 
gdal_edit.py -a_nodata 255    $LC/../../input/landcover/LC2021.tif

#### water occurence from https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GSWE/Aggregated/LATEST/
#### eventualy modify this
#### wget -nd  -r -A .tif https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GSWE/Aggregated/LATEST/${dir}/tiles

WATER=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/water 
#### gdalbuildvrt  -srcnodata 255 -vrtnodata 100  -overwrite   occurrence.vrt  occurrence_download/occurrence-*-??????????.tif 
gdal_translate  -projwin $(getCorners4Gtranslate $WATER/../geomorpho90m/elevation.tif) -tr  0.00083333333333333333333333 0.00083333333333333333333333 -r average  -co COMPRESS=DEFLATE -co ZLEVEL=9   /gpfs/gibbs/pi/hydro/hydro/dataproces/GSW/input/occurrence.vrt   $WATER/occurence.tif  
gdal_edit.py -unsetnodata $WATER/occurence.tif 

gdal_proximity.py -co COMPRESS=DEFLATE -co ZLEVEL=9 -values 0 -distunits PIXEL -ot UInt32 $WATER/occurence.tif $WATER/occurence_proximity.tif

#### FLO1km  https://www.nature.com/articles/sdata201852 

WATER=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/water

gdal_translate  -projwin $(getCorners4Gtranslate $WATER/../geomorpho90m/elevation.tif)  -co COMPRESS=DEFLATE -co ZLEVEL=9 /gpfs/gibbs/pi/hydro/hydro/dataproces/FLOW1k/FLO1K_mean_averaged_1960-2015.tif $WATER/flow_1km.tif

gdalwarp -overwrite  -te $(getCorners4Gwarp $WATER/../geomorpho90m/elevation.tif) -tr 0.00083333333333333333333333 0.00083333333333333333333333 -r max  -co COMPRESS=DEFLATE -co ZLEVEL=9 /gpfs/gibbs/pi/hydro/hydro/dataproces/FLOW1k/FLO1K_mean_averaged_1960-2015.tif $WATER/flow_mean.tif
gdal_edit.py -unsetnodata $WATER/flow_mean.tif

exit

#### uncillary data for other analyisis 

#### ecoregion not used anymore due to few values 
####  wget https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip # just few classis ... not useful 
## ER=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/ecoregions
## gdal_rasterize   -ot Byte   -te  $(getCorners4Gwarp   $ER/../geomorpho90m/elevation.tif )  -l "wwf_terr_ecos"  -a "G200_NUM"   -tr 0.0083333333333333 0.0083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ER/official/wwf_terr_ecos.shp  $ER/wwf_terr_ecosG200_NUM.tif 


##### lat long 
export ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/

apptainer exec  ~/bin/grass8.sif bash <<EOF
grass -f --text -c  -c $ONCHO/input/geomorpho90m/elevation.tif    -e /tmp/mylocation
grass /tmp/mylocation/PERMANENT --exec g.region n=15  s=4 e=15  w=2
grass /tmp/mylocation/PERMANENT --exec r.external  input=$ONCHO/input/geomorpho90m/elevation.tif   output=elv  --overwrite
grass /tmp/mylocation/PERMANENT --exec r.latlong -l  input=elv  output=latiy 
grass /tmp/mylocation/PERMANENT --exec r.latlong     input=elv  output=longx
grass /tmp/mylocation/PERMANENT --exec r.out.gdal -f --o -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Float32 format=GTiff nodata=0 input=longx output=$ONCHO/input/latlon/x.tif
grass /tmp/mylocation/PERMANENT --exec r.out.gdal -f --o -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" type=Float32 format=GTiff nodata=0 input=latiy output=$ONCHO/input/latlon/y.tif
EOF

### road map https://download.geofabrik.de/africa/nigeria.html 

OSM=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/openstreetmap/
cd /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/openstreetmap/
gdal_rasterize  -a_nodata 0 -ot Byte   -te  $(getCorners4Gwarp $OSM/../geomorpho90m/elevation.tif )  -l "gis_osm_roads_free_1"  -burn 1    -tr 0.00083333333333333 0.00083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $OSM/gis_osm_roads_free_1.shp  $OSM/gis_osm_roads_free_1.tif
gdal_edit.py -unsetnodata  $OSM/gis_osm_roads_free_1.tif

gdal_proximity.py -fixed-buf-val 1  -maxdist 5 -nodata 0 -co COMPRESS=DEFLATE -co ZLEVEL=9 -values 0.5 -distunits PIXEL -ot Byte $OSM/gis_osm_roads_free_1.tif   $OSM/gis_osm_roads_free_1_tmp.tif
gdal_edit.py -unsetnodata  $OSM/gis_osm_roads_free_1_tmp.tif

### run in another termina... 
source /gpfs/gibbs/project/sbsc/ga254/conda_envs/gdalcalc_env/lib/python3.10/venv/scripts/common/activate

/gpfs/gibbs/project/sbsc/ga254/conda_envs/gdalcalc_env/bin/gdal_calc.py --overwrite  --NoDataValue=0  --co=COMPRESS=DEFLATE --co=ZLEVEL=9 --co=BIGTIFF=YES  -B $OSM/gis_osm_roads_free_1.tif -A $OSM/gis_osm_roads_free_1_tmp.tif  --outfile=$OSM/gis_osm_roads_free_1_buf.tif    --calc="(A + B )"

### rasterize river network  OSM
OSM=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/openstreetmap
gdal_rasterize  -a_nodata 0 -ot Byte   -te  $(getCorners4Gwarp $OSM/../geomorpho90m/elevation.tif ) -l "gis_osm_water_a_free_1"   -burn 1  -tr 0.00083333333333333 0.00083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $OSM/gis_osm_water_a_free_1.shp  $OSM/gis_osm_water_a_free_1.tif

gdal_rasterize  -a_nodata 0 -ot Byte   -te  $(getCorners4Gwarp $OSM/../geomorpho90m/elevation.tif ) -l "gis_osm_waterways_free_1" -burn 1  -tr 0.00083333333333333 0.00083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $OSM/gis_osm_waterways_free_1.shp  $OSM/gis_osm_waterways_free_1.tif

### rasterize river network  NGA 

NGA=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/nga/
gdal_rasterize  -a_nodata 0 -ot Byte   -te  $(getCorners4Gwarp $NGA/../geomorpho90m/elevation.tif )  -l "NGA_rvrsl_1m_dcw"  -burn 1    -tr 0.00083333333333333 0.00083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $NGA/NGA_rvrsl_1m_dcw.shp  $NGA/NGA_rvrsl_1m_dcw.tif 

gdal_rasterize  -a_nodata 0 -ot Byte   -te  $(getCorners4Gwarp $NGA/../geomorpho90m/elevation.tif )  -l "NGA_rvrsl_1m_esri"  -burn 1    -tr 0.00083333333333333 0.00083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $NGA/NGA_rvrsl_1m_esri.shp  $NGA/NGA_rvrsl_1m_esri.tif 

