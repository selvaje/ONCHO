#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 10:00:00 
#SBATCH -o /gpfs/scratch60/fas/sbsc/ga254/stdout/sc01_crop_input.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/sbsc/ga254/stderr/sc01_crop_input.sh.%J.err
#SBATCH --job-name=sc01_crop_input.sh
#SBATCH --mem=40G

source ~/bin/gdal3
source ~/bin/pktools

##### sbatch /vast/palmer/home.grace/ga254/scripts/ONCHO/sc01_crop_input.sh

MERIT=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT/geomorphometry_90m_wgs84
ONCHO=/gpfs/loomis/project/sbsc/ga254/dataproces/ONCHO

# for TOPO in geom aspect aspect-sine cti dev_scale dxx dy eastness pcurv roughness slope tcurv tri aspect-cosine convergence dev_magnitude dx dxy dyy elev-stdev northness rough-magnitude rough-scale spi tpi vrm ; do 
#     gdalbuildvrt -overwrite $MERIT/${TOPO}/all_${TOPO}_90M.vrt $MERIT/${TOPO}/${TOPO}_90M_???????.tif 
#     gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $MERIT/${TOPO}/all_${TOPO}_90M.vrt $ONCHO/input/geomorpho90m/${TOPO}.tif 
# done
 
# mv $ONCHO/input/geomorpho90m/elev-stdev.tif   $ONCHO/input/geomorpho90m/elev_stdev.tif 
# mv $ONCHO/input/geomorpho90m/aspect-sine.tif   $ONCHO/input/geomorpho90m/aspect_sine.tif 
# mv $ONCHO/input/geomorpho90m/aspect-cosine.tif   $ONCHO/input/geomorpho90m/aspect_cosine.tif 
# mv $ONCHO/input/geomorpho90m/rough-magnitude.tif   $ONCHO/input/geomorpho90m/rough_magnitude.tif 
# mv $ONCHO/input/geomorpho90m/rough-scale.tif   $ONCHO/input/geomorpho90m/rough_scale.tif 
# rm -f $ONCHO/input/geomorpho90m/aspect.tif $ONCHO/input/geomorpho90m/geom.tif $ONCHO/input/geomorpho90m/cti.tif $ONCHO/input/geomorpho90m/spi.tif


# pkreclass -c -9999 -r 0 -co COMPRESS=DEFLATE -co ZLEVEL=9  -i   $ONCHO/input/geomorpho90m/aspect_cosine.tif  -o  $ONCHO/input/geomorpho90m/aspect_cosine_tmp.tif 
# mv $ONCHO/input/geomorpho90m/aspect_cosine_tmp.tif   $ONCHO/input/geomorpho90m/aspect_cosine.tif 

# pkreclass -c -9999 -r 0 -co COMPRESS=DEFLATE -co ZLEVEL=9  -i   $ONCHO/input/geomorpho90m/aspect_sine.tif  -o  $ONCHO/input/geomorpho90m/aspect_sine_tmp.tif 
# mv $ONCHO/input/geomorpho90m/aspect_sine_tmp.tif   $ONCHO/input/geomorpho90m/aspect_sine.tif 
# gdal_translate -projwin 2 15 15 4  -co COMPRESS=DEFLATE -co ZLEVEL=9 /gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT/input_tif/all_tif.vrt  $ONCHO/input/geomorpho90m/elevation.tif

HYDRO=/gpfs/gibbs/pi/hydro/hydro/dataproces/MERIT_HYDRO/hydrography90m_v.1.0

for HYDROD in   flow.index  r.stream.distance  r.stream.slope  r.watershed  ; do 
    for DIR in $(ls $HYDRO/$HYDROD ) ; do 
     gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $HYDRO/$HYDROD/$DIR/$(basename $DIR _tiles20d).vrt /tmp/$(basename $DIR _tiles20d).tif 
     pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m /tmp/$(basename $DIR _tiles20d).tif -msknodata -9999  -nodata 0 -i  /tmp/$(basename $DIR _tiles20d).tif   -o $ONCHO/input/hydrography90m/$(basename $DIR _tiles20d).tif 
     rm /tmp/$(basename $DIR _tiles20d).tif 
     gdal_edit.py -a_nodata -9999 $ONCHO/input/hydrography90m/$(basename $DIR _tiles20d).tif 
    done 
done 

rm -f $ONCHO/input/hydrography90m/channel_*.tif $ONCHO/input/hydrography90m/order_*.tif $ONCHO/input/hydrography90m/{segment.tif,regional_unit.tif,basin.tif,depression.tif,direction.tif,sub_catchment.tif,outlet.tif}

pkreclass -c -9999999  -r -2147483648  -co COMPRESS=DEFLATE -co ZLEVEL=9  -i   $ONCHO/input/hydrography90m/accumulation.tif  -o  $ONCHO/input/hydrography90m/accumulation_tmp.tif 
mv $ONCHO/input/hydrography90m/accumulation_tmp.tif    $ONCHO/input/hydrography90m/accumulation.tif 

pkreclass -c -9999999  -r -2147483648 -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $ONCHO/input/hydrography90m/slope_curv_max_dw_cel.tif -o $ONCHO/input/hydrography90m/slope_curv_max_dw_cel_tmp.tif
mv $ONCHO/input/hydrography90m/slope_curv_max_dw_cel_tmp.tif  $ONCHO/input/hydrography90m/slope_curv_max_dw_cel.tif

pkreclass -c -9999999  -r -2147483648 -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $ONCHO/input/hydrography90m/slope_curv_min_dw_cel.tif -o $ONCHO/input/hydrography90m/slope_curv_min_dw_cel_tmp.tif
mv $ONCHO/input/hydrography90m/slope_curv_min_dw_cel_tmp.tif  $ONCHO/input/hydrography90m/slope_curv_min_dw_cel.tif

pkreclass -c -9999999  -r -2147483648 -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $ONCHO/input/hydrography90m/slope_grad_dw_cel.tif -o $ONCHO/input/hydrography90m/slope_grad_dw_cel_tmp.tif
mv $ONCHO/input/hydrography90m/slope_grad_dw_cel_tmp.tif  $ONCHO/input/hydrography90m/slope_grad_dw_cel.tif 

#### remove some tif


exit 

CHELSA=/gpfs/gibbs/pi/hydro/hydro/dataproces/CHELSA/climatologies/bio

for CHELSAD  in $CHELSA/CHELSA_bio*_1981-2010_V.2.1.tif  ; do
     filename=$(basename $CHELSAD _1981-2010_V.2.1.tif )
  #   gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $CHELSAD  $ONCHO/input/chelsa/${filename}.tif
     gdal_translate -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/input/chelsa/${filename}.tif  $ONCHO/input/chelsa/${filename}_r.tif
done

SOILT=/gpfs/gibbs/pi/hydro/hydro/dataproces/SOILTEMP/input

for SOIL  in $( ls  $SOILT/*.tif | grep -v msk )  ; do
filename=$(basename  $SOIL .tif)
# gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $SOIL  $ONCHO/input/soiltemp/${filename}.tif
# pkgetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -min  -9999999  -max 9999999 -data 1 -nodata 0   -i $ONCHO/input/soiltemp/${filename}.tif -o $ONCHO/input/soiltemp/${filename}_msk01.tif
# pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $ONCHO/input/soiltemp/${filename}_msk01.tif -msknodata 0 -nodata -9999 -i $ONCHO/input/soiltemp/${filename}.tif -o $ONCHO/input/soiltemp/${filename}_msk.tif
# gdal_fillnodata.py  -co COMPRESS=DEFLATE -co ZLEVEL=9 -nomask  -md 40  -si 1  $ONCHO/input/soiltemp/${filename}_msk.tif $ONCHO/input/soiltemp/${filename}.tif 
gdal_translate -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ONCHO/input/soiltemp/${filename}.tif  $ONCHO/input/soiltemp/${filename}_r.tif
done

gdalbuildvrt -overwrite -separate $ONCHO/input/soiltemp/all_tif_msk01.vrt   $ONCHO/input/soiltemp/*_msk01.tif
pkstatprofile -ot Byte  -of GTiff  -co COMPRESS=DEFLATE -co ZLEVEL=9 -f sum  -i $ONCHO/input/soiltemp/all_tif_msk01.vrt -o $ONCHO/input/soiltemp/soiltmp_msk_tmp.tif 
pkgetmask -ot Byte   -co COMPRESS=DEFLATE -co ZLEVEL=9  -min   0.5  -max 50  -data 1  -nodata 0   -i $ONCHO/input/soiltemp/soiltmp_msk_tmp.tif   -o $ONCHO/input/soiltemp/soiltmp_msk.tif 
rm $ONCHO/input/soiltemp/soiltmp_msk_tmp.tif $ONCHO/input/soiltemp/*_msk01.tif   $ONCHO/input/soiltemp/all_tif_msk01.vrt 

SOILGRIDS=/gpfs/gibbs/pi/hydro/hydro/dataproces/SOILGRIDS

# for SOIL in $SOILGRIDS/*_acc/*_WeigAver.vrt ; do
#     filename=$(basename  $SOIL  .vrt)
#     gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $SOIL  $ONCHO/input/soilgrids/${filename}_acc.tif
#     gdal_edit.py -a_nodata -9999  $ONCHO/input/soilgrids/soilgrids_msk.tif  # to fake the no-data
# done

for SOIL in $(ls $SOILGRIDS/*/*_WeigAver.tif |  grep -v  -e _acc -e out_  )  ; do
filename=$(basename  $SOIL  .tif)
# gdal_translate -projwin 2 15 15 4 -co COMPRESS=DEFLATE -co ZLEVEL=9   $SOIL  $ONCHO/input/soilgrids/${filename}_tmp.tif
# gdal_fillnodata.py  -co COMPRESS=DEFLATE -co ZLEVEL=9 -nomask  -md 40  -si 1  $ONCHO/input/soilgrids/${filename}_tmp.tif  $ONCHO/input/soilgrids/${filename}.tif
rm -f $ONCHO/input/soilgrids/${filename}_tmp.tif
gdal_translate -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ONCHO/input/soilgrids/${filename}.tif $ONCHO/input/soilgrids/${filename}_r.tif
done

pkgetmask -ot Byte -co COMPRESS=DEFLATE -co ZLEVEL=9 -min 65534 -max 65536 -data 0 -nodata 1 -i $ONCHO/input/soilgrids/AWCtS_WeigAver.tif -o $ONCHO/input/soilgrids/soilgrids_msk.tif
gdal_edit.py -a_nodata 0 $ONCHO/input/soilgrids/soilgrids_msk.tif


##### population density 

# alternative https://ghsl.jrc.ec.europa.eu/download.php?ds=pop 

# https://data.grid3.org/maps/GRID3::grid3-nigeria-gridded-population-estimates-version-2-0/about 
# cd /gpfs/loomis/project/sbsc/ga254/dataproces/ONCHO/input/population 
# wget https://wopr.worldpop.org/download/495
# unzip 495

POP=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/population
# https://ghsl.jrc.ec.europa.eu/download.php?ds=pop

wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R9_C19.zip
wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C19.zip
wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C20.zip
wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2022A/GHS_POP_E2020_GLOBE_R2022A_54009_100/V1-0/tiles/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C19.zip

unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C19.zip  
unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R8_C20.zip  
unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R9_C19.zip
unzip GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_R9_C20.zip


gdalbuildvrt -overwrite $POP/pop_100m.vrt  $POP/GHS_POP_E2020_GLOBE_R2022A_54009_100_V1_0_*.tif

gdalwarp -overwrite -ot Int16  -te $(getCorners4Gwarp $POP/../geomorpho90m/elevation.tif ) -t_srs EPSG:4326 -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9 $POP/pop_100m.vrt $POP/GHSpop_90m.tif
gdal_edit.py -a_nodata -9999 $POP/GHSpop_90m.tif   

# pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m /tmp/$(basename $DIR _tiles20d).tif -msknodata -9999  -nodata 0 -i  /tmp/$(basename $DIR _tiles20d).tif   -o $ONCHO/input/hydrography90m/$(basename $DIR _tiles20d).tif 


### livestock fatto il download a mano https://dataverse.harvard.edu/dataverse/glw_3


LS=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/livestock
for file in $LS/5_??_2010_Da.tif  ; do 
gdal_translate -ot Int32  -projwin  $(getCorners4Gtranslate  $LS/../geomorpho90m/elevation.tif )  -tr 0.00083333333333333 0.00083333333333333 -r bilinear  -co COMPRESS=DEFLATE -co ZLEVEL=9 $file  $LS/$(basename $file .tif )_90m.tif 
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


#  wget https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip 

ER=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/ecoregions
gdal_rasterize   -ot Byte   -te  $(getCorners4Gwarp   $ER/../geomorpho90m/elevation.tif )  -l "wwf_terr_ecos"  -a "G200_NUM"   -tr 0.0083333333333333 0.0083333333333333  -co COMPRESS=DEFLATE -co ZLEVEL=9  $ER/official/wwf_terr_ecos.shp  $ER/wwf_terr_ecosG200_NUM.tif 

