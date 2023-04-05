#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc09_downloadPoint_PredExtraction.sh.%j.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc09_downloadPoint_PredExtraction.sh.%j.err
#SBATCH --job-name=sc09_downloadPoint_PredExtraction.sh
#SBATCH --mem=40G

###### sbatch  /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc09_downloadPoint_PredExtraction.sh

export ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
cd $ONCHO/vector

module purge
module  load miniconda/4.10.3
conda activate gdx_env
source  /gpfs/gibbs/project/sbsc/ga254/conda_envs/gdx_env/lib/python3.1/venv/scripts/common/activate

python /vast/palmer/home.grace/ga254/scripts/ONCHO/gdx_download.py clean_data_nigeria_project.csv  b15f7670-64e1-4bfd-b6b1-5797e149514c /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/clean_data_nigeria_project.csv

conda  deactivate
module purge

source ~/bin/gdal3
source ~/bin/pktools
source ~/bin/grass78m

tr '\r' '\n' < $ONCHO/vector/clean_data_nigeria_project.csv | awk '!/^[[:space:]]*$/' | grep -v "94716"  > $ONCHO/vector/clean_data_nigeria_project_unix.csv

echo "x y pa" > $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt
awk -F "," '{ gsub("\"","") ; if (NR>1) { if ($5=="Yes")
                                                       print $4 , $3 , 1 
                                                       else 
              print $4 , $3 , 0 }}' $ONCHO/vector/clean_data_nigeria_project_unix.csv   >> $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt

                        # remove point in the sea > 4.2 
awk -F , '{if(NR>1 &&  $(NF-4) > 4.2 ) print $(NF-5)  , $(NF-4) , 1 }'   $ONCHO/vector/ecoregions_adaleke.csv | grep -v "0 0 1" | grep -v "3.05555555555556" >> $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt
awk -F , '{if(NR>10) print $(NF-5)  , $(NF-4) , 1 }'  $ONCHO/vector/Simulium_data_Akindele.csv  >>  $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt
awk -F , '{if ($8=="Present") {PA=1} else {PA=0}  ;  if(NR>1) print $5  , $4 , PA  }'  $ONCHO/vector/Blackfly_data_arimoro.csv >> $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt

awk  '{ if (NR>1) print $1 , $2 }'      $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt   >  $ONCHO/vector/NigeriaHabitatSites_x_y.txt

awk  '{ if (NR>1) print $1 , $2 , $3 }'  $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt   >  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_noheader.txt
rm -fr  $ONCHO/vector/NigeriaHabitatSites_x_y_pa.gpkg
pkascii2ogr   -n "PA"  -a_srs epsg:4326 -f GPKG  -i $ONCHO/vector/NigeriaHabitatSites_x_y_pa_noheader.txt -o $ONCHO/vector/NigeriaHabitatSites_x_y_pa.gpkg
gdal_rasterize  -co COMPRESS=DEFLATE -co ZLEVEL=9  -ot Byte -a_nodata 255   -a "PA" -tr 0.000833333333333 0.000833333333333  -te 2 4 15 15  -a_srs EPSG:4326 $ONCHO/vector/NigeriaHabitatSites_x_y_pa.gpkg  $ONCHO/vector/NigeriaHabitatSites_x_y_pa.tif


echo "ncols        15600"                  >  $ONCHO/input/id/raster_ID.asc 
echo "nrows        13200"                  >> $ONCHO/input/id/raster_ID.asc 
echo "xllcorner    2"                      >> $ONCHO/input/id/raster_ID.asc 
echo "yllcorner    4"                      >> $ONCHO/input/id/raster_ID.asc 
echo "cellsize     0.000833333333333"      >> $ONCHO/input/id/raster_ID.asc 

awk ' BEGIN {  
for (row=1 ; row<=13200 ; row++)  { 
     for (col=1 ; col<=15600 ; col++) { 
         printf ("%i " ,  col+(row-1)*15600  ) } ; printf ("\n")  }}' >> $ONCHO/input/id/raster_ID.asc 
 
# transform the created arcinfo ascii grid in a tif.
 
gdal_translate   -a_srs EPSG:4326 -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/input/id/raster_ID.asc  $ONCHO/input/id/raster_ID.tif 
rm   $ONCHO/input/id/raster_ID.asc

gdallocationinfo -geoloc -wgs84 -valonly  $ONCHO/input/id/raster_ID.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y.txt  > $ONCHO/vector/NigeriaHabitatSites_x_y_ID.txt 
paste -d " "  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_noheader.txt   $ONCHO/vector/NigeriaHabitatSites_x_y_ID.txt | sort -k 4,4   > $ONCHO/vector/NigeriaHabitatSites_x_y_pa_ID.txt 
bash   /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/average.sh $ONCHO/vector/NigeriaHabitatSites_x_y_pa_ID.txt  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_ID_uniq.txt  <<EOF 
n
4
8
EOF

awk '{if($3==0) {PA=0} else {PA=1}  ;   print $1 , $2 , PA  }'   $ONCHO/vector/NigeriaHabitatSites_x_y_pa_ID_uniq.txt > $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_noheader.txt

rm $ONCHO/vector/NigeriaHabitatSites_x_y_ID.txt $ONCHO/vector/NigeriaHabitatSites_x_y_pa_ID_uniq.txt

awk  '{ print $1 , $2 }'      $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_noheader.txt  >  $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt

echo "x y pa" >  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header.txt
awk  '{ print $1 , $2 , $3 }' $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_noheader.txt  >>  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header.txt

echo "x y" >  $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_header.txt
awk  '{ print $1 , $2}' $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_noheader.txt  >>  $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_header.txt


#### crate 0 proximit and 1 proximity

pkascii2ogr -n "PA" -a_srs epsg:4326 -f GPKG -i $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_noheader.txt -o $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq.gpkg
 
# grass78  -f -text --tmp-location  -c $ONCHO/input/geomorpho90m/elevation.tif     <<'EOF'

# for radius in 0.2 0.4 0.5 0.6 0.8 1 ; do 
# r.external  input=$ONCHO/input/geomorpho90m/elevation.tif output=elv --overwrite
# g.region res=0:00:30
# v.in.ogr  input=$ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq.gpkg output=pa   where="PA = 0" --o
# v.kernel input=pa output=kernel   radius=$radius   multiplier=1 --o
# g.region res=0:00:03
# r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=-9999 type=Float32 format=GTiff input=kernel output=$ONCHO/input/pa_kernel/kernel0_tmp.tif
# pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $ONCHO/input/geomorpho90m/elevation.tif  -msknodata -9999 -nodata -9999 -i $ONCHO/input/pa_kernel/kernel0_tmp.tif -o $ONCHO/input/pa_kernel/kernel0_R$radius.tif
# g.region res=0:00:30
# v.in.ogr  input=$ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq.gpkg output=pa   where="PA = 1" --o 
# v.kernel input=pa output=kernel   radius=$radius  multiplier=1 --o
# g.region res=0:00:03
# r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=-9999 type=Float32 format=GTiff input=kernel output=$ONCHO/input/pa_kernel/kernel1_tmp.tif
# pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $ONCHO/input/geomorpho90m/elevation.tif  -msknodata -9999 -nodata -9999 -i $ONCHO/input/pa_kernel/kernel1_tmp.tif -o $ONCHO/input/pa_kernel/kernel1_R$radius.tif
# done 
# EOF

rm $ONCHO/input/pa_kernel/kernel0_tmp.tif $ONCHO/input/pa_kernel/kernel1_tmp.tif

### start the predictors  extraction 

rm -f $ONCHO/vector/pred_*.txt 

for var in geomorpho90m hydrography90m ; do 

gdalbuildvrt -separate  -overwrite $ONCHO/input/$var/all_tif.vrt $(ls $ONCHO/input/$var/*.tif | grep -v -e _acc -e _msk)
BB=$(ls $ONCHO/input/$var/*.tif | grep -v -e _acc -e _msk  | wc -l )

for var1 in $( ls $ONCHO/input/$var/*.tif  | grep -v -e _acc -e _msk )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}.txt
echo "" >> $ONCHO/vector/pred_${var}.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/$var/all_tif.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt   |   awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}.txt

done   

for var in chelsa soilgrids soiltemp ; do 

gdalbuildvrt -separate  -overwrite $ONCHO/input/$var/all_tif.vrt $(ls $ONCHO/input/$var/*_r.tif | grep -v -e _acc -e _msk)
BB=$(ls $ONCHO/input/$var/*_r.tif | grep -v -e _acc -e _msk  | wc -l )

for var1 in $( ls $ONCHO/input/$var/*_r.tif  | grep -v -e _acc -e _msk )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}.txt
echo "" >> $ONCHO/vector/pred_${var}.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/$var/all_tif.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt | awk -v BB=$BB   |   'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}.txt

done #### close the for var in chelsa soilgrids soiltemp 

var=soilgrids
gdalbuildvrt -separate  -overwrite $ONCHO/input/$var/all_tif_acc.vrt $(ls $ONCHO/input/$var/*.tif | grep -e _acc | grep -v _r.tif  )
BB=$(ls $ONCHO/input/$var/*.tif | grep -e _acc | grep -v _r.tif | wc -l )

for var1 in $( ls $ONCHO/input/$var/*.tif  | grep -e _acc | grep -v _r.tif )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}_acc.txt
echo "" >> $ONCHO/vector/pred_${var}_acc.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/$var/all_tif_acc.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}_acc.txt

#####vegetation index 

echo "ndti" > $ONCHO/vector/pred_ndti.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndti.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_ndti.txt
echo "ndvi" > $ONCHO/vector/pred_ndvi.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndvi.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_ndvi.txt
echo "ndwi" > $ONCHO/vector/pred_ndwi.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndwi.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_ndwi.txt

#### livestock 

gdalbuildvrt -separate  -overwrite $ONCHO/input/livestock/all_tif.vrt  $ONCHO/input/livestock/LS_??.tif 
BB=$(ls  $ONCHO/input/livestock/LS_??.tif  | wc -l )

for var in  $ONCHO/input/livestock/LS_??.tif; do echo -n $(basename $var .tif)" " ; done > $ONCHO/vector/pred_LS.txt
echo "" >> $ONCHO/vector/pred_LS.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/livestock/all_tif.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_LS.txt


##### population 

echo "GHSpop_90m" > $ONCHO/vector/pred_pop.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/population/GHSpop_90m.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_pop.txt

##### water 
echo "occurence" > $ONCHO/vector/pred_occ.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/water/occurence.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_occ.txt

echo "occurence_proximity" > $ONCHO/vector/pred_occ_prox.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/water/occurence.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_occ_prox.txt

echo "flow_mean" > $ONCHO/vector/pred_flow.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/water/flow_mean.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_flow.txt 

####  landcover 

echo "LC2021" > $ONCHO/vector/pred_landcover.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/landcover/LC2021.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_landcover.txt

##### ecorigion 

echo "ER2017" > $ONCHO/vector/pred_ecoreg.txt
gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input/ecoregions/ERcase 
2017.tif   < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_ecoreg.txt

# #### kernel  
for radius in 3 ; do 
echo "kernel1_R$radius" > $ONCHO/vector/pred_kernel1_R$radius.txt
gdallocationinfo -geoloc -wgs84 -valonly    $ONCHO/input/pa_kernel/kernel1_R$radius.tif   < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_kernel1_R$radius.txt

echo "kernel0_R$radius" > $ONCHO/vector/pred_kernel0_R$radius.txt
gdallocationinfo -geoloc -wgs84 -valonly    $ONCHO/input/pa_kernel/kernel0_R$radius.tif   < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt  >> $ONCHO/vector/pred_kernel0_R$radius.txt
done 


#### lat long  not ingested as predictor because already in the x y pa table 

#### merging predictors 

paste -d " " $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header.txt $ONCHO/vector/pred_*.txt |  sed  's,  , ,g' > $ONCHO/vector/x_y_pa_predictors.txt
awk '{ print $1="", $2="", $0  }' $ONCHO/vector/x_y_pa_predictors.txt |  sed  's,    ,,g'  > $ONCHO/vector/x_y_pa_predictors4R.txt


