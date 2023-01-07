#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc10_downloadPoint_RFmodel_mlr3spatial.sh.%A_%a.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc10_downloadPoint_RFmodel_mlr3spatial.sh.%A_%a.err
#SBATCH --job-name=sc10_downloadPoint_RFmodel_mlr3spatial.sh
#SBATCH --mem=100G
#SBATCH --array=1
## array=1-41

##### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc10_downloadPoint_RFmodel_rfClasic.sh

#  x 2 15 
#  y 4 15                                                                                                              remove the first row that includes only see area
# for x in $(seq 2 2 14) ; do for y in $(seq 4 2 14 ) ; do echo $x $(expr $x + 2 ) $y $(expr $y + 2 ) ; done ; done  | awk '{ if (NR>1) print  }'  >   $ONCHO/vector/tile_list.txt  

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
cd $ONCHO/vector

# SLURM_ARRAY_TASK_ID=1

geo_string=$(head  -n  $SLURM_ARRAY_TASK_ID $ONCHO/vector/tile_list.txt   | tail  -1 )
export xmin=$( echo $geo_string | awk '{  print $1 }' ) 
export xmax=$( echo $geo_string | awk '{  print $2 }' ) 
export ymin=$( echo $geo_string | awk '{  print $3 }' ) 
export ymax=$( echo $geo_string | awk '{  print $4 }' ) 

echo geo_string  =  $xmin  $xmax $ymin $ymax
echo prediction_${xmin}_${ymin}.tif


if [ $SLURM_ARRAY_TASK_ID -eq 1  ] ; then

rm -f $ONCHO/vector/data*.RData
rm -f $ONCHO/vector/allVar.mod.rf.txt
rm -f $ONCHO/vector/importance_allVar.txt
rm -f $ONCHO/prediction/prediction_*_*.tif
rm -f $ONCHO/prediction/*.tif.aux.xml

# module purge
# module load miniconda/4.10.3
# conda activate gdx_env
# source  /gpfs/gibbs/project/sbsc/ga254/conda_envs/gdx_env/lib/python3.1/venv/scripts/common/activate

# python /vast/palmer/home.grace/ga254/scripts/ONCHO/gdx_download.py clean_data_nigeria_project.csv  b15f7670-64e1-4bfd-b6b1-5797e149514c /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/clean_data_nigeria_project.csv

# conda  deactivate
# module purge
source ~/bin/gdal3
source ~/bin/pktools

# tr '\r' '\n' < $ONCHO/vector/clean_data_nigeria_project.csv | awk '!/^[[:space:]]*$/' | grep -v "94716"  > $ONCHO/vector/clean_data_nigeria_project_unix.csv

# echo "x y pa" > $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt
# awk -F "," '{ gsub("\"","") ; if (NR>1) { if ($5=="Yes")
#                                                        print $4 , $3 , 1 
#                                                        else 
#                                                        print $4 , $3 , 0 }}' $ONCHO/vector/clean_data_nigeria_project_unix.csv   >> $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt
# awk  '{ if (NR>1) print $1 , $2 }'      $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt   >  $ONCHO/vector/NigeriaHabitatSites_x_y.txt

# awk  '{ if (NR>1) print $1 , $2 , $3 }'      $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt   >  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_noheader.txt
# rm -fr  $ONCHO/vector/NigeriaHabitatSites_x_y_pa.gpkg
# pkascii2ogr   -n "PA"  -a_srs epsg:4326 -f GPKG  -i $ONCHO/vector/NigeriaHabitatSites_x_y_pa_noheader.txt -o $ONCHO/vector/NigeriaHabitatSites_x_y_pa.gpkg
# gdal_rasterize  -co COMPRESS=DEFLATE -co ZLEVEL=9  -ot Byte -a_nodata 255   -a "PA" -tr 0.000833333333333 0.000833333333333  -te 2 4 15 15  -a_srs EPSG:4326 $ONCHO/vector/NigeriaHabitatSites_x_y_pa.gpkg  $ONCHO/vector/NigeriaHabitatSites_x_y_pa.tif

# rm -f $ONCHO/vector/pred_*.txt 

# for var in geomorpho90m hydrography90m ; do 

# gdalbuildvrt -separate  -overwrite $ONCHO/inputND/$var/all_tif.vrt $(ls $ONCHO/inputND/$var/*.tif | grep -v -e _acc -e _msk)
# BB=$(ls $ONCHO/inputND/$var/*.tif | grep -v -e _acc -e _msk  | wc -l )

# for var1 in $( ls $ONCHO/inputND/$var/*.tif  | grep -v -e _acc -e _msk )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}.txt
# echo "" >> $ONCHO/vector/pred_${var}.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/inputND/$var/all_tif.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}.txt

# done 

# for var in chelsa soilgrids soiltemp ; do 

# gdalbuildvrt -separate  -overwrite $ONCHO/inputND/$var/all_tif.vrt $(ls $ONCHO/inputND/$var/*_r.tif | grep -v -e _acc -e _msk)
# BB=$(ls $ONCHO/inputND/$var/*_r.tif | grep -v -e _acc -e _msk  | wc -l )

# for var1 in $( ls $ONCHO/inputND/$var/*_r.tif  | grep -v -e _acc -e _msk )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}.txt
# echo "" >> $ONCHO/vector/pred_${var}.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/inputND/$var/all_tif.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}.txt

# done ### close the for var in chelsa soilgrids soiltemp 

# var=soilgrids
# gdalbuildvrt -separate  -overwrite $ONCHO/inputND/$var/all_tif_acc.vrt $(ls $ONCHO/inputND/$var/*.tif | grep -e _acc | grep -v _r.tif  )
# BB=$(ls $ONCHO/inputND/$var/*.tif | grep -e _acc | grep -v _r.tif | wc -l )

# for var1 in $( ls $ONCHO/inputND/$var/*.tif  | grep -e _acc | grep -v _r.tif )  ; do echo -n $(basename  $var1 .tif)" " ; done > $ONCHO/vector/pred_${var}_acc.txt
# echo "" >> $ONCHO/vector/pred_${var}_acc.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/inputND/$var/all_tif_acc.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_${var}_acc.txt

# ### vegetation index 

# echo "ndti" > $ONCHO/vector/pred_ndti.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndti.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y.txt  >> $ONCHO/vector/pred_ndti.txt
# echo "ndvi" > $ONCHO/vector/pred_ndvi.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndvi.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y.txt  >> $ONCHO/vector/pred_ndvi.txt
# echo "ndwi" > $ONCHO/vector/pred_ndwi.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/input_orig/sentinel/tif_wgs84_90m_index/ndwi.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y.txt  >> $ONCHO/vector/pred_ndwi.txt

# ### livestock 

# gdalbuildvrt -separate  -overwrite $ONCHO/inputND/livestock/all_tif.vrt  $ONCHO/inputND/livestock/LS_??.tif 
# BB=$(ls  $ONCHO/inputND/livestock/LS_??.tif  | wc -l )

# for var in  $ONCHO/inputND/livestock/LS_??.tif; do echo -n $(basename $var .tif)" " ; done > $ONCHO/vector/pred_LS.txt
# echo "" >> $ONCHO/vector/pred_LS.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/inputND/livestock/all_tif.vrt < $ONCHO/vector/NigeriaHabitatSites_x_y.txt | awk -v BB=$BB 'ORS=NR%BB?FS:RS' >> $ONCHO/vector/pred_LS.txt


# ### population 

# echo "GHSpop_90m" > $ONCHO/vector/pred_pop.txt
# gdallocationinfo -geoloc -wgs84 -valonly $ONCHO/inputND/population/GHSpop_90m.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y.txt  >> $ONCHO/vector/pred_pop.txt

# paste -d " " $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt $ONCHO/vector/pred_*.txt |  sed  's,  , ,g' > $ONCHO/vector/x_y_pa_predictors.txt
# awk '{ print $1="", $2="", $0  }' $ONCHO/vector/x_y_pa_predictors.txt |  sed  's,    ,,g'  > $ONCHO/vector/x_y_pa_predictors4R.txt

# #### fi ### remove with the exit 
# #### exit 

echo "#######################################################"
echo "############### FIRST RF #############################"
echo "#######################################################"

module load R/4.1.0-foss-2020b

# # # see http://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/CompareRandomForestPackages.html
# # # R  --vanilla --no-readline   -q  <<'EOF'  this is not working with ranger 

# first RF for sorting the most important variables

Rscript --vanilla  -e '
library(psych)
library(randomForest)
set.seed(1)

table = read.table("x_y_pa_predictors4R.txt", header = TRUE, sep = " ")
table$geom = as.factor(table$geom)
table$pa = as.factor(table$pa)

des.table = describe(table)

write.table(des.table, "stat_allVar.txt", quote = FALSE  )

# mod.rfP = ranger( pa ~ . , data=table ,   importance=TRUE, proximity=TRUE)
# print(mod.rfP)

mod.rfR = randomForest( pa ~ . , data=table ,   importance=TRUE, proximity=TRUE)
print(mod.rfR)

# impP=as.data.frame(importance(mod.rfP))
impR=as.data.frame(importance(mod.rfR))
save.image("data0.RData")

# impP.s = impP[order(impP$"importance(mod.rfP)",decreasing=TRUE), , drop = FALSE]
impR.s = impR[order(impR$MeanDecreaseAccuracy,decreasing=TRUE), , drop = FALSE] 

# impP.s
impR.s

write.table(impR.s, "importance_allVar.txt", quote = FALSE  )
s.mod.rfR = capture.output(mod.rfR)
write.table(s.mod.rfR, "allVar.mod.rf.txt", quote = FALSE , row.names = FALSE )

save.image("data0.RData")
'
 

IMPN=40
rm -f $ONCHO/vector/x_y_pa_*_tmp.txt
for COLNAME in  $(awk -v IMPN=$IMPN  '{if (NR>1 && NR<=IMPN) print $1  }' $ONCHO/vector/importance_allVar.txt) ; do
awk  -v COLNAME=$COLNAME ' { if (NR==1){ for (col=1;col<=NF;col++) { if ($col==COLNAME) {colprint=col; print $colprint}}} else {print $colprint }}' $ONCHO/vector/x_y_pa_predictors4R.txt  > $ONCHO/vector/x_y_pa_${COLNAME}_tmp.txt
done 

for COLNAME in  $(awk -v IMPN=$IMPN '{if (NR>1 && NR<=IMPN) print $1  }' $ONCHO/vector/importance_allVar.txt) ; do
    echo $COLNAME $(gdalinfo $ONCHO/inputND/*/$COLNAME.tif | grep "NoData"  | awk -F =  '{ print $2  }' )
done > $ONCHO/vector/NoData_predictors4R_select.txt

paste -d " " <(awk '{print $3}' $ONCHO/vector/NigeriaHabitatSites_x_y_pa.txt) $ONCHO/vector/x_y_pa_*_tmp.txt > $ONCHO/vector/x_y_pa_predictors4R_select.txt
rm -f $ONCHO/vector/x_y_pa_*_tmp.txt

#### module load GCCcore/11.2.0 to install rgdal
#### options(menu.graphics=FALSE) for no gui 
echo "#######################################################"
echo "############### SECOND RF #############################"
echo "#######################################################"

#### https://mlr3spatial.mlr-org.com/articles/mlr3spatial.html

# training RF and save the model for make the prediction later on 

Rscript --vanilla  -e '
library("randomForest")
library("stars")
library("terra")
library("future")
set.seed(1)

table = read.table("x_y_pa_predictors4R_select.txt", header = TRUE, sep = " ")
table$pa = as.factor(table$pa)

mod.rfR = randomForest( pa ~ . , data=table ,   importance=TRUE, proximity=TRUE)
print(mod.rfR)

# impP=as.data.frame(importance(mod.rfP))
impR=as.data.frame(importance(mod.rfR))
save.image("data0.RData")

# impP.s = impP[order(impP$"importance(mod.rfP)",decreasing=TRUE), , drop = FALSE]
impR.s = impR[order(impR$MeanDecreaseAccuracy,decreasing=TRUE), , drop = FALSE] 

# impP.s
impR.s


OOB.votes <- predict (mod.rfR,table,type="prob");
OOB.pred <- OOB.votes[,2];

pred.obj <- prediction (OOB.pred,table$pa);

RP.perf <- performance(pred.obj, "rec","prec");
plot (RP.perf);





write.table(impR.s, "importance_selVar.txt", quote = FALSE  )
s.mod.rfR = capture.output(mod.rfR)
write.table(s.mod.rfR, "sellVar.mod.rf.txt", quote = FALSE , row.names = FALSE )

save.image("data1.RData")
'

else 
  sleep 600 ### 3000  in case of runnging th gdallocation  
fi ### close the first array loop 

module purge
module load netCDF/4.7.4-gompi-2020b
source ~/bin/gdal3
module load R/4.1.0-foss-2020b

echo geo_string  =  $xmin  $xmax $ymin $ymax

### make the RF prediction 

Rscript --vanilla  -e   '
library("randomForest")
#library("Rcpp")
library("raster")
# library("future")

set.seed(1)

xmin <- as.numeric(Sys.getenv("xmin"))
xmax <- as.numeric(Sys.getenv("xmax"))
ymin <- as.numeric(Sys.getenv("ymin"))
ymax <- as.numeric(Sys.getenv("ymax"))

xmin
xmax
ymin
ymax

load("data1.RData")

print ("start the prediction")
# bb = st_bbox ( c(xmin = xmin , xmax =  xmax , ymin =  ymin, ymax =  ymax  )  , crs = 4326 ) 

for (var in names(table)[-1] ) {
print(var)
raster  = raster::raster (Sys.glob(paste0("../inputND/*/",var,".tif")  ))
# raster_crop = raster[bb] 

print(raster)
assign(paste0(var) , raster)
}
rm (raster)
gc() ; gc()

print ("make stack layer")
stack = get(names(table)[2] )
stack

save.image(paste0("data2_",xmin,"_",ymin,".RData"))

for (var in names(table)[-2:-1] ) {
print(var)
stack =  stack(stack,get(var))
}

print ( c(xmin, xmax  , ymin , ymax) )

save.image(paste0("data2_",xmin,"_",ymin,".RData"))

extent  <- extent( xmin , xmax , ymin , ymax )
extent  
env=crop(stack,extent)

print ("stack  info")
env
print ("env  str")
str(env)
 
# rm (stack)
gc() ;  gc() ;  

save.image(paste0("data2_",xmin,"_",ymin,".RData"))

env_predR = predict(env , mod.rfR, type="response" )
raster::writeRaster (env_predR, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction/predictionR_",xmin,"_",ymin,".tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Byte")

save.image(paste0("data3_",xmin,"_",ymin,".RData"))

'
rm -f  $ONCHO/prediction/predictionP_${xmin}_${ymin}.tif.aux.xml  $ONCHO/prediction/predictionR_${xmin}_${ymin}.tif.aux.xml

if [ $SLURM_ARRAY_TASK_ID -eq 41  ] ; then
sleep 2000
module purge
source ~/bin/gdal3
source ~/bin/pktools

rm -f $ONCHO/prediction/prediction_all.* $ONCHO/prediction/prediction_all_1km.*  $ONCHO/prediction/prediction?_all.* $ONCHO/prediction/prediction?_all_1km.*

gdalbuildvrt  $ONCHO/prediction/predictionP_all.vrt $ONCHO/prediction/predictionP_[0-9]*_[0-9]*.tif
gdal_translate  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/prediction/predictionP_all.vrt  $ONCHO/prediction/predictionP_all.tif 

gdalbuildvrt  $ONCHO/prediction/predictionR_all.vrt $ONCHO/prediction/predictionR_[0-9]*_[0-9]*.tif
gdal_translate  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/prediction/predictionR_all.vrt  $ONCHO/prediction/predictionR_all.tif 

pksetmask -m $ONCHO/input/hydrography90m/accumulation.tif -co COMPRESS=DEFLATE -co ZLEVEL=9 -msknodata -2147483648 -nodata -9999 -i $ONCHO/prediction/predictionP_all.tif -o $ONCHO/prediction/predictionP_all_msk.tif 

pksetmask -m $ONCHO/input/hydrography90m/accumulation.tif -co COMPRESS=DEFLATE -co ZLEVEL=9 -msknodata -2147483648 -nodata 255  -i $ONCHO/prediction/predictionR_all.tif -o $ONCHO/prediction/predictionR_all_msk.tif 

cd $ONCHO/prediction/
rm -f gdaltindex $ONCHO/prediction/all_tif_shp.*
gdaltindex $ONCHO/prediction/all_tif_shp.shp  predictionR_[0-9]*_[0-9]*.tif

gdal_translate -tr 0.00833333333333 0.00833333333333 -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 $ONCHO/prediction/predictionP_all.tif $ONCHO/prediction/predictionP_all_1km.tif
pksetmask -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9   -msknodata -9999 -nodata  -9999 -i $ONCHO/prediction/predictionP_all_1km.tif -o   $ONCHO/prediction/predictionP_all_1km_msk.tif
# gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9 -scale 0.344 0.966 0.01 1 $ONCHO/prediction/prediction_all_1km_msk.tif $ONCHO/prediction/prediction_all_1km_msk_s.tif

gdal_translate -tr 0.00833333333333 0.00833333333333 -r mode -co COMPRESS=DEFLATE -co ZLEVEL=9 $ONCHO/prediction/predictionR_all.tif $ONCHO/prediction/predictionR_all_1km.tif
pksetmask -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9   -msknodata -9999 -nodata 255  -i $ONCHO/prediction/predictionR_all_1km.tif -o   $ONCHO/prediction/predictionR_all_1km_msk.tif
# gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9 -scale 0.344 0.966 0.01 1 $ONCHO/prediction/prediction_all_1km_msk.tif $ONCHO/prediction/prediction_all_1km_msk_s.tif


fi
