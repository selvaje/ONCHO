#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc10_downloadPoint_RFmodel_mlr3terra.sh.%A_%a.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc10_downloadPoint_RFmodel_mlr3terra.sh.%A_%a.err
#SBATCH --job-name=sc10_downloadPoint_RFmodel_mlr3terra.sh
#SBATCH --mem=40G
#SBATCH --array=1-41
## array=1-41

#####  for seed in $(seq 1 20);  do sbatch  --export=seed=$seed /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc10_downloadPoint_RFmodel_mlr3terra.sh ; done 

#  x 2 15 
#  y 4 15                                                                                                              remove the first row that includes only see area
# for x in $(seq 2 2 14) ; do for y in $(seq 4 2 14 ) ; do echo $x $(expr $x + 2 ) $y $(expr $y + 2 ) ; done ; done  | awk '{ if (NR>1) print  }'  >   $ONCHO/vector/tile_list.txt  

module load StdEnv

echo "seed = $seed "

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

echo $seed 

if [ $SLURM_ARRAY_TASK_ID -eq 1  ] ; then

rm -f $ONCHO/vector_seed$seed/data*.RData
rm -f $ONCHO/vector_seed$seed/allVar.mod.rf.txt
rm -f $ONCHO/vector_seed$seed/importance_allVar.txt
rm -f $ONCHO/prediction_seed$seed/prediction_*_*.tif
rm -f $ONCHO/prediction_seed$seed/*.tif.aux.xml

source ~/bin/gdal3
source ~/bin/pktools

echo "#######################################################"
echo "############### FIRST RF #############################"
echo "#######################################################"

module --ignore_cache load R/4.1.0-foss-2020b

# # # see http://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/CompareRandomForestPackages.html
# # # R  --vanilla --no-readline   -q  <<'EOF'  this is not working with ranger 

# first RF for sorting the most important variables

Rscript --vanilla  -e '
library(ranger)
library(psych)

seed <- as.numeric(Sys.getenv("seed"))
set.seed(seed)

table = read.table("x_y_pa_predictors4R.txt", header = TRUE, sep = " ")
des.table = describe(table)

write.table(des.table, "stat_allVar.txt", quote = FALSE  )

mod.rfP = ranger( pa ~ . , table ,   probability = TRUE  ,  classification=TRUE ,   importance="permutation")
print(mod.rfP)
#print(mod.rfP$oob_error)

mod.rfR = ranger( pa ~ . , table , probability = FALSE  ,  classification=TRUE ,   importance="permutation")
print(mod.rfR)
# print(mod.rfR$oob_error)


impP=as.data.frame(importance(mod.rfP))
impR=as.data.frame(importance(mod.rfR))
save.image(paste0("../vector_seed",seed,"/data0.RData"))

impP.s = impP[order(impP$"importance(mod.rfP)",decreasing=TRUE), , drop = FALSE]
impR.s = impR[order(impR$"importance(mod.rfR)",decreasing=TRUE), , drop = FALSE]

impP.s
impR.s

write.table(impP.s, paste0("../vector_seed",seed,"/importanceP_allVar.txt"), quote = FALSE  )
s.mod.rfP = capture.output(mod.rfP)
write.table(s.mod.rfP, paste0("../vector_seed",seed,"/allVarP.mod.rf.txt"), quote = FALSE , row.names = FALSE )

write.table(impR.s, paste0("../vector_seed",seed,"/importanceR_allVar.txt"), quote = FALSE  )
s.mod.rfR = capture.output(mod.rfR)
write.table(s.mod.rfR, paste0("../vector_seed",seed,"/allVarR.mod.rf.txt"), quote = FALSE , row.names = FALSE )

save.image(paste0("../vector_seed",seed,"/data0.RData"))
'

# 

IMPN=40
rm -f $ONCHO/vector_seed$seed/x_y_pa_*_tmp.txt
for COLNAME in  $(awk -v IMPN=$IMPN  '{if (NR>1 && NR<=IMPN) print $1  }' $ONCHO/vector_seed$seed/importanceR_allVar.txt) ; do
awk  -v COLNAME=$COLNAME ' { if (NR==1){ for (col=1;col<=NF;col++) { if ($col==COLNAME) {colprint=col; print $colprint}}} else {print $colprint }}' $ONCHO/vector/x_y_pa_predictors4R.txt  > $ONCHO/vector_seed$seed/x_y_pa_${COLNAME}_tmp.txt
done 

for COLNAME in  $(awk -v IMPN=$IMPN '{if (NR>1 && NR<=IMPN) print $1  }' $ONCHO/vector_seed$seed/importanceR_allVar.txt) ; do
    echo $COLNAME $(gdalinfo $ONCHO/input/*/$COLNAME.tif | grep "NoData"  | awk -F =  '{ print $2  }' )
done > $ONCHO/vector_seed$seed/predictors4R_select.txt

paste -d " " <(awk '{print $3}' $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header.txt ) $ONCHO/vector_seed$seed/x_y_pa_*_tmp.txt > $ONCHO/vector_seed$seed/x_y_pa_predictors4R_select.txt
rm -f $ONCHO/vector_seed$seed/x_y_pa_*_tmp.txt


#### module load GCCcore/11.2.0 to install rgdal
#### options(menu.graphics=FALSE) for no gui 
echo "#######################################################"
echo "############### SECOND RF #############################"
echo "#######################################################"

#### https://mlr3spatial.mlr-org.com/articles/mlr3spatial.html

# training RF and save the model for make the prediction later on 

Rscript --vanilla  -e '
library("rlang")
library(mlr3spatiotempcv)  # spatio-temporal resampling 
library(mlr3tuning)        # hyperparameter tuning package
library("mlr3learners")
library("ranger")
library("stars")
library("terra")
library("future")
seed <- as.numeric(Sys.getenv("seed"))
set.seed(seed)

table = read.table(paste0("../vector_seed",seed,"/x_y_pa_predictors4R_select.txt"), header = TRUE, sep = " ")
table$pa = as.factor(table$pa)  # usefull for il backend
summary(table)

backend = as_data_backend(table)    # this is just table for the learner

task = as_task_classif(backend, target = "pa")
print(task)

learnerP = lrn("classif.ranger" , predict_type = "prob" , importance = "permutation"  )    ### https://mlr3extralearners.mlr-org.com/articles/learners/list_learners.html 
learnerP$parallel_predict = TRUE  
print(learnerP)
learnerP$train(task)  # usefull to obtain the $model
learnerP$model

learnerR = lrn("classif.ranger" , predict_type = "response" , importance = "permutation"  )    ### https://mlr3extralearners.mlr-org.com/articles/learners/list_learners.html 
learnerR$parallel_predict = TRUE  
print(learnerR)
learnerR$train(task)
learnerR$model

conf.matrix = learnerR$model$confusion.matrix
pred.error =  learnerR$model$prediction.error 

write.table(pred.error , paste0("../vector_seed",seed,"/pred_error.txt"))
write.table(conf.matrix, paste0("../vector_seed",seed,"/conf_matrix.txt"))

write.table(capture.output(learnerR$model), paste0("../vector_seed",seed,"/selVarR.mod.rf.txt"), quote = FALSE , row.names = FALSE )
write.table(capture.output(learnerP$model), paste0("../vector_seed",seed,"/selVarP.mod.rf.txt"), quote = FALSE , row.names = FALSE )

save.image(paste0("../vector_seed",seed,"/data1.RData"))
'
else 
sleep 1000 
fi ### close the first array loop 

module purge
module --ignore_cache load netCDF/4.7.4-gompi-2020b
source ~/bin/gdal3
module --ignore_cache load R/4.1.0-foss-2020b

echo geo_string  =  $xmin  $xmax $ymin $ymax

### make the RF prediction 

Rscript --vanilla  -e   '
library("rlang")
library("mlr3")
library("mlr3spatial")
library("mlr3learners")
library("ranger")
library("stars")
library("terra")
library("future")

seed <- as.numeric(Sys.getenv("seed"))
set.seed(seed)

xmin <- as.numeric(Sys.getenv("xmin"))
xmax <- as.numeric(Sys.getenv("xmax"))
ymin <- as.numeric(Sys.getenv("ymin"))
ymax <- as.numeric(Sys.getenv("ymax"))

xmin
xmax
ymin
ymax

load(paste0("../vector_seed",seed,"/data1.RData"))

print ("define response variable")
#### training ml 
learnerP$train(task)    #### define response variable 
learnerR$train(task)    #### define response variable 
print(learnerP)
print(learnerP)

print ("start the prediction")
# bb = st_bbox ( c(xmin = xmin , xmax =  xmax , ymin =  ymin, ymax =  ymax  )  , crs = 4326 ) 

for (var in names(table)[-1] ) {
print(var)
raster  = terra::rast (Sys.glob(paste0("../input/*/",var,".tif")  ) )

print(raster)
assign(paste0(var) , raster)
}
rm (raster)
gc() ; gc()

print ("make stack layer")
stack = get(names(table)[2] )
stack

for (var in names(table)[-2:-1] ) {
print(var)
stack =  c(stack,get(var))
}

print ( c(xmin, xmax  , ymin , ymax) )

# summary(stack)

extent  <- terra::ext( xmin , xmax , ymin , ymax )
extent
env=terra::crop(stack,extent)

print ("env  info")
setMinMax(env)
print ("env  str")
class(env)
str(env)
 
# rm (stack)
gc() ;  gc() ;  
print ("create the table")
                      
save.image(paste0("../vector_seed",seed,"/data2_",xmin,"_",ymin,".RData"))
env_predP = terra::predict(env,  model =  learnerP,  predict_type = "prob" ,  fun = predict )

terra::writeRaster (env_predP, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction_seed",seed,"/prediction_seed",seed,"P_",xmin,"_",ymin,".tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Float32")

env_predC = terra::predict(env , model = learnerR,   predict_type = "response" , fun = predict )

terra::writeRaster (env_predC, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction_seed",seed,"/prediction_seed",seed,"R_",xmin,"_",ymin,".tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Byte")

save.image(paste0("../vector_seed",seed,"/data2_",xmin,"_",ymin,".RData"))

'
rm -f  $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_${xmin}_${ymin}.tif.aux.xml  $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_${xmin}_${ymin}.tif.aux.xml

if [ $SLURM_ARRAY_TASK_ID -eq 41  ] ; then
sleep 2000
module purge
source ~/bin/gdal3
source ~/bin/pktools

rm -f $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all.* $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all_1km.*  $ONCHO/prediction_seed${seed}/prediction_seed${seed}?_all.* $ONCHO/prediction_seed${seed}/prediction_seed${seed}?_all_1km.*

gdalbuildvrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.vrt $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_[0-9]*_[0-9]*.tif
gdal_translate  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.vrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.tif 

gdalbuildvrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.vrt $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_[0-9]*_[0-9]*.tif
gdal_translate  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.vrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.tif 

pksetmask -m $ONCHO/input/hydrography90m/accumulation.tif -co COMPRESS=DEFLATE -co ZLEVEL=9 -msknodata -2147483648 -nodata -9999 -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.tif -o $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_msk.tif 

pksetmask -m $ONCHO/input/hydrography90m/accumulation.tif -co COMPRESS=DEFLATE -co ZLEVEL=9 -msknodata -2147483648 -nodata 255  -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.tif -o $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_msk.tif 

cd $ONCHO/prediction_seed${seed}/
rm -f gdaltindex $ONCHO/prediction_seed${seed}/all_tif_shp.*
gdaltindex $ONCHO/prediction_seed${seed}/all_tif_shp.shp  prediction_seed${seed}R_[0-9]*_[0-9]*.tif

gdal_translate -tr 0.00833333333333 0.00833333333333 -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.tif $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_1km.tif

pksetmask -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9   -msknodata -9999 -nodata  -9999 -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_1km.tif -o   $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_1km_msk.tif
# gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9 -scale 0.344 0.966 0.01 1 $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all_1km_msk.tif $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all_1km_msk_s.tif

gdal_translate -tr 0.00833333333333 0.00833333333333 -r mode -co COMPRESS=DEFLATE -co ZLEVEL=9 $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.tif $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_1km.tif
pksetmask -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9   -msknodata -9999 -nodata 255  -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_1km.tif -o   $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_1km_msk.tif
# gdal_translate -co COMPRESS=DEFLATE -co ZLEVEL=9 -scale 0.344 0.966 0.01 1 $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all_1km_msk.tif $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all_1km_msk_s.tif

rm $ONCHO/prediction_seed${seed}/*.tif.aux.xml
fi
