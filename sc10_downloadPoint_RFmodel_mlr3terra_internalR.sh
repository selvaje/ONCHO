#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc10_downloadPoint_RFmodel_mlr3terra.sh.%A_%a.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc10_downloadPoint_RFmodel_mlr3terra.sh.%A_%a.err
#SBATCH --job-name=sc10_downloadPoint_RFmodel_mlr3terra.sh
#SBATCH --mem=80G
#SBATCH --array=1
## array=1-41

#####  for seed in $(seq 1 100);  do sbatch  --export=seed=$seed /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc10_downloadPoint_RFmodel_mlr3terra.sh ; done 

#  x 2 15 
#  y 4 15                                                                                                              remove the first row that includes only see area
# for x in $(seq 2 2 14) ; do for y in $(seq 4 2 14 ) ; do echo $x $(expr $x + 2 ) $y $(expr $y + 2 ) ; done ; done  | awk '{ if (NR>1) print  }'  >   $ONCHO/vector/tile_list.txt  

module load StdEnv

echo "seed = $seed "

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
cd $ONCHO/vector

SLURM_ARRAY_TASK_ID=1
export seed=1

geo_string=$(head  -n  $SLURM_ARRAY_TASK_ID $ONCHO/vector/tile_list.txt   | tail  -1 )
export xmin=$( echo $geo_string | awk '{  print $1 }' ) 
export xmax=$( echo $geo_string | awk '{  print $2 }' ) 
export ymin=$( echo $geo_string | awk '{  print $3 }' ) 
export ymax=$( echo $geo_string | awk '{  print $4 }' ) 

echo geo_string  =  $xmin  $xmax $ymin $ymax
echo prediction_${xmin}_${ymin}.tif

echo $seed 

rm -f $ONCHO/vector_seed$seed/data*.RData
rm -f $ONCHO/vector_seed$seed/*.txt
rm -f $ONCHO/prediction_seed$seed/prediction_*_*.tif
rm -f $ONCHO/prediction_seed$seed/*.tif.aux.xml

# if [ $SLURM_ARRAY_TASK_ID -eq 1  ] ; then

module load GDAL/3.6.2-foss-2022b 
module load PKTOOLS/2.6.7.6-foss-2020b

echo "#######################################################"
echo "############### FIRST RF #############################"
echo "#######################################################"

##  module --ignore_cache load R/4.1.0-foss-2020b
module load R/4.3.0-foss-2022b   

# # # see http://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/CompareRandomForestPackages.html
# # # R  --vanilla --no-readline   -q  <<'EOF'  this is not working with ranger 

# first RF for sorting the most important variables
### Rscript --vanilla --verbose   -e '

Rscript --vanilla --verbose /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc10_downloadPoint_RFmodel_mlr3terra_R1.sh

# > /vast/palmer/scratch/sbsc/ga254/stderr/sc10_downloadPoint_RFmodel_mlr3terra.sh_errpr.err  2>&1 

Rscript --vanilla --verbose   -e '
options(echo=T)

library("randomForest")
library("varSelRF")

seed <- as.numeric(Sys.getenv("seed"))
set.seed(1)

table = read.table("x_y_pa_predictors4R_flypresent.txt", header = TRUE, sep = " ")
table$ER2017 =   as.factor(table$ER2017)
table$LC2021 =   as.factor(table$LC2021)
table$pa =    as.factor(table$pa)

des.table = summary(table)
write.table(des.table, "stat_allVar.txt", quote = FALSE  )

# # nomralization 
# library("broman")
# table_norm = as.data.frame(normalize(table[-c(1,1)]))
# var_names = names(table[-c(1,1)])
# names(table_norm) <- var_names
# table_norm = cbind(table[c(1)] , table_norm )

# table_norm$ER2017 =   table$ER2017
# table_norm$LC2021 =   table$LC2021
# table_norm$pa =       table$pa
# table = table_norm 


# t <- tuneRF(table[,-1], table[,1],  stepFactor = 1,  plot = TRUE, ntreeTry = 500, trace = TRUE,  improve = 0.01)  

print("variable selection base on varSelRF") 
rf.vs =  varSelRF(table[-c(1,1)]  , table$pa , ntree = 500,  mtryFactor=11 ,    ntreeIterat = 500, vars.drop.frac = 0.1)
rf.vs
table.rf.vs = subset(table, select = rf.vs$selected.vars)
table.rf.vs$pa = table$pa 

print("start the ranger implementation")

library(ranger)
# library(tuneRanger) ## for now not implemented

print("fit a model with variable selection")
mod.rfP.vs = ranger( pa ~ . , table.rf.vs  , probability = TRUE  ,  classification=TRUE ,   importance="permutation")
mod.rfR.vs = ranger( pa ~ . , table.rf.vs  , probability = FALSE ,  classification=TRUE ,   importance="permutation")

impP.vs=as.data.frame(importance(mod.rfP.vs))
impR.vs=as.data.frame(importance(mod.rfR.vs))
save.image(paste0("../vector_seed",seed,"/data0.RData"))

impP.vs.s = impP.vs[order(impP.vs$"importance(mod.rfP.vs)",decreasing=TRUE), , drop = FALSE]
impR.vs.s = impR.vs[order(impR.vs$"importance(mod.rfR.vs)",decreasing=TRUE), , drop = FALSE]

write.table(impP.vs.s, paste0("../vector_seed",seed,"/importanceP_selvsVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfP = capture.output(mod.rfP.vs)
write.table(s.mod.rfP, paste0("../vector_seed",seed,"/selvsVarP.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

write.table(impR.vs.s, paste0("../vector_seed",seed,"/importanceR_selvsVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfR = capture.output(mod.rfR.vs)
write.table(s.mod.rfR, paste0("../vector_seed",seed,"/selvsVarR.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

print("fit a model with all variables")
mod.rfP.all = ranger( pa ~ . , table  , probability = TRUE   ,  classification=TRUE ,   importance="permutation")
mod.rfR.all = ranger( pa ~ . , table  , probability = FALSE  ,  classification=TRUE ,   importance="permutation")

impP.all=as.data.frame(importance(mod.rfP.all))
impR.all=as.data.frame(importance(mod.rfR.all))
save.image(paste0("../vector_seed",seed,"/data0.RData"))

impP.all.s = impP.all[order(impP.all$"importance(mod.rfP.all)",decreasing=TRUE), , drop = FALSE]
impR.all.s = impR.all[order(impR.all$"importance(mod.rfR.all)",decreasing=TRUE), , drop = FALSE]

write.table(impP.all.s, paste0("../vector_seed",seed,"/importanceP_allVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfP.all = capture.output(mod.rfP.all)
write.table(s.mod.rfP.all, paste0("../vector_seed",seed,"/allVarP.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

write.table(impR.all.s, paste0("../vector_seed",seed,"/importanceR_allVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfR.all = capture.output(mod.rfR.all)
write.table(s.mod.rfR.all, paste0("../vector_seed",seed,"/allVarR.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

print("spatial validation procedure")

library(rlang)
library(mlr3spatiotempcv)  # spatio-temporal resampling 
library(mlr3tuning)        # hyperparameter tuning package
library("mlr3learners")
library("mlr3misc")
library("stars")
library("terra")
library("future")
library("blockCV")
library("sf")
######################### 

tablexy = read.table("NigeriaHabitatSites_x_y_uniq_header.txt", header = TRUE, sep = " ")
table.rf.vs$x = tablexy$x
table.rf.vs$y = tablexy$y

# # create task  for cross validation    https://mlr3book.mlr-org.com/special.html  
task = mlr3spatiotempcv::TaskClassifST$new(
   id = "identifier_table.rf.vs",
   backend = mlr3::as_data_backend(table.rf.vs), 
   target = "pa",
   positive = "1",
   coordinate_names = c("x", "y"),
   coords_as_features = FALSE,
   crs = 4326
   )
print(task)

# backend = as_data_backend(table.rf.vs)    # this is just table for the learner  
# task = as_task_classif(backend, target = "pa" ,  positive = "1" )

### https://mlr3extralearners.mlr-org.com/articles/learners/list_learners.html 
learnerP = lrn("classif.ranger", predict_types = "prob" ,  predict_type = "prob" , importance = "permutation", properties = "twoclass" ) 
learnerP$parallel_predict = TRUE  
print(learnerP)
learnerP$train(task)  # usefull to obtain the $model
learnerP$model

learnerR = lrn("classif.ranger", predict_types = "response", predict_type = "response"  ,  importance = "permutation", properties = "twoclass" )
learnerR$parallel_predict = TRUE  
print(learnerR)
learnerR$train(task)
learnerR$model

#### get  list of resamplint tecnique via  as.data.table(mlr_resamplings)
#### https://geocompr.robinlovelace.net/spatial-cv.html
#### resampling methods ---------------------------------------------------------
#### https://mlr-org.com/resamplings.html    https://mlr3spatiotempcv.mlr-org.com/articles/spatiotemp-viz.html 
#### as.data.table(mlr_resamplings)

rsmp_repeated_cv             = mlr3::rsmp("repeated_cv", folds = 10, repeats = 10)
rsmp_repeated_spcv_coords    = mlr3::rsmp("repeated_spcv_coords", folds = 10, repeats = 10)
## rsmp_repeated_spcv_env       = mlr3::rsmp("repeated_spcv_env",  folds = 10, repeats = 10, features= "CHELSA_bio4_r" )
## rsmp_repeated_spcv_disc      = mlr3::rsmp("repeated_spcv_disc", folds = 10, repeats = 10 ,  radius = 10L, buffer = 10L)

save.image(paste0("../vector_seed",seed,"/data1.RData"))

rsmp_repeated_cv$instantiate(task)
rsmp_repeated_spcv_coords$instantiate(task)
## rsmp_repeated_spcv_env$instantiate(task)
## rsmp_repeated_spcv_disc$instantiate(task)

# https://mlr3spatiotempcv.mlr-org.com/articles/mlr3spatiotempcv.html
# reduce verbosity 
## lgr::get_logger("mlr3")$set_threshold("warn")
# run spatial cross-validation and save it to resample result rf  (rr_rfc)

rr_repeated_cv             = mlr3::resample(task = task, learner = learnerP, resampling = rsmp_repeated_cv            , store_models = TRUE )
rr_repeated_spcv_coords    = mlr3::resample(task = task, learner = learnerP, resampling = rsmp_repeated_spcv_coords   , store_models = TRUE )
## rr_repeated_spcv_env    = mlr3::resample(task = task, learner = learnerR, resampling = rsmp_repeated_spcv_env      , store_models = TRUE )
## rr_repeated_spcv_disc      = mlr3::resample(task = task, learner = learnerR, resampling = rsmp_repeated_spcv_disc     , store_models = TRUE )


# compute the classification accuracy  as a data.table  # https://mlr3.mlr-org.com/reference/mlr_measures.html#ref-examples 

score_repeated_cv                      = rr_repeated_cv$score(measure = mlr3::msr("classif.acc"))
score_repeated_spcv_coords    = rr_repeated_spcv_coords$score(measure = mlr3::msr("classif.acc"))
## score_repeated_spcv_env       = rr_repeated_spcv_env$score(measure = mlr3::msr("classif.acc"))
## score_repeated_spcv_disc        = rr_repeated_spcv_disc$score(measure = mlr3::msr("classif.acc"))

save.image(paste0("../vector_seed",seed,"/data1.RData"))

# keep only the columns you need
score_spcv_rfc              =  as.data.frame(score_repeated_cv$iteration)
score_spcv_rfc$cv           =      score_repeated_cv$classif.acc
score_spcv_rfc$spcv_coords  =      score_repeated_spcv_coords$classif.acc 
##  score_spcv_rfc$spcv_env     =  score_repeated_spcv_env$classif.acc 
##  score_spcv_rfc$spcv_disc    =  score_repeated_spcv_disc$classif.acc 

write.table(score_spcv_rfc, paste0("../vector_seed",seed,"/cvP_selvsVar_seed",seed,".txt"), quote = FALSE , row.names = FALSE )
save.image(paste0("../vector_seed",seed,"/data1.RData"))
q()
'

Rscript --vanilla --verbose   -e '
library(ranger)
seed <- as.numeric(Sys.getenv("seed"))
load(paste0("../vector_seed",seed,"/data1.RData"))
impP=as.data.frame(importance(learnerP$model))
impR=as.data.frame(importance(learnerR$model))

impP.s = impP[order(impP$"importance(learnerP$model)",decreasing=TRUE), , drop = FALSE]
impR.s = impR[order(impR$"importance(learnerR$model)",decreasing=TRUE), , drop = FALSE]

write.table(impP.s, paste0("../vector_seed",seed,"/importanceP_selvsVar_seed",seed,".txt"), quote = FALSE  )
write.table(impR.s, paste0("../vector_seed",seed,"/importanceR_selvsVar_seed",seed,".txt"), quote = FALSE  )

conf.matrix = learnerR$model$confusion.matrix
pred.error =  learnerR$model$prediction.error 

write.table(pred.error , paste0("../vector_seed",seed,"/pred_error_seed",seed,".txt"))
write.table(conf.matrix, paste0("../vector_seed",seed,"/conf_matrix_seed",seed,".txt"))

write.table(capture.output(learnerR$model), paste0("../vector_seed",seed,"/selvsVarR.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )
write.table(capture.output(learnerP$model), paste0("../vector_seed",seed,"/selvsVarP.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

save.image(paste0("../vector_seed",seed,"/data1.RData"))
q()
'

else 
sleep 1500
echo "no sleep" 
fi ### close the first array loop 

exit 

module purge
module --ignore_cache load netCDF/4.7.4-gompi-2020b
module load GDAL/3.6.2-foss-2022b 
# module --ignore_cache load R/4.1.0-foss-2020b
module load R/4.3.0-foss-2022b   

echo geo_string  =  $xmin  $xmax $ymin $ymax

echo " make the RF prediction "

Rscript --vanilla  -e   '
library("rlang")
library("mlr3")
library("mlr3spatial")
library("mlr3learners")
library("mlr3spatiotempcv")
library(ranger)
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
print(learnerR)

print ("start the prediction")
# bb = st_bbox ( c(xmin = xmin , xmax =  xmax , ymin =  ymin, ymax =  ymax  )  , crs = 4326 ) 

table.rf.vs = subset(table.rf.vs, select =  -c(pa , x , y )) ### remove pa x y 

for (var in  names( table.rf.vs) )  {
print(var)
raster  = terra::rast (Sys.glob(paste0("../input/*/",var,".tif")  ) )

print(raster)
assign(paste0(var) , raster)
}
rm (raster)
gc() ; gc()

print ("make stack layer")
stack = get(names(table.rf.vs)[1] )
stack

for (var in names(table.rf.vs)[-1] ) {
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

env
 
# rm (stack)
gc() ;  gc() ;  
print ("create the table")
                      
save.image(paste0("../vector_seed",seed,"/data2_",xmin,"_",ymin,".RData"))
env_predP = terra::predict(env,  model =  learnerP,   predict_type = "prob" ,  fun = predict )

terra::writeRaster (env_predP, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction_seed",seed,"/prediction_seed",seed,"P_",xmin,"_",ymin,".tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Float32" , NAflag=-9999)

env_predR = terra::predict(env , model = learnerR,   predict_type = "response" , fun = predict )

terra::writeRaster (env_predR, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction_seed",seed,"/prediction_seed",seed,"R_",xmin,"_",ymin,".tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Byte" , NAflag=255)

# predict on table 
# newdata = as.data.frame(as.matrix(env))
# colSums(is.na(newdata))  # 0 NAs
# but assuming there were 0s results in a more generic approach
# ind = rowSums(is.na(newdata)) == 0
# tmp = learnerR$predict_newdata(newdata = newdata , task = task  , predict_type = "response" , fun = predict )  
# newdata$pred = data.table::as.data.table(tmp)[["response"]]
# pred_2 = env$elevation
# now fill the raster with the predicted values
# pred_2[] = newdata$pred
# check if terra and our manual prediction is the same
## all(value(env_predR - pred_2) == 0)

save.image(paste0("../vector_seed",seed,"/data2_",xmin,"_",ymin,".RData"))

'

rm -f  $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_${xmin}_${ymin}.tif.aux.xml  $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_${xmin}_${ymin}.tif.aux.xml

if [ $SLURM_ARRAY_TASK_ID -eq 41  ] ; then
sleep 2000
module purge

module load GDAL/3.6.2-foss-2022b                                                                                                                                                           
module load PKTOOLS/2.6.7.6-foss-2020b

rm -f $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all.* $ONCHO/prediction_seed${seed}/prediction_seed${seed}_all_1km.*  $ONCHO/prediction_seed${seed}/prediction_seed${seed}?_all.* $ONCHO/prediction_seed${seed}/prediction_seed${seed}?_all_1km.*

gdalbuildvrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.vrt $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_[0-9]*_[0-9]*.tif
gdal_translate  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.vrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.tif 

gdalbuildvrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.vrt $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_[0-9]*_[0-9]*.tif
gdal_translate -ot Byte  -co COMPRESS=DEFLATE -co ZLEVEL=9   $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.vrt  $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.tif 

pksetmask -m $ONCHO/input/hydrography90m/accumulation.tif -co COMPRESS=DEFLATE -co ZLEVEL=9 -msknodata -2147483648 -nodata -9999 -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.tif -o $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_msk.tif 

pksetmask -ot Byte   -m $ONCHO/input/hydrography90m/accumulation.tif -co COMPRESS=DEFLATE -co ZLEVEL=9 -msknodata -2147483648 -nodata 255  -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.tif -o $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_msk.tif 

cd $ONCHO/prediction_seed${seed}/
rm -f gdaltindex $ONCHO/prediction_seed${seed}/all_tif_shp.*
gdaltindex $ONCHO/prediction_seed${seed}/all_tif_shp.shp  prediction_seed${seed}R_[0-9]*_[0-9]*.tif

gdal_translate -tr 0.00833333333333 0.00833333333333 -r average -co COMPRESS=DEFLATE -co ZLEVEL=9 $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all.tif $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_1km.tif

pksetmask -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9   -msknodata -9999 -nodata  -9999 -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_1km.tif -o   $ONCHO/prediction_seed${seed}/prediction_seed${seed}P_all_1km_msk.tif

gdal_translate  -ot Byte -tr 0.00833333333333 0.00833333333333 -r mode -co COMPRESS=DEFLATE -co ZLEVEL=9 $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all.tif $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_1km.tif
pksetmask -ot Byte  -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9   -msknodata -9999 -nodata 255  -i $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_1km.tif -o   $ONCHO/prediction_seed${seed}/prediction_seed${seed}R_all_1km_msk.tif


rm $ONCHO/prediction_seed${seed}/*.tif.aux.xml
fi


exit 

# for controlling 
for seq  in $(seq 1 100) ; do echo $seq $( pkstat --hist -i prediction_seed${seq}/prediction_seed${seq}R_all_1km_msk.tif   | grep -v " 0" | grep -v "255 " ) ; done  | grep -e nan -e FileOpenError  | awk '{printf ("%i " , $1) }'

# check if all the tiles have been done correctly 
for seq  in $(seq 1 100) ; do echo $seq $( pkstat --hist -i prediction_seed${seq}/prediction_seed${seq}R_all_1km_msk.tif   | grep -v " 0" | grep -e "255 " ) ; done | grep -v 134625

for seq  in prediction_seed??/prediction_seed8R_all_1km_msk.tif  ; do echo $seq $( pkstat --hist -i $file  | grep -v " 0" | grep -e "255 " ) ; done | grep -v 134625

grep OOB   vector_seed*/allVarR.mod.rf.txt    | awk ' { sum = sum + $4 } END {print sum / 100}' 
grep OOB   vector_seed*/selvsVarR.mod.rf.txt  | awk ' { sum = sum + $4 } END {print sum / 100}' 
grep OOB   vector_seed*/selvrVarR.mod.rf.txt  | awk ' { sum = sum + $4 } END {print sum / 100}' 


#### 
gdalbuildvrt  $ONCHO/prediction_all/prediction_all.vrt $ONCHO/prediction_*/prediction_seed*R_all_msk.tif 
pksetmask -m  ../input/geomorpho90m/slope.tif  -co COMPRESS=DEFLATE -co ZLEVEL=9 


# 
### check if the R map  has been created 
for seq  in $(seq 1 100) ; do ll  prediction_seed${seq}/prediction_seed${seq}R_all_1km_msk.tif 2>> /tmp/test.txt  ; done
awk '{ gsub("seed", " " )  ; gsub("R_all", " " ) ;  printf ("%i " , $6)  }'  /tmp/test.txt

#### 




exit

install.packeges("randomForest" , "varSelRF" ,"broman" , "ranger" , "rlang" , "mlr3" , "mlr3spatiotempcv"  , "mlr3tuning" , "mlr3learners" ,  "mlr3misc" , "stars", "terra", "future", "blockCV", "sf","mlr3spatial" ,  "stars")
