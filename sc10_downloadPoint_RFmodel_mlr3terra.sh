#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc10_downloadPoint_RFmodel_mlr3terra.sh.%A_%a.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc10_downloadPoint_RFmodel_mlr3terra.sh.%A_%a.err
#SBATCH --job-name=sc10_downloadPoint_RFmodel_mlr3terra.sh
#SBATCH --mem=80G
#SBATCH --array=1-41
## array=1-41

#####  for seed in $(seq 1 100);  do sbatch  --export=seed=$seed /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc10_downloadPoint_RFmodel_mlr3terra.sh ; done 

#  x 2 15 
#  y 4 15                                                                                                              remove the first row that includes only see area
# for x in $(seq 2 2 14) ; do for y in $(seq 4 2 14 ) ; do echo $x $(expr $x + 2 ) $y $(expr $y + 2 ) ; done ; done  | awk '{ if (NR>1) print  }'  >   $ONCHO/vector/tile_list.txt  

module load StdEnv

echo "seed = $seed "

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
cd $ONCHO/vector

##  SLURM_ARRAY_TASK_ID=1
export seed=$seed

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

if [ $SLURM_ARRAY_TASK_ID -eq 1  ] ; then

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

Rscript --vanilla --verbose /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc10_downloadPoint_RFmodel_mlr3terra_R1.sh $seed 

else 
sleep 1500
echo "no sleep" 
fi ### close the first array loop 


module purge
module --ignore_cache load netCDF/4.7.4-gompi-2020b
module load GDAL/3.6.2-foss-2022b 
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
for seq  in $(seq 1 100) ; do echo $seq $( pkstat --hist -i prediction_seed${seq}/prediction_seed${seq}R_all_1km_msk.tif   | grep -v " 0" | grep -e "255 " ) ; done | grep -v 134625  | awk '{printf ("%i " , $1) }'

for seq  in prediction_seed??/prediction_seed8R_all_1km_msk.tif  ; do echo $seq $( pkstat --hist -i $file  | grep -v " 0" | grep -e "255 " ) ; done | grep -v 134625

grep OOB   vector_seed*/allVarR.mod.rf.txt    | awk ' { sum = sum + $4 } END {print sum / 100}' 
grep OOB   vector_seed*/selvsVarR.mod.rf_seed*.txt  | awk ' { sum = sum + $4 } END {print sum / 100}' 
grep OOB   vector_seed*/selvsVarP.mod.rf_seed*.txt  | awk ' { sum = sum + $6 } END {print sum / 100}' 


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
