#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc32_prevmap_model.sh.%A_%a.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc32_prevmap_model.sh.%A_%a.err
#SBATCH --job-name=sc32_prevmap_model.sh
#SBATCH --mem=200G
#SBATCH --array=2
## array=1-138   2 for testing 

#####   sbatch  /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc32_prevmap_model.sh  

#  x 2 15 
#  y 4 15                                                                                                              remove the first row that includes only see area
# for x in $(seq 2 1 14) ; do for y in $(seq 4 1 14 ) ; do echo $x $(expr $x + 1 ) $y $(expr $y + 1 ) ; done ; done   }'  >   $ONCHO/vector/tile_list1d.txt  
# remove manualy line with full sea area 1  2 12  13  23  $ONCHO/vector/tile_list1d.txt 

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

# if [ $SLURM_ARRAY_TASK_ID -eq 1  ] ; then
source ~/bin/gdal3
source ~/bin/pktools

echo "#######################################################"
echo "############### FIRST MLE  #############################"
echo "#######################################################"

#### module --ignore_cache load R/4.1.0-foss-2020b

# # # see http://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/CompareRandomForestPackages.html
# # # R  --vanilla --no-readline   -q  <<'EOF'  this is not working with ranger 

# first RF for sorting the most important variables

# Rscript --vanilla --verbose   -e '

# library(terra)
# library(PrevMap); 
# library(splancs)
# library(geoR)

# xypaprob  = read.table("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/x_y_pa_prob.txt", header = TRUE, sep = " ")

# fit.MLE <- linear.model.MLE(  pa~prob  ,  coords=~x+y,  start.cov.pars = c(10,0.1), data=xypaprob  ,kappa=0.5)

# save.image("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/mle/dataMLE0.RData")

# q()
# '
# else 
# sleep 300 
# fi ### close the first array loop 

module purge
module --ignore_cache load netCDF/4.7.4-gompi-2020b
source ~/bin/gdal3
module --ignore_cache load R/4.1.0-foss-2020b

echo geo_string  =  $xmin  $xmax $ymin $ymax

echo " make the RF prediction "

Rscript --vanilla  -e   '

library(terra) 
library(PrevMap)
library(splancs)
library(geoR)

xmin <- as.numeric(Sys.getenv("xmin"))
xmax <- as.numeric(Sys.getenv("xmax"))
ymin <- as.numeric(Sys.getenv("ymin"))
ymax <- as.numeric(Sys.getenv("ymax"))

xmin
xmax
ymin
ymax


load("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/mle/dataMLE0.RData")

extent   <- terra::ext( xmin , xmax , ymin , ymax )
extent

x=terra::crop(terra::rast("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/latlon/x.tif"), extent )
y=terra::crop(terra::rast("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/latlon/y.tif"), extent )

coormatrix=as.matrix(c(x,y))     
colnames(coormatrix) = c("x","y")


prob  = as.data.frame(terra::crop(terra::rast ("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif"),extent))
colnames(prob) = c("prob")
prob[is.na(prob)] <- -1 

print ("predict ghe mle model")

pred.MLE <- spatial.pred.linear.MLE(fit.MLE, 
                                    predictors = head (prob , 1000000 ) , 
                                    grid.pred  = head (coormatrix , 1000000 ) , 
                                    type="marginal",
                                    n.sim.prev = 1,
                                    standard.errors = FALSE)

str(pred.MLE)

pred.MLE.rast  = rast(matrix(pred.MLE$prevalence$predictions, ncol=1000, nrow=1000))
terra::writeRaster(pred.MLE.rast, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/mle/predtile/mlepprb",xmin,"_",ymin,".tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Float32" , NAflag=-9999)


save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/mle/dataMLE2_",xmin,"_",ymin,".RData"))

'

exit 

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


for HEAD in $(seq 1 143 ) ; do 

geo_string=$(head  -n  $HEAD $ONCHO/vector/tile_list1d.txt   | tail  -1 )
export xmin=$( echo $geo_string | awk '{  print $1 }' )
export xmax=$( echo $geo_string | awk '{  print $2 }' )
export ymin=$( echo $geo_string | awk '{  print $3 }' )
export ymax=$( echo $geo_string | awk '{  print $4 }' )

gdalbuildvrt -overwrite  -te $xmin $ymin $xmax $ymax test.vrt    /gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif 
echo $HEAD    $(pkinfo  -mm -i test.vrt )
done | grep "max -9999" | awk '{ print $1 }' 

