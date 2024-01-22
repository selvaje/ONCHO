#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 4:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc14_extract_confusionMatrix_clean.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc14_extract_confusionMatrix_clean.sh.%J.err
#SBATCH --job-name=sc14_extract_confusionMatrix_clean.sh 
#SBATCH --mem=10G


#####   sbatch   /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc14_extract_confusionMatrix_clean.sh


module load StdEnv

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
cd $ONCHO/vector


module load GDAL/3.6.2-foss-2022b

echo "obs pred" > NigeriaHabitatSites_pa_obs_pred_uniq_noheader_noclean.txt 
paste -d " " <(awk '{print $3}' NigeriaHabitatSites_x_y_pa_uniq_noheader.txt ) \
      <(gdallocationinfo -valonly -geoloc ../prediction_all/prediction_seedP_all_msk_mean_nokernel.tif  < <(awk '{print $1,$2}' NigeriaHabitatSites_x_y_pa_uniq_noheader.txt )) | awk '{ if($2<0.5) {print $1, 0 } else {print $1, 1 }  }'  >> NigeriaHabitatSites_pa_obs_pred_uniq_noheader_noclean.txt 

echo "obs pred" > abridged_consolidated_pa_obs_pred_uniq_noheader_flypresent.txt
paste -d " " <(awk '{print $3}' abridged_consolidated_x_y_pa_uniq_noheader_flypresent.txt) \
      <(gdallocationinfo -valonly -geoloc ../prediction_all/prediction_seedP_all_msk_mean_flypresent.tif  < <(awk '{print $1,$2}'    abridged_consolidated_x_y_pa_uniq_noheader_flypresent.txt ) ) | awk '{ if($2<0.5) {print $1, 0 } else {print $1, 1 }  }'  >> abridged_consolidated_pa_obs_pred_uniq_noheader_flypresent.txt

echo "obs pred" > abridged_consolidated_pa_obs_pred_uniq_noheader_adultpresent.txt
paste -d " " <(awk '{print $3}' abridged_consolidated_x_y_pa_uniq_noheader_adultpresent.txt) \
      <(gdallocationinfo -valonly -geoloc ../prediction_all/prediction_seedP_all_msk_mean_1km_adultpresent.tif  < <(awk '{print $1,$2}'    abridged_consolidated_x_y_pa_uniq_noheader_adultpresent.txt ))  | awk '{ if($2<0.5) {print $1, 0 } else {print $1, 1 }  }' >> abridged_consolidated_pa_obs_pred_uniq_noheader_adultpresent.txt

echo "obs pred" > abridged_consolidated_pa_obs_pred_uniq_noheader_larvaepresent.txt
paste -d " " <(awk '{print $3}' abridged_consolidated_x_y_pa_uniq_noheader_larvaepresent.txt) \
      <(gdallocationinfo -valonly -geoloc ../prediction_all/prediction_seedP_all_msk_mean_1km_larvaepresent.tif < <(awk '{print $1,$2}'    abridged_consolidated_x_y_pa_uniq_noheader_larvaepresent.txt )) | awk '{ if($2<0.5) {print $1, 0 } else {print $1, 1 }  }' >> abridged_consolidated_pa_obs_pred_uniq_noheader_larvaepresent.txt

module load R/4.3.0-foss-2022b   

Rscript --vanilla --verbose   -e '
library(caret)

fly = read.table("NigeriaHabitatSites_pa_obs_pred_uniq_noheader_noclean.txt", header = TRUE, sep = " ")
fly$pred = as.factor(fly$pred)
fly$obs  = as.factor(fly$obs)
# Creating confusion matrix
matrix_fly <- confusionMatrix(fly$pred, fly$obs)
capture.output( matrix_fly , file = "NigeriaHabitatSites_pa_obs_pred_uniq_noheader_noclean_confusionMatrix.txt") 

fly = read.table("abridged_consolidated_pa_obs_pred_uniq_noheader_flypresent.txt", header = TRUE, sep = " ")
fly$pred = as.factor(fly$pred)
fly$obs  = as.factor(fly$obs)
# Creating confusion matrix
matrix_fly <- confusionMatrix(fly$pred, fly$obs)
capture.output( matrix_fly , file = "abridged_consolidated_pa_obs_pred_uniq_noheader_flypresent_confusionMatrix.txt") 

fly = read.table("abridged_consolidated_pa_obs_pred_uniq_noheader_adultpresent.txt", header = TRUE, sep = " ")
fly$pred = as.factor(fly$pred)
fly$obs  = as.factor(fly$obs)
# Creating confusion matrix
matrix_fly <- confusionMatrix(fly$pred, fly$obs)
capture.output( matrix_fly , file = "abridged_consolidated_pa_obs_pred_uniq_noheader_adultpresent_confusionMatrix.txt") 

fly = read.table("abridged_consolidated_pa_obs_pred_uniq_noheader_larvaepresent.txt", header = TRUE, sep = " ")
fly$pred = as.factor(fly$pred)
fly$obs  = as.factor(fly$obs)
# Creating confusion matrix
matrix_fly <- confusionMatrix(fly$pred, fly$obs)
capture.output( matrix_fly , file = "abridged_consolidated_pa_obs_pred_uniq_noheader_larvaepresent_confusionMatrix.txt") 

q()
'

