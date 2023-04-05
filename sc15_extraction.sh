#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc15_extraction.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc15_extraction.sh.%J.err
#SBATCH --job-name=sc15_extraction.sh 
#SBATCH --mem=6G

#####  /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc15_extraction.sh 

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
source ~/bin/gdal3
source ~/bin/pktools

GDAL_CACHEMAX=30000

##### probability for NOEC
cd $ONCHO/
echo "prob_mean"  >   vector/23_1.23_EntoSurveillanceSites_Adeleke_mean.csv
gdallocationinfo -valonly -wgs84 $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif   < <(awk -F "," '{ if (NR>1) print $5 , $4  }'  vector/23_1.23_EntoSurveillanceSites_Adeleke.csv)  >>   vector/23_1.23_EntoSurveillanceSites_Adeleke_mean.csv
echo "prob_stdev" >   vector/23_1.23_EntoSurveillanceSites_Adeleke_stdev.csv
gdallocationinfo -valonly -wgs84   $ONCHO/prediction_all/prediction_seedP_all_msk_stdev.tif < <(awk -F "," '{ if (NR>1) print $5 , $4  }'  vector/23_1.23_EntoSurveillanceSites_Adeleke.csv) >>   vector/23_1.23_EntoSurveillanceSites_Adeleke_stdev.csv  

paste -d ","  $ONCHO/vector/23_1.23_EntoSurveillanceSites_Adeleke.csv $ONCHO/vector/23_1.23_EntoSurveillanceSites_Adeleke_mean.csv $ONCHO/vector/23_1.23_EntoSurveillanceSites_Adeleke_stdev.csv  > $ONCHO/vector/23_1.23_EntoSurveillanceSites_Adeleke_mean_stdev.csv

rm $ONCHO/vector/23_1.23_EntoSurveillanceSites_Adeleke_mean.csv  $ONCHO/vector/23_1.23_EntoSurveillanceSites_Adeleke_stdev.csv

###### prob for uniq 

echo "prob_mean"   >   $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_mean.txt
gdallocationinfo -valonly -wgs84 $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif   < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt   >>   $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_mean.txt 

echo "prob_stdev"  >   $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_stdev.txt
gdallocationinfo -valonly -wgs84 $ONCHO/prediction_all/prediction_seedP_all_msk_stdev.tif  < $ONCHO/vector/NigeriaHabitatSites_x_y_uniq_noheader.txt   >>   $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_stdev.txt 


paste -d " "   $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header.txt  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_mean.txt  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_stdev.txt  >  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_mean_stdev_uniq_header.txt 

rm  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_mean.txt  $ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq_header_stdev.txt 
