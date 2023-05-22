#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc12_residuals.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc12_residuals.sh.%J.err
#SBATCH --job-name=sc12_residuals.sh
#SBATCH --mem=10G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc12_residuals.sh

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO
source ~/bin/gdal3
source ~/bin/pktools
module load R/4.1.0-foss-2020b

cd $ONCHO/

echo "x y pa prob"  >  $ONCHO/vector/x_y_pa_prob.txt
paste -d " "    <( awk '{if(NR>1) print $1, $2 , $3 }' $ONCHO/vector/x_y_pa_predictors.txt)  \
<(gdallocationinfo -valonly -geoloc $ONCHO/prediction_all/prediction_seedP_all_msk_mean.tif  < <(awk '{if(NR>1) print $1, $2 }' $ONCHO/vector/x_y_pa_predictors.txt)) >> $ONCHO/vector/x_y_pa_prob.txt 


Rscript --vanilla  -e '

library(gstat)  
library(sp)

table = read.table("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/x_y_pa_prob.txt", header = TRUE, sep = " ")
table$residual = table$pa - table$prob
table$random <- sample(100, size = nrow(table), replace = TRUE)

coordinates(table)= ~ x+y
proj4string(table) = CRS("+proj=longlat +datum=WGS84")

var_residual=variogram(residual~1, data=table)
var_random=variogram(random~1, data=table)

pdf("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/variogram_residual.pdf" , width=6, height=5  )
plot(var_residual , ylim=c(0.14,0.20))
dev.off()

# pdf("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/variogram_random.pdf" , width=6, height=5  )
# plot(var_random)
# dev.off()


q()

'


