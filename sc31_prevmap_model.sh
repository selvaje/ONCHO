#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /vast/palmer/scratch/sbsc/ga254/stdout/sc31_prevmap_model.sh.%J.out 
#SBATCH -e /vast/palmer/scratch/sbsc/ga254/stderr/sc31_prevmap_model.sh.%J.err
#SBATCH --job-name=sc31_prevmap_model.sh
#SBATCH --mem=300G

#### sbatch /vast/palmer/home.grace/ga254/scripts_gitreps/ONCHO/sc31_prevmap_model.sh 

ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO


module  load R/4.1.0-foss-2020b

# # # see http://www.css.cornell.edu/faculty/dgr2/_static/files/R_html/CompareRandomForestPackages.html
# # # R  --vanilla --no-readline   -q  <<'EOF'  this is not working with ranger 

# first RF for sorting the most important variables

Rscript --vanilla --verbose   -e '

library(terra)
library(PrevMap); 
library(splancs)
library(geoR)

xypa = read.table("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/x_y_pa_prob.txt", header = TRUE, sep = " ")
xypa$prob <- NULL

fit.MLE <- linear.model.MLE(  pa~1 ,  coords=~x+y,  start.cov.pars = c(10,0.1), xypa ,kappa=0.5)
# diagnostic<-variog.diagnostic.lm(fit.MLE)

x=rast("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/latlon/x_1km.tif"  )
y=rast("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/latlon/y_1km.tif"  )
coormatrix=as.matrix(c(x,y))     
colnames(coormatrix) = c("x","y")


save.image("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/prevmap/prevmap.Data")



pred.MLE1 <- spatial.pred.linear.MLE(fit.MLE,
                                    grid.pred = head(coormatrix, 100000  ),
                                    type="marginal",
                                    n.sim.prev = 1,
                                    standard.errors = FALSE)

pred.MLE2 <- spatial.pred.linear.MLE(fit.MLE,
                                    grid.pred = tail(coormatrix, 1029600  ),
                                    type="marginal",
                                    n.sim.prev = 1,
                                    standard.errors = FALSE)

save.image("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/prevmap/prevmap.RData")

pred.MLE.rast  = rast( matrix ( c(pred.MLE1$prevalence$predictions, pred.MLE2$prevalence$predictions), ncol=1560, nrow=1320))

terra::writeRaster (pred.MLE.rast, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/prevmap/linearmodelMLE.tif"), gdal=c("COMPRESS=DEFLATE","ZLEVEL=9"), overwrite=TRUE , datatype="Float32" , NAflag=-9999)
save.image("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/prevmap/prevmap.RData")

'


