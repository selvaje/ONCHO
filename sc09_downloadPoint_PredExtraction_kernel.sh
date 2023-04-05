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

source ~/bin/gdal3
source ~/bin/pktools
source ~/bin/grass78m
 
grass78  -f -text --tmp-location  -c $ONCHO/input/geomorpho90m/elevation.tif     <<'EOF'

for radius in 1.5 2 2.5 3 ; do 
r.external  input=$ONCHO/input/geomorpho90m/elevation.tif output=elv --overwrite
g.region res=0:00:30
v.in.ogr  input=$ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq.gpkg output=pa   where="PA = 0" --o
v.kernel input=pa output=kernel   radius=$radius   multiplier=1 --o
g.region res=0:00:03
r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=-9999 type=Float32 format=GTiff input=kernel output=$ONCHO/input/pa_kernel/kernel0_tmp.tif
pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $ONCHO/input/geomorpho90m/elevation.tif  -msknodata -9999 -nodata -9999 -i $ONCHO/input/pa_kernel/kernel0_tmp.tif -o $ONCHO/input/pa_kernel/kernel0_R$radius.tif
g.region res=0:00:30
v.in.ogr  input=$ONCHO/vector/NigeriaHabitatSites_x_y_pa_uniq.gpkg output=pa   where="PA = 1" --o 
v.kernel input=pa output=kernel   radius=$radius  multiplier=1 --o
g.region res=0:00:03
r.out.gdal --o -f -c -m createopt="COMPRESS=DEFLATE,ZLEVEL=9" nodata=-9999 type=Float32 format=GTiff input=kernel output=$ONCHO/input/pa_kernel/kernel1_tmp.tif
pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m $ONCHO/input/geomorpho90m/elevation.tif  -msknodata -9999 -nodata -9999 -i $ONCHO/input/pa_kernel/kernel1_tmp.tif -o $ONCHO/input/pa_kernel/kernel1_R$radius.tif
done 
EOF

rm $ONCHO/input/pa_kernel/kernel0_tmp.tif $ONCHO/input/pa_kernel/kernel1_tmp.tif

