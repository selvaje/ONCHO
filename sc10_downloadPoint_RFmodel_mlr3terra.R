
library("randomForest")
library("varSelRF")

seed=1 
set.seed(seed)

table = read.table("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/x_y_pa_predictors4R.txt", header = TRUE, sep = " ")
table$ER2017 =   as.factor(table$ER2017)
table$LC2021 =   as.factor(table$LC2021)
table$pa     =   as.factor(table$pa)

des.table = summary(table)
write.table(des.table, "/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/stat_allVar.txt", quote = FALSE  )

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

#### variable selection base on varSelRF 
rf.vs =  varSelRF(table[-c(1,1)]  , table$pa , ntree = 500,  mtryFactor=11 ,    ntreeIterat = 500, vars.drop.frac = 0.1)
rf.vs
table.rf.vs = subset(table, select = rf.vs$selected.vars)
table.rf.vs$pa = table$pa 

#####


library("ranger")
# library(tuneRanger) ## for now not implemented

# fit a model with variable selection 
mod.rfP.vs = ranger( pa ~ . , table.rf.vs  , probability = TRUE  ,  classification=TRUE ,   importance="permutation")
mod.rfR.vs = ranger( pa ~ . , table.rf.vs  , probability = FALSE ,  classification=TRUE ,   importance="permutation")

impP.vs=as.data.frame(importance(mod.rfP.vs))
impR.vs=as.data.frame(importance(mod.rfR.vs))
save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data0.RData"))

impP.vs.s = impP.vs[order(impP.vs$"importance(mod.rfP.vs)",decreasing=TRUE), , drop = FALSE]
impR.vs.s = impR.vs[order(impR.vs$"importance(mod.rfR.vs)",decreasing=TRUE), , drop = FALSE]

write.table(impP.vs.s, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/importanceP_selvsVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfP = capture.output(mod.rfP.vs)
write.table(s.mod.rfP, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/selvsVarP.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

write.table(impR.vs.s, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/importanceR_selvsVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfR = capture.output(mod.rfR.vs)
write.table(s.mod.rfR, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/selvsVarR.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

# fit a model with all variables
mod.rfP.all = ranger( pa ~ . , table  , probability = TRUE   ,  classification=TRUE ,   importance="permutation")
mod.rfR.all = ranger( pa ~ . , table  , probability = FALSE  ,  classification=TRUE ,   importance="permutation")

impP.all=as.data.frame(importance(mod.rfP.all))
impR.all=as.data.frame(importance(mod.rfR.all))
save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data0.RData"))

impP.all.s = impP.all[order(impP.all$"importance(mod.rfP.all)",decreasing=TRUE), , drop = FALSE]
impR.all.s = impR.all[order(impR.all$"importance(mod.rfR.all)",decreasing=TRUE), , drop = FALSE]

write.table(impP.all.s, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/importanceP_allVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfP.all = capture.output(mod.rfP.all)
write.table(s.mod.rfP.all, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/allVarP.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

write.table(impR.all.s, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/importanceR_allVar_seed",seed,".txt"), quote = FALSE  )
s.mod.rfR.all = capture.output(mod.rfR.all)
write.table(s.mod.rfR.all, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/allVarR.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )


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

tablexy = read.table("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector/NigeriaHabitatSites_x_y_uniq_header.txt", header = TRUE, sep = " ")
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

save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data1.RData"))

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

save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data1.RData"))

# keep only the columns you need
score_spcv_rfc              =  as.data.frame(score_repeated_cv$iteration)
score_spcv_rfc$cv           =      score_repeated_cv$classif.acc
score_spcv_rfc$spcv_coords  =      score_repeated_spcv_coords$classif.acc 
##  score_spcv_rfc$spcv_env     =  score_repeated_spcv_env$classif.acc 
##  score_spcv_rfc$spcv_disc    =  score_repeated_spcv_disc$classif.acc 

write.table(score_spcv_rfc, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/cvP_selvsVar_seed",seed,".txt"), quote = FALSE , row.names = FALSE )
save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data1.RData"))
q()



library(ranger)
seed <- as.numeric(Sys.getenv("seed"))
load(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data1.RData"))
impP=as.data.frame(importance(learnerP$model))
impR=as.data.frame(importance(learnerR$model))

impP.s = impP[order(impP$"importance(learnerP$model)",decreasing=TRUE), , drop = FALSE]
impR.s = impR[order(impR$"importance(learnerR$model)",decreasing=TRUE), , drop = FALSE]

write.table(impP.s, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/importanceP_selvsVar_seed",seed,".txt"), quote = FALSE  )
write.table(impR.s, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/importanceR_selvsVar_seed",seed,".txt"), quote = FALSE  )

conf.matrix = learnerR$model$confusion.matrix
pred.error =  learnerR$model$prediction.error 

write.table(pred.error , paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/pred_error_seed",seed,".txt"))
write.table(conf.matrix, paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/conf_matrix_seed",seed,".txt"))

write.table(capture.output(learnerR$model), paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/selvsVarR.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )
write.table(capture.output(learnerP$model), paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/selvsVarP.mod.rf_seed",seed,".txt"), quote = FALSE , row.names = FALSE )

save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data1.RData"))



library("rlang")
library("mlr3")
library("mlr3spatial")
library("mlr3learners")
library("mlr3spatiotempcv")
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

load(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data1.RData"))

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
raster  = terra::rast (Sys.glob(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/input/*/",var,".tif")  ) )

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
                      
save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data2_",xmin,"_",ymin,".RData"))
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

save.image(paste0("/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO/vector_seed",seed,"/data2_",xmin,"_",ymin,".RData"))

