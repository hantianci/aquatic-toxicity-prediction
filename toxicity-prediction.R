# download the data 
data = read.csv("/Users/tapanasovich/Downloads/qsar_fish_toxicity.csv",  
header = FALSE, sep = ";") 

# data type assesment 
pairs(data)

library(gbm)

# create hyperparameter grid
hyper_grid <- expand.grid(
  shrinkage = c(0.01,0.005,0.001),
  interaction.depth = c(2, 4, 5),
  n.minobsinnode = 5,
  bag.fraction = c(.65, 0.8,1), 
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# total number of combinations
nrow(hyper_grid)
## [1] 27
 

# grid search 
for(i in 1:nrow(hyper_grid)) {
  
# reproducibility
set.seed(123)
  
# train model
gbm.tune <- gbm(
    formula = V7 ~ .,
    distribution = "gaussian",
    data = data,
    n.trees = 10000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    cv.folds = 10,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
# add min training error and trees to grid
hyper_grid$optimal_trees[i] <- which.min(gbm.tune$cv.error)

hyper_grid$min_RMSE[i] <- (min(gbm.tune$cv.error))
}
library(ggplot2)
library(dplyr)
 
hyper_grid %>% 
dplyr::arrange(min_RMSE) %>%
head(10)
  
# for reproducibility
set.seed(123)

# train GBM model
gbm.fit.final <- gbm(formula = V7 ~ .,distribution = "gaussian",
  data = data,n.trees = 3700, interaction.depth = 5,
  shrinkage = 0.005,n.minobsinnode = 5, bag.fraction = .65, 
  train.fraction = 1,n.cores = NULL, # will use all cores by default
  verbose = FALSE
  )  
 
par(mar = c(5, 8, 1, 1))
summary(gbm.fit.final, 
  method = relative.influence, # also can use permutation.test.gbm
  las = 2
  )
 

fit<-gbm.fit.final$fit
residuals<-data$V7- fit

library(ggplot2)
p1<- ggplot(data, aes( fit,  residuals)) +
  geom_point()
p1<-p1 + geom_smooth(method = "loess")      
p2<- ggplot(data, aes( fit,  residuals^2)) +
  geom_point()
p2<-p2 + geom_smooth(method = "loess")      

p3 <- ggplot(data, aes(sample=residuals))
p3<-p3 + stat_qq() + stat_qq_line()
grid.arrange(p1, p2, p3, ncol = 2)
   
library(pdp) 

vp1<-  gbm.fit.final %>%
partial(pred.var = "V1", n.trees = gbm.fit.final$n.trees, 
grid.resolution = 100) %>%
autoplot(rug = TRUE, train = data) 
 
vp2<- gbm.fit.final %>%
partial(pred.var = "V2", n.trees = gbm.fit.final$n.trees, 
grid.resolution = 100) %>%
autoplot(rug = TRUE, train = data) 
 
vp3<-  gbm.fit.final %>%
partial(pred.var = "V3", n.trees = gbm.fit.final$n.trees, 
grid.resolution = 100) %>%
autoplot(rug = TRUE, train = data)  
 
vp4<-  gbm.fit.final %>%
partial(pred.var = "V4", n.trees = gbm.fit.final$n.trees, 
grid.resolution = 100) %>%
autoplot(rug = TRUE, train = data) 
 
vp5<- gbm.fit.final %>%
partial(pred.var = "V5", n.trees = gbm.fit.final$n.trees, 
grid.resolution = 100) %>%
autoplot(rug = TRUE, train = data) 
 
vp6<-  gbm.fit.final %>%
partial(pred.var = "V6", n.trees = gbm.fit.final$n.trees, 
grid.resolution = 100) %>%
autoplot(rug = TRUE, train = data)  
grid.arrange(vp1, vp2, vp3,  vp4, vp5, vp6,   ncol = 3) 
  
#__________________________________________

library(vip) 
library(lattice)
library(viridis)

# Compute partial dependence data for lstat and rm
pd <- partial( gbm.fit.final, pred.var = c("V2", "V5"), n.trees = gbm.fit.final$n.trees  )

# Default PDP
pdp1 <- plotPartial(pd)

# Add contour lines and use a different color palette
rwb <- colorRampPalette(c("red", "white", "blue"))
pdp2 <- plotPartial(pd, contour = TRUE, col.regions = rwb)

# 3-D surface
pdp3 <- plotPartial(pd, levelplot = FALSE, zlab = "V7", colorkey = TRUE, 
                    screen = list(z = -20, x = -60))

# Figure 5
grid.arrange(pdp1, pdp2, pdp3, ncol = 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(lmSubsets)
MOS1_best <- lmSelect(V7 ~ poly(V1,3)+V1*V4+V1*V5+poly(V2,3)+ poly(V3,3) + poly(V6,3) +V2*V4+V2*V5+V3*V4+V3*V5, data=data,  penalty = "BIC", nbest = 20) 
  
coef( MOS1_best, size = 12)
 
## summary statistics 
summary(MOS1_best ) 
## visualize 
plot(MOS1_best )

#  poly(V1, 3) 2
# V4          try as factor                
# V5            0+1 vs rest                    
#  poly(V2, 3)    3         
#  poly(V3, 1)     1                
#  poly(V6, 3)     3          
#  V1:V4                       
#  V5:V2                       

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# GAM
library(mgcv)

function (formula, data, nfold = 10, debug.level = 0, 
          printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)
{
  if (is.null(cvparts)) {
    set.seed(seed)
    cvparts <- sample(1:nfold, nrow(data), replace = TRUE)
  }
  folds <- unique(cvparts)
  khat <- hat <- numeric(nrow(data))
  scale.gam <- summary(gam(formula, data = data))$scale
  for (i in folds) {
    trainrows <- cvparts != i
    testrows <- cvparts == i
    elev.gam <- gam(formula, data = data[trainrows, ], 
                    gamma = gamma)
    hat[testrows] <- predict(elev.gam, newdata = data[testrows,
                                                      ], select = TRUE)
    res <- residuals(elev.gam)
  }
  y <- eval(formula[[2]], envir = as.data.frame(data))
  res <- y - hat
  cvscale <- sum(res^2)/length(res)
  prntvec <- c(GAMscale = scale.gam, `CV-mse-GAM ` = cvscale)
  if (printit)
    print(round(prntvec, 4))
  invisible(list(fitted = hat, resid = res, cvscale = cvscale,
                 scale.gam = scale.gam))
}

CVgam01=CVgam(V7 ~ V1 + V2 + V3 +V4+V5+V6,   data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)

CVgam01=CVgam(V7 ~ poly(V1, 2)  + poly(V2, 3)  + V3   +V4+as.factor(V5>1)*V2
  + poly(V6, 3) +V1 *as.factor(V4>0)    ,    data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)

lmfit=lm(V7 ~ poly(V1, 2)  + poly(V2, 3)  + V3   +V4+as.factor(V5>1)*V2
  + poly(V6, 3) +V1 *as.factor(V4>0)    ,    data=data) 
  
par(mfrow = c(2,2))
plot(lmfit)

#~~~~~~~~~~~~~~~~~~~~~~~~
library(mgcv)

CVscales =array(rep(0,25*25),dim=c(5,5,5,5))

for (i in 6:10)
for (j  in 6:10)
for (l  in 6:10)
for (v  in 6:10)
{{{{
CVgam01=CVgam(V7 ~ s(V1,k=i) +s(V2,k=j) + s(V3,k=l) +s(V6,k=v) ,   data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)
CVscales[i-5,j-5,l-5,v-5] = CVgam01$cvscale
}}}}

which(CVscales== min(CVscales), arr.ind = TRUE)

CVgam01=CVgam(V7 ~ s(V1,k=10) +s(V2,k=6) + s(V3,k=6) +s(V6,k=9) ,   data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)

CVgam01=CVgam(V7 ~ s(V1,k=10) +s(V2,k=6) + s(V3,k=6) +s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1),   data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)

b <- gam(V7 ~ s(V1,k=10) +s(V2,k=6) + s(V3,k=6) +s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1),
method="REML", data=data, select=TRUE) 
summary(b)

plot(b,pages=1,residuals=TRUE)

#~~~~~~~
CVgam01=CVgam(V7 ~ s(V1,k=10) +s(V2,k=6) + V3+s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1) +s(V1,by=as.factor(V4>0)) +s(V2,by=as.factor(V5>1))    , 
  data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)

CVgam01=CVgam(V7 ~ s(V1,k=10) +s(V2,k=6) + V3+s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1) +s(V2,by=as.factor(V5>1))    , 
  data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)


b <- gam(V7 ~ s(V1,k=10) +s(V2,k=6) + V3+s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1) +s(V2,by=as.factor(V5>1))   ,
method="REML", data=data, select=TRUE) 
summary(b)
plot(b,pages=1,residuals=TRUE)

CVgam01=CVgam(V7 ~ s(V1,k=10) +s(V2,k=6) + V3+s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1) +V2*as.factor(V5>1)    , 
  data=data, nfold = 10, debug.level = 0, method = "GCV.Cp", printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)


b <- gam(V7 ~ s(V1,k=10) +s(V2,k=6) + V3+s(V6,k=9) 
+as.factor(V4>0) + as.factor(V5>1) +V2*as.factor(V5>1)   ,
method="REML", data=data, select=TRUE) 
summary(b)
plot(b,pages=1,residuals=TRUE)


residuals<-b$residuals

library(ggplot2)
p1<- ggplot(data, aes( fit,  residuals)) +
  geom_point()
p1<-p1 + geom_smooth(method = "loess")      
p2<- ggplot(data, aes( fit,  residuals^2)) +
  geom_point()
p2<-p2 + geom_smooth(method = "loess")      

p3 <- ggplot(data, aes(sample=residuals))
p3<-p3 + stat_qq() + stat_qq_line()
grid.arrange(p1, p2, p3, ncol = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`











