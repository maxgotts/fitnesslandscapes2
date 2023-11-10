library(futureInterpolate)
# source("~/Desktop/Research/Future Interpolation/futureInterpolate/R/futureInterpolation.R")
library(dplyr)

path0 <- "~/Desktop/Research/Future\ Interpolation/Bobwhite/bobwhite_data/"
nestData <- read.csv(paste0(path0, "OA_KFswqnoboNest.csv"))

nest <- data.frame(birdid=unique(nestData$birdid))
rownames(nest) <- nest$birdid
for (bird in unique(nestData$birdid)) {
  sub <- filter(nestData, birdid==bird)
  nest[bird, "sex"] <- unique(sub$sex)
  nest[bird, "area"] <- unique(sub$area)
  for (var in c("shrub","tree","ag","agid","cg","mg","ng","nggr","ngid","ngpb","ngrpb","swq_d2tree","swq_d2shrub","nestsurv")) {
    nest[bird, var] <- mean(sub[,var])
  }
}

axes <- c("shrub","tree","ag","agid","cg","mg","ng","swq_d2tree","swq_d2shrub")
resp <- "nestsurv"

pca <- perform.PCA(nest, INCLUDE=axes)
nest <- cbind(nest, pca$columns)

tps <- TPS(nest, AXES="PCA", Z=resp, OUTPUT=c("contour","model", "matrix"))

tps$model() %>% summary
tps$matrix() %>% head


