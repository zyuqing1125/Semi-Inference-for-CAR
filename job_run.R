source("/home/simu_updated_empirical_bootstrap_s2_conditionalX2_gexp.R")

library(mgcv)
library(MASS)
library(splines2)



args <- commandArgs(TRUE)
seed <- as.numeric(args[[1]])

res1<- delta(seed)

save(res1,file=paste("/home/simu_updated_empirical_bootstrap_s2_conditionalX2_gexp/seed",seed,".RData",sep=""))