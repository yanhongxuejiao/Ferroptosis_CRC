library(survival)
library(randomForestSRC)
library(ggRandomForests)
library(dplyr)
source("RSFBE_boot.R") # load packed function for running dynamic/tiered RSF-BE with 1000 bootstraps

mdata.all <- read.csv("male GSE39582.csv")
fdata.all <- read.csv("female GSE39582.csv")
ferrgene <-  colnames(mdata.all)[c(12:420)]

mdata.all <-  subset(mdata.all,stage !=0 & !is.na(age))
mdata.all$stage<-ifelse(mdata.all$stage == 1|mdata.all$stage== 2, 'early', 'late')

fdata.all <-  subset(fdata.all,stage !=0 & !is.na(age))
fdata.all$stage<-ifelse(fdata.all$stage == 1|fdata.all$stage== 2, 'early', 'late')

factor.cols <- c("stage","KRAS_mut")
mdata.all[factor.cols] <- lapply(mdata.all[factor.cols], as.factor)  
fdata.all[factor.cols] <- lapply(fdata.all[factor.cols], as.factor)  

mdata.all$os5.months <- ifelse(mdata.all$os.time>=60, 60, mdata.all$os.time)
mdata.all$os5.status <- ifelse(mdata.all$os.time>60, 0, mdata.all$os.status)
fdata.all$os5.months <- ifelse(fdata.all$os.time>=60, 60, fdata.all$os.time)
fdata.all$os5.status <- ifelse(fdata.all$os.time>60, 0, fdata.all$os.status)

mdata.all$pfs5.months <- ifelse(mdata.all$rfs.time>=60, 60, mdata.all$rfs.time)
mdata.all$pfs5.status <- ifelse(mdata.all$rfs.time>60, 0, mdata.all$rfs.status)
fdata.all$pfs5.months <- ifelse(fdata.all$rfs.time>=60, 60, fdata.all$rfs.time)
fdata.all$pfs5.status <- ifelse(fdata.all$rfs.time>60, 0, fdata.all$rfs.status)
mdata <- mdata.all[,c("geo_accession", "sex", "os5.months","os5.status","pfs5.months","pfs5.status","age","stage","KRAS_mut",ferrgene)]
fdata <- fdata.all[,c("geo_accession", "sex","os5.months","os5.status","pfs5.months","pfs5.status","age","stage","KRAS_mut",ferrgene)]

data.all<-rbind(mdata,fdata)

# Dynamic RSF-BE with 1000 bootstraps using parallel processing
## An example of applying RSF-BE on males in GSE39582
covar <- c("age","stage") # covariates enforced to be included in the model
count <- length(ferrgene)
B<-1000 # 1000 bootstraps
ntree <- 300				
library(parallel)
mdata.os <- subset(mdata.all, !is.na(os5.months) & os5.months>0& !is.na(os5.status)&os5.status!="")
ml<-mdata.os[,c("os5.months","os5.status","age","stage",ferrgene)]
nsplit <- round(nrow(ml)/15)	# number of node splits
nodesize <- round(nrow(ml)/length(ml$os5.status[ml$os5.status==1])/2) # node size

num_cores <- detectCores() - 1  
cl <- makeCluster(num_cores)
clusterEvalQ(cl, library(randomForestSRC))
clusterExport(cl, list("ml", "covar", "nsplit", "nodesize", "ntree","count"))

sleep_for_a_minute <- function() { Sys.sleep(60) }
start_time <- Sys.time()

results <- parLapply(cl, 1:B, rsfbe_parallel)

end_time <- Sys.time()
print(end_time - start_time)

tr.cind <- unlist(lapply(results, `[[`, 1))
varb.store.m <- unlist(lapply(results, `[[`, 2))
stopCluster(cl)

freq.geo.m <- as.data.frame(table(varb.store.m))
colnames(freq.geo.m)[1] <-"Gene"

write.csv(freq.geo.m,"freq.GSE39582.male.csv", row.names = F)
write.csv(data.frame(c=tr.cind),"cindex.GSE39582.male.csv", row.names = F)



