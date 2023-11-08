### The function to perform dynamic RSF-BE with 1000 bootstraps using parallel processing 
rsfbe_parallel <- function(seed) {
  # If the only remaining columns are covariates, exit the loop
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  # resampling with replacement
  ml.boot <- ml[boot.indices, ]
  ml.temp <- ml.boot
  
  rsf.err <- rep(0,count)				# list to save error rates
  rm.var <- rep(0,count)			# list to save removed variables
  
  vrb.list <- as.list(rep(NA, count))
  for(k in 1 : count){
    if(all(names(ml.temp) %in% c("os5.months","os5.status","age","stage"))){break}
    # if(all(names(ml.temp) %in% c("pfs5.months","pfs5.months","age","stage"))){break}
    
    # grow the RSF
    set.seed(seed)
    rsf.out <- rfsrc(Surv(time = os5.months, event = os5.status) ~ ., data = ml.temp, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
    # rsf.out <- rfsrc(Surv(time = pfs5.months, event = pfs5.status) ~ ., data = ml.temp, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
    
    # save error rate of current RSF
    rsf.err[k] <- rsf.out$err.rate[ntree]
    # get list of variables ordered by their minimal depth values
    v.max <- max.subtree(rsf.out)
    d <- sort(round(v.max$order[, 1], 3))
    
    vrb.list[[k]]<-colnames(ml.temp)[!colnames(ml.temp)%in% c("os5.months","os5.status","age","stage")]
    # vrb.list[[k]]<-colnames(ml.temp)[!colnames(ml.temp)%in% c("pfs5.months","pfs5.status","age","stage")]
    
    # dynamic backward elimination
    nToRemove <- ifelse(length(vrb.list[[k]]) > 200, 100,
                        ifelse(length(vrb.list[[k]]) > 50, 10, 1))
    # get the features with worst minimal depth values
    outvars <- NULL
    m <- 0

    while(length(outvars) < nToRemove && m < length(d)){
      outvar <- names(d[length(d)-m])
      # If outvar exists and is not a covariate, then add to outvars list
      if(length(outvar) > 0 && !(outvar %in% covar)){
        outvars <- c(outvars, outvar)
      }
      m <- m + 1
    }
    # save name of removed variable in a list
    rm.var[k]<-paste(outvars, collapse=", ")
    ml.temp<-ml.temp[,!(names(ml.temp) %in% outvars)]
    print(colnames(ml.temp))
    # remove old RSF to ensure free working memory 
    rsf.out <- NULL
  }
  
  output.os.f<-as.data.frame(cbind(rm.var,rsf.err))
  
  min_err<-min(output.os.f$rsf.err[which(output.os.f$rsf.err>0)])
  min_err_varlist.ind<-which(output.os.f$rsf.err==min_err)
  varb<-vrb.list[min_err_varlist.ind][[1]][-c(1:2)] #variable list with highest C index
  # print(min_err)
  tr.cind<-1-as.numeric(min_err)
  output <- list("tr.cind" = tr.cind, "varb.store.m" = varb)
  return(output)
}


### The functions to generate data for drawing partial dependence plots
### based on dynamic RSF-BE with 1000 bootstraps using parallel processing 
rsfbe.plot.os5 <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = os5.months, event = os5.status) ~., 
                   data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=36,smooth.lines = F)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=60,smooth.lines = F)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
rsfbe.plot.os <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = os.months, event = os.status) ~., 
                   data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=time1,smooth.lines = smooth.lines,npts =npts)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=time2,smooth.lines = smooth.lines,npts =npts)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
rsfbe.plot.rfs <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = rfs.months, event = rfs.status) ~., 
                   data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=36,smooth.lines = F)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=60,smooth.lines = F)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
rsfbe.plot.pfs5 <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = pfs5.months, event = pfs5.status) ~., data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=36,smooth.lines = F)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=60,smooth.lines = F)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}