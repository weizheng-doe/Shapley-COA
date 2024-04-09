
#####################################################################
#####################Code for Airport game with d=101################
#####################################################################
library(ggplot2) 
library(MASS)
library(foreach) 
library(doParallel)
library(stringr)

####================================Setting========================================
##------------------1. weights----------------------------
#d=1.1
weights_air<-rep(c(1:10),c(8,10,7,13,12,11,10,15,10,5))

##------------------2.true value----------------------------
realsh_air<-rep(0,101)
realsh_air[1]<-1/(101+1-1)
for(i in 2:101){
  realsh_air[i]<-realsh_air[i-1]+(weights_air[i]-weights_air[i-1])/(101+1-i)
}
realsh_air

####================================Functions========================================
##---------------------1.value function-----------------------
val.air<-function(set,weights_air){
  value<-max(weights_air[set])
  return(value)
}

##---------------------2.estimate of shapley---------------------
#---------------------------2.1 SRS---------------------------
est.shsrs.air<-function(d,n,val,weights_air){
  sh<-rep(0,d)
  for(l in 1:n){
    perml<-sample(d); preC<-0;
    for(i in 1:d){
      delta<-val(perml[1:i],weights_air)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

#---------------------------2.2 LS----------------------------
onels<-function(d){
  sq <- matrix(nrow=d, ncol=d) 
  sq[1,]<-c(1,sample(d-1)+1);
  for (r in 2:d) {
    sq[r,]<-c(sq[r-1,-1],sq[r-1,1])
  }
  pcol<-sample(d);prow<-c(1,sample(c(2:d)));
  sq<-sq[,pcol]; sq<-sq[prow,];
  return(sq)
}

est.shls.air<-function(d,n,val,weights_air){
  sh<-rep(0,d); k<-n/d;
  for(t in 1:k){
    lst<-onels(d);
    for(l in 1:d){
      perml<-lst[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i],weights_air)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

#---------------------------2.3 COA---------------------------
onecoa<-function(d){
  firstline<-c(0,sample(d-1)); cz<-matrix(0,d-1,d); D<-matrix(0,d*(d-1),d)
  for(i in 1:(d-1)){
    cz[i,]<-(firstline*i)%%d
  }
  for(j in 0:(d-1)){
    D[(j*(d-1)+1):((j+1)*(d-1)),]<-(cz+j)%%d
  }
  D<-D+1
  return(D)
}

est.shcoa.air<-function(d,n,val,weights_air){
  sh<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  for(t in 1:k){
    coat<-onecoa(d);
    for(l in 1:m){
      perml<-coat[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i],weights_air)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

#---------------------------2.4 StrRS-------------------------
structed.perm<-function(permatrix,jcom,d){
  nr<-nrow(permatrix); t<-nr/d
  jpermatrix<-matrix(0,nrow=nr,ncol=d)
  for(i in 1:nr){
    tt<-ceiling(i/t); perm<-permatrix[i,]
    id<-which(perm==jcom)
    tcom<-perm[tt]
    perm[tt]<-jcom
    perm[id]<-tcom
    jpermatrix[i,]<-perm
  }
  return(jpermatrix)
}

est.shssrs.air<-function(d,n,val,weights_air){
  sh<-rep(0,d); groupsize<-n/d
  srsperm<-matrix(0,nrow=n,ncol=d)
  for(i in 1:n){
    srsperm[i,]<-sample(d)
  }
  for(j in 1:d){
    jperm<-structed.perm(srsperm,j,d)
    for(i in 1:groupsize){
      sh[j]<-sh[j]+val(jperm[i,1],weights_air)
    }
    for(i in (groupsize+1):n){
      loc<-ceiling(i/groupsize)
      sh[j]<-sh[j]+val(jperm[i,1:loc],weights_air)-val(jperm[i,1:(loc-1)],weights_air)
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}


#---------------------------2.5 Summary----------------------------
strict.est.sh.air<-function(d,n,samplemethod,val,weights_air){ 
  if(samplemethod=='SRS'){return(est.shsrs.air(d,n,val,weights_air))}
  else if(samplemethod=='LS'){return(est.shls.air(d,n,val,weights_air))}
  else if(samplemethod=='COA') {return(est.shcoa.air(d,n,val,weights_air))}
  else {return(est.shssrs.air(d,n,val,weights_air))}
}

#---------------------------2.6 Parallel computing---------------------
strict.air.bx<-function(real,d,n,samplemethod,val,rt,weights_air,cores){
  nc<-d+4; nn<-length(n); nall<-nn*rt; nother<-nn;
  resultall<-matrix(0,nall,nc);  
  for(i in 1:nn){
    cores<-min(detectCores(logical=F),cores)
    cl <- makeCluster(cores)
    registerDoParallel(cl,cores=sepnum)
    result<-foreach(j=1:rt, .combine= "rbind",.export =c("strict.est.sh.air","est.shsrs.air","est.shssrs.air","est.shls.air","est.shcoa.air","onels","onecoa","structed.perm")) %dopar%
      {
        t0<-Sys.time()
        sh<-strict.est.sh.air(d,n[i],samplemethod,val,weights_air); 
        mse<-sum((sh-real)^2)/d
        s <- c(n[i],0,as.numeric(Sys.time()-t0,units="secs"),mse,sh)
        return(s)
      }
    stopCluster(cl)
    unregister_dopar <- function() {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    result[,2]<-1:rt
    resultall[((i-1)*rt+1):(i*rt),]<-result
  }
  colnames(resultall)<-c('size','rt','time','mse',1:d)
  return(data.frame(resultall,method=samplemethod))
}



####================================Simulations================================
d_air<-101
n_air<-seq(101*100,101*100*10,101*100)
n_nair<-length(n_air)

##----------------------1.srs---------------
n_mz<-20;n_zs<-1000/n_mz;


for(i in 1:n_nair){
  assign(paste("result_airport_srs_", i, sep=""), c())
}
time_airport_srs<-matrix(0,nrow=n_zs,ncol=n_nair)
for(ii in 1:n_zs){
  for(jj in 1:n_nair){
    print(paste(paste("round:",ii,sep=""),paste("; m:",n_air[jj],sep=""),sep=""))
    time0_srs<-Sys.time()
    t_srs<-strict.air.bx(realsh_air,d_air,n_air[jj],'SRS',val.air,n_mz,weights_air,20)
    t_srs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_airport_srs_", jj, sep=""), rbind(get(paste("result_airport_srs_", jj, sep="")),t_srs))
    time_airport_srs[ii,jj]<-as.numeric(Sys.time()-time0_srs,units="secs")
    tsave_airport_srs<-get(paste("result_airport_srs_", jj, sep=""))
    save(tsave_airport_srs,file=eval(paste(paste("result_airport_srs_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_airport_srs[ii,jj],sep=""))
  }
}
result_air_srs<-c()
for(i in 1:n_nair){
  result_air_srs<-rbind(result_air_srs,get(paste("result_airport_srs_", i, sep="")))
}


##----------------------2.strrs---------------
n_mz<-20;n_zs<-1000/n_mz;


for(i in 1:n_nair){
  assign(paste("result_airport_strrs_", i, sep=""), c())
}
time_airport_strrs<-matrix(0,nrow=n_zs,ncol=n_nair)
for(ii in 1:n_zs){
  for(jj in 1:n_nair){
    print(paste(paste("round:",ii,sep=""),paste("; m:",n_air[jj],sep=""),sep=""))
    time0_strrs<-Sys.time()
    t_strrs<-strict.air.bx(realsh_air,d_air,n_air[jj],'StrRS',val.air,n_mz,weights_air,20)
    t_strrs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_airport_strrs_", jj, sep=""), rbind(get(paste("result_airport_strrs_", jj, sep="")),t_strrs))
    time_airport_strrs[ii,jj]<-as.numeric(Sys.time()-time0_strrs,units="secs")
    tsave_airport_strrs<-get(paste("result_airport_strrs_", jj, sep=""))
    save(tsave_airport_strrs,file=eval(paste(paste("result_airport_strrs_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_airport_strrs[ii,jj],sep=""))
  }
}
result_air_strrs<-c()
for(i in 1:n_nair){
  result_air_strrs<-rbind(result_air_strrs,get(paste("result_airport_strrs_", i, sep="")))
}

save(result_air_strrs,file="result/result_airport_strrs.RData")
save(time_airport_strrs,file="result/time_airport_strrs.RData")


##----------------------3.ls---------------
n_mz<-20;n_zs<-1000/n_mz;


for(i in 1:n_nair){
  assign(paste("result_airport_ls_", i, sep=""), c())
}
time_airport_ls<-matrix(0,nrow=n_zs,ncol=n_nair)
for(ii in 1:n_zs){
  for(jj in 1:n_nair){
    print(paste(paste("round:",ii,sep=""),paste("; m:",n_air[jj],sep=""),sep=""))
    time0_ls<-Sys.time()
    t_ls<-strict.air.bx(realsh_air,d_air,n_air[jj],'LS',val.air,n_mz,weights_air,20)
    t_ls[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_airport_ls_", jj, sep=""), rbind(get(paste("result_airport_ls_", jj, sep="")),t_ls))
    time_airport_ls[ii,jj]<-as.numeric(Sys.time()-time0_ls,units="secs")
    tsave_airport_ls<-get(paste("result_airport_ls_", jj, sep=""))
    save(tsave_airport_ls,file=eval(paste(paste("result_airport_ls_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_airport_ls[ii,jj],sep=""))
  }
}
result_air_ls<-c()
for(i in 1:n_nair){
  result_air_ls<-rbind(result_air_ls,get(paste("result_airport_ls_", i, sep="")))
}

save(result_air_ls,file="result/result_airport_ls.RData")
save(time_airport_ls,file="result/time_airport_ls.RData")


##----------------------4.coa---------------
n_mz<-20;n_zs<-1000/n_mz;


for(i in 1:n_nair){
  assign(paste("result_airport_coa_", i, sep=""), c())
}
time_airport_coa<-matrix(0,nrow=n_zs,ncol=n_nair)
for(ii in 1:n_zs){
  for(jj in 1:n_nair){
    print(paste(paste("round:",ii,sep=""),paste("; m:",n_air[jj],sep=""),sep=""))
    time0_coa<-Sys.time()
    t_coa<-strict.air.bx(realsh_air,d_air,n_air[jj],'COA',val.air,n_mz,weights_air,20)
    t_coa[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_airport_coa_", jj, sep=""), rbind(get(paste("result_airport_coa_", jj, sep="")),t_coa))
    time_airport_coa[ii,jj]<-as.numeric(Sys.time()-time0_coa,units="secs")
    tsave_airport_coa<-get(paste("result_airport_coa_", jj, sep=""))
    save(tsave_airport_coa,file=eval(paste(paste("result_airport_coa_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_airport_coa[ii,jj],sep=""))
  }
}
result_air_coa<-c()
for(i in 1:n_nair){
  result_air_coa<-rbind(result_air_coa,get(paste("result_airport_coa_", i, sep="")))
}




####================================Analysis================================
##----------------------1.analysis------------------------
allresult_air<-rbind(result_air_srs,result_air_strrs,result_air_ls,result_air_coa)
size_air<-as.factor(c('10100','20200','30300','40400','50500','60600','70700','80800','90900','101000'))
cols_air<- c("LS" = "#00CC66", "SRS" = "#FF6666", "StrRS" = "#AB82FF", "COA" = "#6699FF")
shape_air<-c("LS" = 15, "SRS" =16, "StrRS" = 17, "COA" =18)

#---------------------2.1 Figure 2(a): estimation------------------------

allresult_air_n10100<-allresult_air[allresult_air$size==10100,]
allest_air_n10100<-c(); nc_air<-ncol(allresult_air_n10100);
for(i in 1:d_air){
  t<-allresult_air_n10100[,c(4+i,nc_air)]; names(t)[1]<-'estimate';
  allest_air_n10100<-rbind(allest_air_n10100,data.frame(t,variable=paste('x',i,sep='')))
}
tnames<-c()
for(i in 1:d_air){tnames<-c(tnames,paste('x',i,sep=''))}
allest_air_n10100_partial<-allest_air_n10100[allest_air_n10100$variable%in%c('x1','x9','x19','x26','x39','x51','x62','x72','x87','x97'),]
allest_air_n10100_partial<-cbind(allest_air_n10100_partial,real=rep(realsh_air[c(1,9,19,26,39,51,62,72,87,97)],each=1000*4))
ggplot(data=allest_air_n10100_partial,aes(x=variable, y=estimate, fill=factor(method,level=c('SRS','StrRS','LS','COA')),color=factor(method,level=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=c('x1','x9','x19','x26','x39','x51','x62','x72','x87','x97'),labels=c('x1','x9','x19','x26','x39','x51','x62','x72','x87','x97'))+
  scale_color_manual(values=cols_air)+scale_fill_manual(values=cols_air)+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.2,0.75),text = element_text(size = 17),axis.text.x = element_text(size = 10))

#---------------------1.2 Figure 2(b): square loss------------------------
allresult_air$size<-as.character(allresult_air$size)
ggplot(data=allresult_air,aes(x=factor(size,level=size_air), y=d_air*mse, fill=factor(method,level=c('SRS','StrRS','LS','COA')),color=factor(method,level=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_air)+
  scale_color_manual(values=cols_air)+scale_fill_manual(values=cols_air)+
  labs(x='Sample size',y='Squared loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.85,0.7),text = element_text(size = 17),axis.text.x = element_text(size = 10))



#---------------------2.3 Table 1: time------------------------
time_air<-matrix(nrow=4,ncol=10)
time_air[1,]<-apply(time_airport_srs,2,mean)
time_air[2,]<-apply(time_airport_strrs,2,mean)
time_air[3,]<-apply(time_airport_ls,2,mean)
time_air[4,]<-apply(time_airport_coa,2,mean)
colnames(time_air)<-as.factor(c('10100','20200','30300','40400','50500','60600','70700','80800','90900','101000'))
rownames(time_air)<-c('SRS','StrRS','LS','COA')
time_air

time_air_ratio<-t(round(apply(time_air,1,function(x) x/time_air[4,]),2))
time_air_ratio

