################################################################
###Code for Simulation(1) in SM: Electoral college withd=51#####
################################################################
library(openxlsx) 
library(ggplot2) 
library(MASS)
library(foreach)
library(doParallel)
library(stringr)
####===================================Data=========================================
# '2020voting.xlsx' can be created according to https://www.archives.gov/electoral-college/2020, with two columns corresponding to names and votes of 51 states respectively.
#full data
data_2020voting <- read.xlsx("2020voting.xlsx")
# View(data_2020voting)

data_sspi<-as.vector(data_2020voting[,2])


####================================Functions========================================
##---------------------1.value function-----------------------
val.sspi<-function(set,votes){
  if((sum(votes[set]))>269){value<-1;}
  else{value<-0;}
  return(value)
}

##---------------------2.estimate of shapley---------------------
#---------------------------2.1 SRS---------------------------
est.shsrs.sspi<-function(d,n,val,votes){
  sh<-rep(0,d)
  for(l in 1:n){
    perml<-sample(d); preC<-0;
    for(i in 1:d){
      if(preC==1){break;}
      delta<-val(perml[1:i],votes)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
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

est.shls.sspi<-function(d,n,val,votes){
  sh<-rep(0,d); k<-n/d;
  for(t in 1:k){
    lst<-onels(d);
    for(l in 1:d){
      perml<-lst[l,]; preC<-0;
      for(i in 1:d){
        if(preC==1){break;}
        delta<-val(perml[1:i],votes)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
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

est.shcoa.sspi<-function(d,n,val,votes){
  dc<-d+2
  sh<-rep(0,d); k<-n/dc/(dc-1);m<-dc*(dc-1);
  for(t in 1:k){
    coat<-onecoa(dc);
    for(l in 1:m){
      perml<-setdiff(coat[l,],c(52,53)); preC<-0;
      for(i in 1:d){
        if(preC==1){break;}
        delta<-val(perml[1:i],votes)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
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

est.shssrs.sspi<-function(d,n,val,votes){
  sh<-rep(0,d); groupsize<-n/d
  srsperm<-matrix(0,nrow=n,ncol=d)
  for(i in 1:n){
    srsperm[i,]<-sample(d)
  }
  for(j in 1:d){
    jperm<-structed.perm(srsperm,j,d)
    for(i in (groupsize+1):n){
      loc<-ceiling(i/groupsize)
      sh[j]<-sh[j]+val(jperm[i,1:loc],votes)-val(jperm[i,1:(loc-1)],votes)
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}


#---------------------------2.5 Summary----------------------------
strict.est.sh.sspi<-function(d,n,samplemethod,val,votes){ 
  if(samplemethod=='SRS'){return(est.shsrs.sspi(d,n,val,votes))}
  else if(samplemethod=='LS'){return(est.shls.sspi(d,n,val,votes))}
  else if(samplemethod=='COA') {return(est.shcoa.sspi(d,n,val,votes))}
  else {return(est.shssrs.sspi(d,n,val,votes))}
}

#---------------------------2.6 Parallel computing---------------------
strict.sspi.bx<-function(d,n,samplemethod,val,votes,rt,cores){
  nc<-d+3; nn<-length(n); nall<-nn*rt; nother<-nn;
  resultall<-matrix(0,nall,nc);  
  for(i in 1:nn){
    cores<-min(detectCores(logical=F),cores)
    cl <- makeCluster(cores)
    registerDoParallel(cl,cores=sepnum)
    result<-foreach(j=1:rt, .combine= "rbind",.export =c("strict.est.sh.sspi","est.shsrs.sspi","est.shssrs.sspi","est.shls.sspi","est.shcoa.sspi","onels","onecoa","structed.perm")) %dopar%
      {
        t0<-Sys.time()
        sh<-strict.est.sh.sspi(d,n[i],samplemethod,val,votes); 
        s <- c(n[i],0,as.numeric(Sys.time()-t0,units="secs"),sh)
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
  colnames(resultall)<-c('size','rt','time',1:d)
  return(data.frame(resultall,method=samplemethod))
}


####================================Simulations================================
d_sspi<-51
n_sspi_coa<-c(1:3)*53*52
n_sspi<-c(51,51*3,51*6,seq(510,ceiling(n_sspi_coa[3]/d_sspi)*d_sspi,9*d_sspi))
n_nsspi<-length(n_sspi)
n_nsspi_coa<-length(n_sspi_coa)


##------------------1.SRS------------------------
n_mz<-10;n_zs<-1000/n_mz;

for(i in 1:n_nsspi){
  assign(paste("result_sspi_srs_", i, sep=""), c())
}
time_sspi_srs<-matrix(0,nrow=n_zs,ncol=n_nsspi)
for(ii in 1:n_zs){
  for(jj in 1:n_nsspi){
    print(paste(paste("srs-round:",ii,sep=""),paste("; m:",n_sspi[jj],sep=""),sep=""))
    time0_srs<-Sys.time()
    t_srs<-strict.sspi.bx(d_sspi,n_sspi[jj],'SRS',val.sspi,data_sspi,n_mz,10)
    t_srs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_sspi_srs_", jj, sep=""), rbind(get(paste("result_sspi_srs_", jj, sep="")),t_srs))
    time_sspi_srs[ii,jj]<-as.numeric(Sys.time()-time0_srs,units="secs")
    tsave_sspi_srs<-get(paste("result_sspi_srs_", jj, sep=""))
    save(tsave_sspi_srs,file=eval(paste(paste("result_sspi_srs_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_sspi_srs[ii,jj],sep=""))
  }
}
result_sspi_srs<-c()
for(i in 1:n_nsspi){
  result_sspi_srs<-rbind(result_sspi_srs,get(paste("result_sspi_srs_", i, sep="")))
}


##------------------2.strrs------------------------
n_mz<-10;n_zs<-1000/n_mz;

for(i in 1:n_nsspi){
  assign(paste("result_sspi_strrs_", i, sep=""), c())
}
time_sspi_strrs<-matrix(0,nrow=n_zs,ncol=n_nsspi)
for(ii in 1:n_zs){
  for(jj in 1:n_nsspi){
    print(paste(paste("strrs-round:",ii,sep=""),paste("; m:",n_sspi[jj],sep=""),sep=""))
    time0_strrs<-Sys.time()
    t_strrs<-strict.sspi.bx(d_sspi,n_sspi[jj],'StrRS',val.sspi,data_sspi,n_mz,10)
    t_strrs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_sspi_strrs_", jj, sep=""), rbind(get(paste("result_sspi_strrs_", jj, sep="")),t_strrs))
    time_sspi_strrs[ii,jj]<-as.numeric(Sys.time()-time0_strrs,units="secs")
    tsave_sspi_strrs<-get(paste("result_sspi_strrs_", jj, sep=""))
    save(tsave_sspi_strrs,file=eval(paste(paste("result_sspi_strrs_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_sspi_strrs[ii,jj],sep=""))
  }
}
result_sspi_strrs<-c()
for(i in 1:n_nsspi){
  result_sspi_strrs<-rbind(result_sspi_strrs,get(paste("result_sspi_strrs_", i, sep="")))
}

##------------------3.LS------------------------
n_mz<-10;n_zs<-1000/n_mz;

for(i in 1:n_nsspi){
  assign(paste("result_sspi_ls_", i, sep=""), c())
}
time_sspi_ls<-matrix(0,nrow=n_zs,ncol=n_nsspi)
for(ii in 1:n_zs){
  for(jj in 1:n_nsspi){
    print(paste(paste("ls-round:",ii,sep=""),paste("; m:",n_sspi[jj],sep=""),sep=""))
    time0_ls<-Sys.time()
    t_ls<-strict.sspi.bx(d_sspi,n_sspi[jj],'LS',val.sspi,data_sspi,n_mz,10)
    t_ls[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_sspi_ls_", jj, sep=""), rbind(get(paste("result_sspi_ls_", jj, sep="")),t_ls))
    time_sspi_ls[ii,jj]<-as.numeric(Sys.time()-time0_ls,units="secs")
    tsave_sspi_ls<-get(paste("result_sspi_ls_", jj, sep=""))
    save(tsave_sspi_ls,file=eval(paste(paste("result_sspi_ls_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_sspi_ls[ii,jj],sep=""))
  }
}
result_sspi_ls<-c()
for(i in 1:n_nsspi){
  result_sspi_ls<-rbind(result_sspi_ls,get(paste("result_sspi_ls_", i, sep="")))
}


##------------------4.COA------------------------
n_mz<-10;n_zs<-1000/n_mz;

for(i in 1:n_nsspi_coa){
  assign(paste("result_sspi_coa_", i, sep=""), c())
}
time_sspi_coa<-matrix(0,nrow=n_zs,ncol=n_nsspi_coa)
for(ii in 1:n_zs){
  for(jj in 1:n_nsspi_coa){
    print(paste(paste("coa-round:",ii,sep=""),paste("; m:",n_sspi_coa[jj],sep=""),sep=""))
    time0_coa<-Sys.time()
    t_coa<-strict.sspi.bx(d_sspi,n_sspi_coa[jj],'COA',val.sspi,data_sspi,n_mz,10)
    t_coa[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_sspi_coa_", jj, sep=""), rbind(get(paste("result_sspi_coa_", jj, sep="")),t_coa))
    time_sspi_coa[ii,jj]<-as.numeric(Sys.time()-time0_coa,units="secs")
    tsave_sspi_coa<-get(paste("result_sspi_coa_", jj, sep=""))
    save(tsave_sspi_coa,file=eval(paste(paste("result_sspi_coa_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_sspi_coa[ii,jj],sep=""))
  }
}
result_sspi_coa<-c()
for(i in 1:n_nsspi_coa){
  result_sspi_coa<-rbind(result_sspi_coa,get(paste("result_sspi_coa_", i, sep="")))
}



####================================Analysis================================
#-----------------------------0. pseudo real value-----------------------------------

pseudo_sspi<-rep(0,d_sspi)
for(i in 1:n_nsspi){
  for(j in 1:1000){
    pseudo_sspi<-pseudo_sspi+result_sspi_srs[(i-1)*1000+j,3+1:d_sspi]*n_sspi[i]
  }
}
pseudo_sspi<-pseudo_sspi/sum(n_sspi)/1000
pseudo_sspi

#Table S1
sort_sspi<-matrix(nrow=d_sspi,ncol=4)
sort_sspi[,1]<-c(1:d_sspi)
sort_sspi[,2]<-data_2020voting[sort(as.numeric(pseudo_sspi),decreasing = TRUE,index.return=TRUE)$ix,1]
sort_sspi[,3]<-data_2020voting[sort(as.numeric(data_2020voting[,2]),decreasing = TRUE,index.return=TRUE)$ix,2]
sort_sspi[,4]<-as.numeric(round(sort(pseudo_sspi,decreasing = TRUE),4))
colnames(sort_sspi)<-c('Rank','State','Votes','Sh')

#---------------------1. summarize------------------------
cols_sspi<- c("LS" = "#00CC66", "SRS" = "#FF6666", "StrRS" = "#AB82FF", "COA" = "#6699FF")
shape_sspi<-c("LS" = 15, "SRS" =16, "StrRS" = 17, "COA" =18)

allresult_sspi<-rbind(result_sspi_srs,result_sspi_strrs,result_sspi_ls,result_sspi_coa)

size_sspi<-as.factor(c('51','153','306','510','969','1428','1887','2346','2805 (2756)','3264','3723','4182','4641','5100','5559 (5512)','6018','6477',
                       '6936','7395','7854','8313 (8268)'))



allresult_sspi[allresult_sspi$size==2805,]$size<-'2805 (2756)'
allresult_sspi[allresult_sspi$size==2756,]$size<-'2805 (2756)'
allresult_sspi[allresult_sspi$size==5559,]$size<-'5559 (5512)'
allresult_sspi[allresult_sspi$size==5512,]$size<-'5559 (5512)'
allresult_sspi[allresult_sspi$size==8313,]$size<-'8313 (8268)'
allresult_sspi[allresult_sspi$size==8268,]$size<-'8313 (8268)'


#---------------------2. square loss: Figure S1 (b) and (c)------------------------

allmse_sspi<-c()
for(i in 1:nrow(allresult_sspi)){
  allmse_sspi<-c(allmse_sspi,sum((allresult_sspi[i,3+c(1:51)]-pseudo_sspi)^2))
}
allmse_sspi<-data.frame(mse=allmse_sspi,size=allresult_sspi$size,rt=allresult_sspi$rt,method=allresult_sspi$method)
allmse_sspi$method<-as.factor(allmse_sspi$method)
allmse_sspi$size<-as.factor(allmse_sspi$size)

allmse_sspi_partial_1<-allmse_sspi[allmse_sspi$size%in%c('51','153','306','510','969','1428','1887','2346','2805 (2756)'),]
allmse_sspi_partial_1_clear<-allmse_sspi_partial_1[-which(allmse_sspi_partial_1$mse>0.05),]
allmse_sspi_partial_2<-allmse_sspi[allmse_sspi$size%in%c('2805 (2756)','3264','3723','4182','4641','5100','5559 (5512)','6018','6477',
                                                         '6936','7395','7854','8313 (8268)'),]
#Figure S1.(b)
ggplot(data=allmse_sspi_partial_1_clear,aes(x=factor(size,levels=size_sspi), y=mse, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_sspi,labels = function(x) str_wrap(x, width = 1)) +
  scale_color_manual(values=cols_sspi)+scale_fill_manual(values=cols_sspi)+
  labs(x='Sample size',y='Square loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 12))

#Figure S1.(c)
ggplot(data=allmse_sspi_partial_2,aes(x=factor(size,levels=size_sspi), y=mse, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_sspi,labels = function(x) str_wrap(x, width = 1)) +
  scale_color_manual(values=cols_sspi)+scale_fill_manual(values=cols_sspi)+
  labs(x='Sample size',y='Square loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 12))

#---------------------3. estimation: Figure S1.(a)------------------------
allresult_sspi_n2805<-allresult_sspi[allresult_sspi$size=='2805 (2756)',]
allest_sspi_n2805<-c(); nc_sspi<-ncol(allresult_sspi_n2805);
for(i in 1:d_sspi){
  t<-allresult_sspi_n2805[,c(3+i,nc_sspi)]; names(t)[1]<-'estimate';
  allest_sspi_n2805<-rbind(allest_sspi_n2805,data.frame(t,variable=data_2020voting[i,1]))
}
tnames<-c()
for(i in 1:d_sspi){tnames<-c(tnames,data_2020voting[i,1])}


# ¡®name.xlsx¡¯ is needed to great Figure S1.(a), whose first column is same as the first column of ¡®2020voting.xlsx¡¯ and the second column shows abbreviations of names in the first column.
sspi_names <- read.xlsx("name.xlsx",colNames=FALSE)
sspi_names_xx<-as.vector(sspi_names[,2])
sspi_names_xx[8]<-"DC"


ggplot(data=allest_sspi_n2805,aes(x=variable, y=estimate, fill=factor(method,level=c('SRS','StrRS','LS','COA')),color=factor(method,level=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels = sspi_names_xx)+
  scale_color_manual(values=cols_sspi)+scale_fill_manual(values=cols_sspi)+
  labs(x='State',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.95,0.8),text = element_text(size = 17),axis.text.x = element_text(size = 10))

