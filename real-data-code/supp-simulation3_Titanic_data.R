########################################################################
##########Code for Simulation(3) in SM: Titanic data with d=7###########
########################################################################
library(gtools)
library(ggplot2) 
library(MASS)
library(stringr)
library(readxl)

#=========================Train model===========================
# The Titanic data downloaded from https://www.kaggle.com/c/titanic/data records some features of passengers and whether or not they survived the sinking of the Titanic. 
# We choose 7 features (Pclass, Sex, Age, SibSp, Parch, Fare and Embarked) and keep 1043 inputs with complete records in the file ¡®dataset_simulated.xlsx¡¯.
titanic <- read_excel("dataset_simulated.xlsx")
titanic<-as.data.frame(titanic)
ntitanic<-nrow(titanic)
rageage<-max(titanic[,3])-min(titanic[,3])
ragefare<-max(titanic[,6])-min(titanic[,6])
logismodel<-glm(Survived~Pclass+Sex+Age+SibSp+Parch+Fare+Embarked, family=binomial(link = logit),data=titanic)
summary(logismodel)
ybar<-mean(logismodel$fitted.values)

#=========================Functions=============================
#Baseline Shapley
val.BS<-function(sets,xt){
  value<-0
  for(i in 1:ntitanic){
    xnew<-vector(length=d)
    for(j in 1:d){
      if(j %in% sets){xnew[j]<-xt[j]}
      else xnew[j]<-titanic[i,j]
    }
    xnew<-as.vector(unlist(xnew))
    value<-value+predict(logismodel,newdata =data.frame(Pclass=xnew[1],Sex=xnew[2],Age=xnew[3],SibSp=xnew[4],Parch=xnew[5],Fare=xnew[6],Embarked=xnew[7]),type='response')
  }
  value<-value/ntitanic
  return(value)
}

#conditional expectation Shapley
val.CES<-function(sets,xt){
  value<-0; ctid<-c(); nsets<-length(sets);
  for(i in 1:ntitanic){
    isequal<-1
    for(j in 1:nsets){
      if(xt[sets[j]]!=titanic[i,sets[j]]) {isequal<-0;break;}
    }
    if(isequal==1) ctid<-c(ctid,i);
  }
  value<-sum(logismodel$fitted.values[ctid])/length(ctid)
  return(value)
}

#cohort Shapley
val.CS<-function(sets,xt){
  value<-0; ctid<-c();
  for(i in 1:ntitanic){
    isequal<-1
    for(j in intersect(c(1,2,4,5,7),sets)){
      if(xt[j]!=titanic[i,j]) {isequal<-0;break;}
    }
    if((isequal==1)&(3%in%sets)){
      if(abs(xt[3]-titanic[i,3])>( 0.1* rageage)) isequal<-0;
    }
    if((isequal==1)&(6%in%sets)){
      if(abs(xt[6]-titanic[i,6])>( 0.1* ragefare)) isequal<-0;
    }
    if(isequal==1) ctid<-c(ctid,i);
  }
  value<-sum(logismodel$fitted.values[ctid])/length(ctid)
  return(value)
}


real.fash<-function(d,val,xt,ybar){
  perms<-permutations(d,d); n<-factorial(d);
  sh<-rep(0,d);
  for(l in 1:n){
    perml<-perms[l,]; preC<-ybar;
    for(i in 1:d){
      if(i==1){
        delta<-val(perml[1],xt)-preC; sh[perml[1]]<-sh[perml[1]]+delta; preC<-delta+preC;
      }
      else{
        delta<-val(perml[1:i],xt)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}
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

est.fashsrs<-function(d,n,val,xt,ybar){
  sh<-rep(0,d)
  for(l in 1:n){
    perml<-sample(d); preC<-ybar;
    for(i in 1:d){
      if(i==1){
        delta<-val(perml[1],xt)-preC; sh[perml[1]]<-sh[perml[1]]+delta; preC<-delta+preC;
      }
      else{
        delta<-val(perml[1:i],xt)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

est.fashls<-function(d,n,val,xt,ybar){
  sh<-rep(0,d); k<-n/d;
  for(t in 1:k){
    lst<-onels(d);
    for(l in 1:d){
      perml<-lst[l,]; preC<-ybar;
      for(i in 1:d){
        if(i==1){
          delta<-val(perml[1],xt)-preC; sh[perml[1]]<-sh[perml[1]]+delta; preC<-delta+preC;
        }
        else{
          delta<-val(perml[1:i],xt)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
        }
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}


est.fashcoa<-function(d,n,val,xt,ybar){
  sh<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  for(t in 1:k){
    coat<-onecoa(d);
    for(l in 1:m){
      perml<-coat[l,]; preC<-ybar;
      for(i in 1:d){
        if(i==1){
          delta<-val(perml[1],xt)-preC; sh[perml[1]]<-sh[perml[1]]+delta; preC<-delta+preC;
        }
        else{
          delta<-val(perml[1:i],xt)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
        }
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

strict.est.fash<-function(d,n,samplemethod,val,xt,ybar){ 
  if(samplemethod=='SRS'){return(est.fashsrs(d,n,val,xt,ybar))}
  else if(samplemethod=='LS'){return(est.fashls(d,n,val,xt,ybar))}
  else {return(est.fashcoa(d,n,val,xt,ybar))}
}

strictpartial.fa<-function(real,inputid,d,n,samplemethod,val,rt,ybar){
  nc<-d+3; nn<-length(n); nall<-nn*rt; nother<-nn; xt<-titanic[inputid,]
  resultall<-matrix(0,nall,nc);  conave<-matrix(0,nother,nc);conmin<-matrix(0,nother,nc);
  conmax<-matrix(0,nother,nc); conmedian<-matrix(0,nother,nc); varsh<-matrix(0,nother,d+1)
  for(i in 1:nn){
    result<-matrix(0,rt,nc)
    for(j in 1:rt){
      told<-Sys.time()
      sh<-strict.est.fash(d,n[i],samplemethod,val,xt,ybar); time<-Sys.time()-told; mse<-sum((sh-real)^2)/d;
      result[j,]<-c(sh,mse,time,n[i])
    }
    varsh[i,]<-c(apply(result[,1:d],2,var),n[i]);
    conave[i,]<-apply(result,2,mean);
    conmedian[i,]<-apply(result,2,median);
    conmin[i,]<-apply(result,2,min);
    conmax[i,]<-apply(result,2,max);
    resultall[((i-1)*rt+1):(i*rt),]<-result
  }
  colnames(resultall)<-c(1:d,'mse','time','size')
  rownames(conave)<-n; rownames(conmedian)<-n;rownames(conmin)<-n;rownames(conmax)<-n; rownames(varsh)<-n;
  colnames(conave)<-c(1:d,'mse','time','size');colnames(conmedian)<-c(1:d,'mse','time','size');colnames(conmin)<-c(1:d,'mse','time','size');colnames(conmax)<-c(1:d,'mse','time','size');colnames(varsh)<-c(1:d,'size');
  return(list(realvalue=real, result=data.frame(resultall,method=samplemethod,id=inputid), ave=data.frame(conave,method=samplemethod), median=data.frame(conmedian,method=samplemethod), min=data.frame(conmin,method=samplemethod), max=data.frame(conmax,method=samplemethod),var=varsh,method=data.frame(samplemethod)))
}


#============================Simulation==========================
d<-7;n1<-c(c(1:(5*(d-1)))*d);n2<-c(c(1:5)*d*(d-1)); 
#Baseline Shapley
t0<-Sys.time()
real_sh_bs_first<-real.fash(d,val.BS,titanic[1,],ybar)
Sys.time()-t0
#Time difference of 11.9524 hours
allresult_bs_first_srs<-strictpartial.fa(real_sh_bs_first,1,d,n2,'SRS',val.BS,20,ybar)
allresult_bs_first_ls<-strictpartial.fa(real_sh_bs_first,1,d,n2,'LS',val.BS,20,ybar)
allresult_bs_first_coa<-strictpartial.fa(real_sh_bs_first,1,d,n2,'COA',val.BS,20,ybar)


#Conditional expectation Shapley
t0<-Sys.time()
real_sh_ces_first<-real.fash(d,val.CES,titanic[1,],ybar)
t_real_sh_ces_first<-Sys.time()-t0
allresult_ces_first_srs<-strictpartial.fa(real_sh_ces_first,1,d,n2,'SRS',val.CES,20,ybar)
allresult_ces_first_ls<-strictpartial.fa(real_sh_ces_first,1,d,n2,'LS',val.CES,20,ybar)
allresult_ces_first_coa<-strictpartial.fa(real_sh_ces_first,1,d,n2,'COA',val.CES,20,ybar)


#Cohort Shapley
t0<-Sys.time()
real_sh_cs_first<-real.fash(d,val.CS,titanic[1,],ybar)
t_real_sh_cs_first<-Sys.time()-t0
allresult_cs_first_srs<-strictpartial.fa(real_sh_cs_first,1,d,n2,'SRS',val.CS,20,ybar)
allresult_cs_first_ls<-strictpartial.fa(real_sh_cs_first,1,d,n2,'LS',val.CS,20,ybar)
allresult_cs_first_coa<-strictpartial.fa(real_sh_cs_first,1,d,n2,'COA',val.CS,20,ybar)


#==========================================Table S3==========================================
real_titanic<-rbind(real_sh_bs_first,real_sh_ces_first,real_sh_cs_first)


#==========================================Figure S4==========================================
#BS
allresult_bs<-rbind(allresult_bs_first_srs$result,allresult_bs_first_ls$result,allresult_bs_first_coa$result)
allresult_bs$size<-as.character(allresult_bs$size)
allresult_bs$size[allresult_bs$size==42]<-'1m(m-1)'
allresult_bs$size[allresult_bs$size==84]<-'2m(m-1)'
allresult_bs$size[allresult_bs$size==126]<-'3m(m-1)'
allresult_bs$size[allresult_bs$size==168]<-'4m(m-1)'
allresult_bs$size[allresult_bs$size==210]<-'5m(m-1)'
allresult_bs$size<-as.character(allresult_bs$size)
#(a)
allresult_bs_n42<-allresult_bs[allresult_bs$size=='1m(m-1)',]
allest_bs_n42<-c(); nc<-ncol(allresult_bs_n42);
for(i in 1:d){
  t<-allresult_bs_n42[,c(i,(d+1):nc)]; names(t)[1]<-'estimate';
  allest_bs_n42<-rbind(allest_bs_n42,data.frame(t,variable=paste('x',i,sep='')))
}
tnames<-c()
for(i in 1:d){tnames<-c(tnames,paste('x',i,sep=''))}
ggplot(data=allest_bs_n42,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=c('Pclass','Sex','Age','SibSp','Parch','Fare','Embarked'))+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.93,0.28))
#(b)
ggplot(data=allresult_bs,aes(x=size, y=d*mse, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(labels=n2)+
  labs(x='Sample size',y='Square loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))



#CES
allresult_ces<-rbind(allresult_ces_first_srs$result,allresult_ces_first_ls$result,allresult_ces_first_coa$result)
allresult_ces$size<-as.character(allresult_ces$size)
allresult_ces$size[allresult_ces$size==42]<-'1m(m-1)'
allresult_ces$size[allresult_ces$size==84]<-'2m(m-1)'
allresult_ces$size[allresult_ces$size==126]<-'3m(m-1)'
allresult_ces$size[allresult_ces$size==168]<-'4m(m-1)'
allresult_ces$size[allresult_ces$size==210]<-'5m(m-1)'
allresult_ces$size<-as.character(allresult_ces$size)
#(c)
allresult_ces_n42<-allresult_ces[allresult_ces$size=='1m(m-1)',]
allest_ces_n42<-c(); nc<-ncol(allresult_ces_n42);
for(i in 1:d){
  t<-allresult_ces_n42[,c(i,(d+1):nc)]; names(t)[1]<-'estimate';
  allest_ces_n42<-rbind(allest_ces_n42,data.frame(t,variable=paste('x',i,sep='')))
}
tnames<-c()
for(i in 1:d){tnames<-c(tnames,paste('x',i,sep=''))}
ggplot(data=allest_ces_n42,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=c('Pclass','Sex','Age','SibSp','Parch','Fare','Embarked'))+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.93,0.28))
#(d)
ggplot(data=allresult_ces,aes(x=size, y=d*mse, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(labels=n2)+
  labs(x='Sample size',y='Square loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))



#CS
allresult_cs<-rbind(allresult_cs_first_srs$result,allresult_cs_first_ls$result,allresult_cs_first_coa$result)
allresult_cs$size<-as.character(allresult_cs$size)
allresult_cs$size[allresult_cs$size==42]<-'1m(m-1)'
allresult_cs$size[allresult_cs$size==84]<-'2m(m-1)'
allresult_cs$size[allresult_cs$size==126]<-'3m(m-1)'
allresult_cs$size[allresult_cs$size==168]<-'4m(m-1)'
allresult_cs$size[allresult_cs$size==210]<-'5m(m-1)'
allresult_cs$size<-as.character(allresult_cs$size)
#(e)
allresult_cs_n42<-allresult_cs[allresult_cs$size=='1m(m-1)',]
allest_cs_n42<-c(); nc<-ncol(allresult_cs_n42);
for(i in 1:d){
  t<-allresult_cs_n42[,c(i,(d+1):nc)]; names(t)[1]<-'estimate';
  allest_cs_n42<-rbind(allest_cs_n42,data.frame(t,variable=paste('x',i,sep='')))
}
tnames<-c()
for(i in 1:d){tnames<-c(tnames,paste('x',i,sep=''))}
ggplot(data=allest_cs_n42,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=c('Pclass','Sex','Age','SibSp','Parch','Fare','Embarked'))+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.93,0.28))
#(f)
ggplot(data=allresult_cs,aes(x=size, y=d*mse, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(labels=n2)+
  labs(x='Sample size',y='Square loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))
