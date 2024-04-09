#########################################################################
###Code for Simulation(2) in SM: Manufacturing system game with d=5######
#########################################################################
library(gtools)
library(ggplot2) 
library(MASS)
library(stringr)
#=========================Functions=============================
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


est.shsareal<-function(d,model,vary,x.sample,condition.sample, No, Ni){
  sh<-rep(0,d); sh2<-rep(0,d);
  perm<-permutations(d,d); n<-factorial(d);
  for(i in 1:n){
    pi<-perm[i,]; preC<-0;
    for(j in 1:d){
      if(j==d){value<-vary}
      else{
        cvary<-c()
        for(l in 1:No){
          y2<-c(); x_nj<-x.sample(pi[-c(1:j)])
          for(h in 1:Ni){
            x_j<-condition.sample(pi[1:j],x_nj,pi[-c(1:j)]); x_all<-c(x_j,x_nj);
            y2<-c(y2,model(x_all[sort(pi,index.return=TRUE)$ix]))
          }
          cvary<-c(cvary,var(y2))
        }
        value<-mean(cvary); 
      }
      delta<-value-preC; preC<-delta+preC;
      sh[pi[j]]<-sh[pi[j]]+delta; sh2[pi[j]]<-sh2[pi[j]]+delta^2;
    }
  }
  sh<-sh/n; sh2<-sh2/n
  return(list(Shapley = sh, SEShapley = sqrt((sh2 - sh^2)/n)))
}


est.shsasrs<-function(d,n,model,vary,x.sample,condition.sample, No, Ni){
  sh<-rep(0,d); sh2<-rep(0,d);
  for(i in 1:n){
    pi<-sample(d); preC<-0;
    for(j in 1:d){
      if(j==d){value<-vary}
      else{
        cvary<-c()
        for(l in 1:No){
          y2<-c(); x_nj<-x.sample(pi[-c(1:j)])
          for(h in 1:Ni){
            x_j<-condition.sample(pi[1:j],x_nj,pi[-c(1:j)]); x_all<-c(x_j,x_nj);
            y2<-c(y2,model(x_all[sort(pi,index.return=TRUE)$ix]))
          }
          cvary<-c(cvary,var(y2))
        }
        value<-mean(cvary); 
      }
      delta<-value-preC; preC<-delta+preC;
      sh[pi[j]]<-sh[pi[j]]+delta; sh2[pi[j]]<-sh2[pi[j]]+delta^2;
    }
  }
  sh<-sh/n; sh2<-sh2/n
  return(list(Shapley = sh, SEShapley = sqrt((sh2 - sh^2)/n)))
}


est.shsals<-function(d,n,model,vary,x.sample,condition.sample, No, Ni){
  sh<-rep(0,d); sh2<-rep(0,d); k<-n/d;
  for(t in 1:k){
    lst<-onels(d);
    for(i in 1:d){
      pi<-lst[i,]; preC<-0;
      for(j in 1:d){
        if(j==d){value<-vary}
        else{
          cvary<-c()
          for(l in 1:No){
            y2<-c(); x_nj<-x.sample(pi[-c(1:j)])
            for(h in 1:Ni){
              x_j<-condition.sample(pi[1:j],x_nj,pi[-c(1:j)]); x_all<-c(x_j,x_nj);
              y2<-c(y2,model(x_all[sort(pi,index.return=TRUE)$ix]))
            }
            cvary<-c(cvary,var(y2))
          }
          value<-mean(cvary); 
        }
        delta<-value-preC; preC<-delta+preC;
        sh[pi[j]]<-sh[pi[j]]+delta; sh2[pi[j]]<-sh2[pi[j]]+delta^2;
      }
    }
  }
  sh<-sh/n; sh2<-sh2/n
  return(list(Shapley = sh, SEShapley = sqrt((sh2 - sh^2)/n)))
}


est.shsacoa<-function(d,n,model,vary,x.sample,condition.sample, No, Ni){
  sh<-rep(0,d); sh2<-rep(0,d); k<-n/d/(d-1);m<-d*(d-1);
  for(t in 1:k){
    coat<-onecoa(d);
    for(i in 1:m){
      pi<-coat[i,]; preC<-0;
      for(j in 1:d){
        if(j==d){value<-vary}
        else{
          cvary<-c()
          for(l in 1:No){
            y2<-c(); x_nj<-x.sample(pi[-c(1:j)])
            for(h in 1:Ni){
              x_j<-condition.sample(pi[1:j],x_nj,pi[-c(1:j)]); x_all<-c(x_j,x_nj);
              y2<-c(y2,model(x_all[sort(pi,index.return=TRUE)$ix]))
            }
            cvary<-c(cvary,var(y2))
          }
          value<-mean(cvary); 
        }
        delta<-value-preC; preC<-delta+preC;
        sh[pi[j]]<-sh[pi[j]]+delta; sh2[pi[j]]<-sh2[pi[j]]+delta^2;
      }
    }
  }
  sh<-sh/n; sh2<-sh2/n
  return(list(Shapley = sh, SEShapley = sqrt((sh2 - sh^2)/n)))
}


strict.est.shsa<-function(d,n,samplemethod,model,vary,x.sample,condition.sample, No, Ni){ 
  if(samplemethod=='SRS'){return(est.shsasrs(d,n,model,vary,x.sample,condition.sample, No, Ni))}
  else if(samplemethod=='LS'){return(est.shsals(d,n,model,vary,x.sample,condition.sample, No, Ni))}
  else {return(est.shsacoa(d,n,model,vary,x.sample,condition.sample, No, Ni))}
}

strictall.gtsa<-function(real,d,n,samplemethod,model,vary,x.sample,condition.sample, No, Ni,rt){
  nc<-d+3; nn<-length(n); nall<-nn*rt; nother<-nn;
  resultall<-matrix(0,nall,nc);  conave<-matrix(0,nother,nc);conmin<-matrix(0,nother,nc);
  conmax<-matrix(0,nother,nc); conmedian<-matrix(0,nother,nc); varsh<-matrix(0,nother,d+1)
  for(i in 1:nn){
    result<-matrix(0,rt,nc)
    for(j in 1:rt){
      t0<-Sys.time()
      sh<-strict.est.shsa(d,n[i],samplemethod,model,vary,x.sample,condition.sample, No, Ni)$Shapley; time<-Sys.time()-t0; mse<-sum((sh-real)^2)/d;
      result[j,]<-c(sh,mse,time,n[i])
    }
    varsh[i,]<-c(apply(result[,1:d],2,var),n[i]);
    conave[i,]<-apply(result,2,mean);
    conmedian[i,]<-apply(result,2,median);
    conmin[i,]<-apply(result,2,min);
    conmax[i,]<-apply(result,2,max);
    resultall[((i-1)*rt+1):(i*rt),]<-result;
  }
  colnames(resultall)<-c(1:d,'mse','time','size')
  rownames(conave)<-n; rownames(conmedian)<-n;rownames(conmin)<-n;rownames(conmax)<-n; rownames(varsh)<-n;
  colnames(conave)<-c(1:d,'mse','time','size');colnames(conmedian)<-c(1:d,'mse','time','size');colnames(conmin)<-c(1:d,'mse','time','size');colnames(conmax)<-c(1:d,'mse','time','size');colnames(varsh)<-c(1:d,'size');
  return(list(realvalue=real, result=data.frame(resultall,method=samplemethod), ave=data.frame(conave,method=samplemethod), median=data.frame(conmedian,method=samplemethod), min=data.frame(conmin,method=samplemethod), max=data.frame(conmax,method=samplemethod),var=varsh,method=data.frame(samplemethod)))
}


est.vary<-function(d,Nv,model,x.sample){
  y1<-c()
  for(i in 1:Nv){x<-x.sample(c(1:d)); y1<-c(y1,model(x))}
  Ey<-mean(y1); vary<-var(y1)
  return(vary)
}
#=================================Simulation====================================
mu<-c(1.2,1.5,4,1.8,3.6,1.5)
model.msmd5<-function(x){
  nu1<-x[1]; nu2<-0.4*x[1]+x[2]+x[3]+x[5];
  nu3<-0.3*x[1]+0.15*x[4];
  nu4<-0.6*x[1]+0.3*x[4];
  nu5<-x[4];
  nu6<-0.85*x[4]+0.3*x[1];
  y<-(nu1/(mu[1]-nu1)+nu2/(mu[2]-nu2)+nu3/(mu[3]-nu3)+nu4/(mu[4]-nu4)+nu5/(mu[5]-nu5)+nu6/(mu[6]-nu6))*24/(sum(x))
  return(y)
}

x.sample_msmd5<-function(sets){
  m<-length(sets); ss<-sort(sets); x<-c();
  mu_vector<-c(0,0); rho12<-0.5; rho34<--0.5; sigma12<-matrix(c(1,rho12,rho12,1),nrow=2); sigma34<-matrix(c(1,rho34,rho34,1),nrow=2);
  alpha<-10;beta<-19;
  I1<-0; I2<-0;
  if(m>1){
    for(i in 1:(m-1)){
      if((ss[i]==1)&(ss[i+1]==2)){
        t1<-which(sets==1); t2<-which(sets==2);
        z12<-mvrnorm(1,mu_vector,sigma12)
        x[c(t1,t2)]<-0.3*qbeta(pnorm(z12),alpha,beta)+0.5; I1<-1;
      }
      else if ((ss[i]==3)&(ss[i+1]==4)){
        t3<-which(sets==3); t4<-which(sets==4);
        z34<-mvrnorm(1,mu_vector,sigma34)
        x[c(t3,t4)]<-0.3*qbeta(pnorm(z34),alpha,beta)+0.5; I2<-1;
      }
    }
  }  
  for(i in 1:m){
    if((sets[i]==1)&(I1==0)) x[i]<-0.3*rbeta(1,alpha,beta)+0.5;
    if((sets[i]==2)&(I1==0)) x[i]<-0.3*rbeta(1,alpha,beta)+0.5;
    if((sets[i]==3)&(I2==0)) x[i]<-0.3*rbeta(1,alpha,beta)+0.5;
    if((sets[i]==4)&(I2==0)) x[i]<-0.3*rbeta(1,alpha,beta)+0.5;
    if(sets[i]==5) x[i]<-0.3*rbeta(1,alpha,beta)+0.5;
  }
  return(x)
}

condition.sample_msmd5<-function(sets,xn,xn_id){
  m<-length(sets); ss<-sort(sets); x<-c();
  mu_vector<-c(0,0); rho12<-0.5; rho34<--0.5; sigma12<-matrix(c(1,rho12,rho12,1),nrow=2); sigma34<-matrix(c(1,rho34,rho34,1),nrow=2);
  alpha<-10;beta<-19;x1<-xn[which(xn_id==1)];x2<-xn[which(xn_id==2)];x3<-xn[which(xn_id==3)];x4<-xn[which(xn_id==4)];
  I1<-0; I2<-0;
  if(m>1){
    for(i in 1:(m-1)){
      if((ss[i]==1)&(ss[i+1]==2)){
        t1<-which(sets==1); t2<-which(sets==2);
        z12<-mvrnorm(1,mu_vector,sigma12)
        x[c(t1,t2)]<-0.3*qbeta(pnorm(z12),alpha,beta)+0.5; I1<-1;
      }
      else if ((ss[i]==3)&(ss[i+1]==4)){
        t3<-which(sets==3); t4<-which(sets==4);
        z34<-mvrnorm(1,mu_vector,sigma34)
        x[c(t3,t4)]<-0.3*qbeta(pnorm(z34),alpha,beta)+0.5; I2<-1;
      }
    }
  }  
  for(i in 1:m){
    if((sets[i]==1)&(I1==0)){
      z2<-qnorm(pbeta(10*x2/3-5/3,alpha,beta));mu1<-rho12*z2;sigma1<-1-rho12^2;
      z1<-rnorm(1,mu1,sigma1); x[i]<-0.3*qbeta(pnorm(z1),alpha,beta)+0.5
    }
    if((sets[i]==2)&(I1==0)){
      z1<-qnorm(pbeta(10*x1/3-5/3,alpha,beta));mu2<-rho12*z1;sigma2<-1-rho12^2;
      z2<-rnorm(1,mu2,sigma2); x[i]<-0.3*qbeta(pnorm(z2),alpha,beta)+0.5
    }
    if((sets[i]==3)&(I2==0)){
      z4<-qnorm(pbeta(10*x4/3-5/3,alpha,beta));mu3<-rho34*z4;sigma3<-1-rho34^2;
      z3<-rnorm(1,mu3,sigma3); x[i]<-0.3*qbeta(pnorm(z3),alpha,beta)+0.5
    }
    if((sets[i]==4)&(I2==0)){
      z3<-qnorm(pbeta(10*x3/3-5/3,alpha,beta));mu4<-rho34*z3;sigma4<-1-rho34^2;
      z4<-rnorm(1,mu4,sigma4); x[i]<-0.3*qbeta(pnorm(z4),alpha,beta)+0.5
    }
    if(sets[i]==5) x[i]<-0.3*rbeta(1,alpha,beta)+0.5;
  }
  return(x)
}

Nv<-10000; No<-500; Ni<-10; d<-5;
vary_songmsmd5<-0
for(t in 1:10){
  vary_songmsmd5<-vary_songmsmd5+est.vary(d,Nv,model.msmd5,x.sample_msmd5)
}
vary_songmsmd5<-vary_songmsmd5/10; vary_songmsmd5

real_sh_msmd5<-est.shsareal(d,model.msmd5,vary_songmsmd5,x.sample_msmd5,condition.sample_msmd5,No,Ni)
real_sh_msmd5

No<-100; Ni<-2; d<-5;
n1<-c(c(1:(6*(d-1)))*d)
n1
#SRS
strictallresult_msmd5_srs<-strictall.gtsa(real_sh_msmd5$Shapley,d,n1,'SRS',model.msmd5,vary_songmsmd5,x.sample_msmd5,condition.sample_msmd5,No,Ni,50)
#rls
strictallresult_msmd5_rls<-strictall.gtsa(real_sh_msmd5$Shapley,d,n1,'LS',model.msmd5,vary_songmsmd5,x.sample_msmd5,condition.sample_msmd5,No,Ni,50)
#coa
n2<-c(c(1:6)*d*(d-1))
strictallresult_msmd5_coa<-strictall.gtsa(real_sh_msmd5$Shapley,d,n2,'COA',model.msmd5,vary_songmsmd5,x.sample_msmd5,condition.sample_msmd5,No,Ni,50)

allresult_msmd5<-rbind(strictallresult_msmd5_srs$result,strictallresult_msmd5_rls$result,strictallresult_msmd5_coa$result)
#====================================Figure S3(a)==============================================
allresult_msmd5_n20<-allresult_msmd5[allresult_msmd5$size==20,]
allest_msmd5_n20<-c(); nc<-ncol(allresult_msmd5_n20);
for(i in 1:d){
  t<-allresult_msmd5_n20[,c(i,(d+1):nc)]; names(t)[1]<-'estimate';
  allest_msmd5_n20<-rbind(allest_msmd5_n20,data.frame(t,variable=paste('x',i,sep='')))
}
tnames<-c()
for(i in 1:d){tnames<-c(tnames,paste('x',i,sep=''))}
ggplot(data=allest_msmd5_n20,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=tnames)+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))


#====================================Figure S3(b)==========================================
n2<-c(c(1:6)*d*(d-1))
allresult_msmd5_partial<-allresult_msmd5[allresult_msmd5$size%in%n2,]
allresult_msmd5_partial$size[allresult_msmd5_partial$size==20]<-'1m(m-1)'
allresult_msmd5_partial$size[allresult_msmd5_partial$size==40]<-'2m(m-1)'
allresult_msmd5_partial$size[allresult_msmd5_partial$size==60]<-'3m(m-1)'
allresult_msmd5_partial$size[allresult_msmd5_partial$size==80]<-'4m(m-1)'
allresult_msmd5_partial$size[allresult_msmd5_partial$size==100]<-'5m(m-1)'
allresult_msmd5_partial$size[allresult_msmd5_partial$size==120]<-'6m(m-1)'
allresult_msmd5_partial$size<-as.character(allresult_msmd5_partial$size)
allresult_msmd5_partial_clear<-allresult_msmd5_partial[-which.max(allresult_msmd5_partial$mse),]

ggplot(data=allresult_msmd5_partial_clear,aes(x=size, y=d*mse, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(labels=n2)+
  labs(x='Sample size',y='Square loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))

#====================================Table S2==========================================
#n=20
var.mean_msmd5_n20<-matrix(nrow=3,ncol=5)
allresult_msmd5_n20_SRS<-allresult_msmd5_n20[allresult_msmd5_n20$method=='SRS',]
allresult_msmd5_n20_LS<-allresult_msmd5_n20[allresult_msmd5_n20$method=='LS',]
allresult_msmd5_n20_COA<-allresult_msmd5_n20[allresult_msmd5_n20$method=='COA',]
var.mean_msmd5_n20[1,]<-apply(allresult_msmd5_n20_SRS[,1:5],2,var)
var.mean_msmd5_n20[2,]<-apply(allresult_msmd5_n20_LS[,1:5],2,var)
var.mean_msmd5_n20[3,]<-apply(allresult_msmd5_n20_COA[,1:5],2,var)

#n=40
var.mean_msmd5_n40<-matrix(nrow=3,ncol=5)
allresult_msmd5_n40_SRS<-allresult_msmd5_n40[allresult_msmd5_n40$method=='SRS',]
allresult_msmd5_n40_LS<-allresult_msmd5_n40[allresult_msmd5_n40$method=='LS',]
allresult_msmd5_n40_COA<-allresult_msmd5_n40[allresult_msmd5_n40$method=='COA',]
var.mean_msmd5_n40[1,]<-apply(allresult_msmd5_n40_SRS[,1:5],2,var)
var.mean_msmd5_n40[2,]<-apply(allresult_msmd5_n40_LS[,1:5],2,var)
var.mean_msmd5_n40[3,]<-apply(allresult_msmd5_n40_COA[,1:5],2,var)

#n=60
var.mean_msmd5_n60<-matrix(nrow=3,ncol=5)
allresult_msmd5_n60_SRS<-allresult_msmd5_n60[allresult_msmd5_n60$method=='SRS',]
allresult_msmd5_n60_LS<-allresult_msmd5_n60[allresult_msmd5_n60$method=='LS',]
allresult_msmd5_n60_COA<-allresult_msmd5_n60[allresult_msmd5_n60$method=='COA',]
var.mean_msmd5_n60[1,]<-apply(allresult_msmd5_n60_SRS[,1:5],2,var)
var.mean_msmd5_n60[2,]<-apply(allresult_msmd5_n60_LS[,1:5],2,var)
var.mean_msmd5_n60[3,]<-apply(allresult_msmd5_n60_COA[,1:5],2,var)

#n=80
var.mean_msmd5_n80<-matrix(nrow=3,ncol=5)
allresult_msmd5_n80_SRS<-allresult_msmd5_n80[allresult_msmd5_n80$method=='SRS',]
allresult_msmd5_n80_LS<-allresult_msmd5_n80[allresult_msmd5_n80$method=='LS',]
allresult_msmd5_n80_COA<-allresult_msmd5_n80[allresult_msmd5_n80$method=='COA',]
var.mean_msmd5_n80[1,]<-apply(allresult_msmd5_n80_SRS[,1:5],2,var)
var.mean_msmd5_n80[2,]<-apply(allresult_msmd5_n80_LS[,1:5],2,var)
var.mean_msmd5_n80[3,]<-apply(allresult_msmd5_n80_COA[,1:5],2,var)

#n=100
var.mean_msmd5_n100<-matrix(nrow=3,ncol=5)
allresult_msmd5_n100_SRS<-allresult_msmd5_n100[allresult_msmd5_n100$method=='SRS',]
allresult_msmd5_n100_LS<-allresult_msmd5_n100[allresult_msmd5_n100$method=='LS',]
allresult_msmd5_n100_COA<-allresult_msmd5_n100[allresult_msmd5_n100$method=='COA',]
var.mean_msmd5_n100[1,]<-apply(allresult_msmd5_n100_SRS[,1:5],2,var)
var.mean_msmd5_n100[2,]<-apply(allresult_msmd5_n100_LS[,1:5],2,var)
var.mean_msmd5_n100[3,]<-apply(allresult_msmd5_n100_COA[,1:5],2,var)

#n=120
var.mean_msmd5_n120<-matrix(nrow=3,ncol=5)
allresult_msmd5_n120_SRS<-allresult_msmd5_n120[allresult_msmd5_n120$method=='SRS',]
allresult_msmd5_n120_LS<-allresult_msmd5_n120[allresult_msmd5_n120$method=='LS',]
allresult_msmd5_n120_COA<-allresult_msmd5_n120[allresult_msmd5_n120$method=='COA',]
var.mean_msmd5_n120[1,]<-apply(allresult_msmd5_n120_SRS[,1:5],2,var)
var.mean_msmd5_n120[2,]<-apply(allresult_msmd5_n120_LS[,1:5],2,var)
var.mean_msmd5_n120[3,]<-apply(allresult_msmd5_n120_COA[,1:5],2,var)

var.mean_msmd5<-data.frame(rbind(var.mean_msmd5_n20,var.mean_msmd5_n40,var.mean_msmd5_n60,var.mean_msmd5_n80,
                                     var.mean_msmd5_n100,var.mean_msmd5_n120),
                               method=rep(c('SRS','LS','COA'),6),size=rep(c(20,40,60,80,100,120),each=3))
View(var.mean_msmd5)

var_msmd5<-as.matrix(var.mean_msmd5[,1:5])
ratio_var_msmd5<-cbind(var_msmd5,apply(var_msmd5,1,sum))
ratio_var_msmd5[1:3,]<-t(t(ratio_var_msmd5[1:3,])/ratio_var_msmd5[3,])
ratio_var_msmd5[4:6,]<-t(t(ratio_var_msmd5[4:6,])/ratio_var_msmd5[6,])
ratio_var_msmd5[7:9,]<-t(t(ratio_var_msmd5[7:9,])/ratio_var_msmd5[9,])
ratio_var_msmd5[10:12,]<-t(t(ratio_var_msmd5[10:12,])/ratio_var_msmd5[12,])
ratio_var_msmd5[13:15,]<-t(t(ratio_var_msmd5[13:15,])/ratio_var_msmd5[15,])
ratio_var_msmd5[16:18,]<-t(t(ratio_var_msmd5[16:18,])/ratio_var_msmd5[18,])
ratio_var_msmd5<-data.frame(method=rep(c('SRS','LS','COA'),6),size=rep(c(20,40,60,80,100,120),each=3),
                                ratio_var_msmd5)

View(ratio_var_msmd5)