##################################################################
##########Code for Application(1): 911 terrorist with d=69########
##################################################################
library(foreach)
library(doParallel)
library(ggplot2)
library(stringr)
#=========================Functions===========================

val.social<-function(sets,adjacent){
  if(length(sets)==1) val<-0
  else{
    subadjacent<-adjacent[sets,sets]
    nsets<-length(sets)
    A<-diag(1,nsets); B<-matrix(0,nsets,nsets)
    for(l in 1:(nsets-1)){
      A<-A%*%subadjacent
      B<-B+A
    }
    val<-ifelse(sum(B==0)>nsets,0,1)
  }
  return(val)
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

est.soshsrs<-function(d,n,val,adjacent){
  sh<-rep(0,d)
  for(l in 1:n){
    perml<-sample(d); preC<-0;
    for(i in 1:d){
      delta<-val(perml[1:i],adjacent)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

est.soshls<-function(d,n,val,adjacent){
  sh<-rep(0,d); k<-n/d;
  for(t in 1:k){
    lst<-onels(d);
    for(l in 1:d){
      perml<-lst[l,]; preC<-0;
      for(i in 1:d){
        delta<-val(perml[1:i],adjacent)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

est.soshstrrs<-function(d,n,val,adjacent){
  sh<-rep(0,d); groupsize<-n/d
  srsperm<-matrix(0,nrow=n,ncol=d)
  for(i in 1:n){
    srsperm[i,]<-sample(d)
  }
  for(j in 1:d){
    jperm<-structed.perm(srsperm,j,d)
    for(i in (groupsize+1):n){
      loc<-ceiling(i/groupsize)
      sh[j]<-sh[j]+val(jperm[i,1:loc],adjacent)-val(jperm[i,1:(loc-1)],adjacent)
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:d);
  return(sh)
}

strict.est.sosh<-function(d,n,samplemethod,val,adjacent){ 
  if(samplemethod=='SRS'){return(est.soshsrs(d,n,val,adjacent))}
  else if(samplemethod=='LS'){return(est.soshls(d,n,val,adjacent))}
  else if (samplemethod=='StrRS'){return(est.soshstrrs(d,n,val,adjacent))}
}

strict.so.bx<-function(d,n,samplemethod,val,rt,adjacent,cores){
  nc<-d+3; nn<-length(n); nall<-nn*rt; nother<-nn;
  resultall<-matrix(0,nall,nc);  
  for(i in 1:nn){
    cores<-min(detectCores(logical=F),cores)
    cl <- makeCluster(cores)
    registerDoParallel(cl,cores=sepnum)
    result<-foreach(j=1:rt, .combine= "rbind",.export =c("strict.est.sosh","est.soshsrs","est.soshstrrs","est.soshls","onels","structed.perm")) %dopar%
      {
        t0<-Sys.time()
        sh<-strict.est.sosh(d,n[i],samplemethod,val,adjacent); 
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


est.soshcoa.71null69design<-function(d,n,val,adjacent){
  sh<-rep(0,(d-2)); k<-n/d/(d-1);m<-d*(d-1);
  for(t in 1:k){
    coat<-onecoa(d);
    for(l in 1:m){
      perml<-setdiff(coat[l,],c(70,71)); preC<-0;
      for(i in 1:(d-2)){
        delta<-val(perml[1:i],adjacent)-preC; sh[perml[i]]<-sh[perml[i]]+delta; preC<-delta+preC;
      }
    }
  }
  sh<-sh/n;  sh<-t(as.matrix(sh)); colnames(sh)<-c(1:(d-2));
  return(sh)
}


strict.so.coa.bx<-function(d,n,samplemethod,val,rt,adjacent,cores){
  nc<-d+1; nn<-length(n); nall<-nn*rt; nother<-nn;
  resultall<-matrix(0,nall,nc);  
  for(i in 1:nn){
    cores<-min(detectCores(logical=F),cores)
    cl <- makeCluster(cores)
    registerDoParallel(cl,cores=sepnum)
    result<-foreach(j=1:rt, .combine= "rbind",.export =c("est.soshcoa.71null69design","onecoa")) %dopar%
      {
        t0<-Sys.time()
        sh<-est.soshcoa.71null69design(d,n[i],samplemethod,val,adjacent); 
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
  colnames(resultall)<-c('size','rt','time',1:(d-2))
  return(data.frame(resultall,method=samplemethod))
}


top.diffset<-function(t,top,etop){
  nr<-nrow(etop); if(is.null(nr)) nr<-1;
  d<-length(top);nt<-length(t)
  diffmatrix<-matrix(0,nrow=nr,ncol=nt)
  if(nr==1){
    for(j in 1:nt){
      diffmatrix[1,j]<-t[j]-sum(etop[1:t[j]]%in%top[1:t[j]])
    }
  }
  else{
    for(i in 1:nr){
      for(j in 1:nt){
        diffmatrix[i,j]<-t[j]-sum(etop[i,1:t[j]]%in%top[1:t[j]])
      }
    }
  }
  colnames(diffmatrix)<-c(paste('t=',t,sep=""))
  return(diffmatrix)
}


wrong.exactperc<-function(diffresult,case,nrep){
  nr<-nrow(diffresult);nsize<-nr/nrep; nc<-length(case);
  perc<-matrix(0,nrow=nsize,ncol=nc)
  for(i in 1:nsize){
    for(j in 1:nc) perc[i,j]<-sum(diffresult[((i-1)*nrep+1):(i*nrep),j])/case[j]/nrep
  }
  colnames(perc)<-colnames(diffresult[,1:nc])
  return(perc)
}

#=========================Simulation===========================
# ¡®source_target.xlsx¡¯ and ¡¯id_label.xlsx¡¯ are provided by Valdis Krebs and Herbert Hamers, 
# which including the connections and person names of 69 nodes in the social network of the 9/11 terrorist attacks respectively.
# Readers can send us an email to request for data. 

#vertex_name
source_vertex_name<-read_excel("id_label.xlsx")
source_vertex_name<-as.matrix(source_vertex_name)

#edge-vertex
source_from_to <- read_excel("source_target.xlsx")
source_from_to<-as.matrix(source_from_to)
source_vertex_edge<-matrix(NA,nrow=69,ncol=25)
nedge<-rep(0,69)
for(i in 1:69){
  for(j in 1:nrow(source_from_to)){
    if(source_from_to[j,1]==i) {
      nedge[i]<-nedge[i]+1
      source_vertex_edge[i,nedge[i]]<-source_from_to[j,2]
    }
    if(source_from_to[j,2]==i) {
      nedge[i]<-nedge[i]+1
      source_vertex_edge[i,nedge[i]]<-source_from_to[j,1]
    }
  }
}
source_vertex_edge<-cbind(c(1:69),source_vertex_edge)

#adjacent_matrix
source_nv<-nrow(source_vertex_edge)
source_adjacent_matrix<-matrix(0,nrow=source_nv,ncol=source_nv)
for(i in 1:source_nv){
  source_adjacent_matrix[i,source_vertex_edge[i,!is.na(source_vertex_edge[i,])]]<-1
}
diag(source_adjacent_matrix)<-0
colnames(source_adjacent_matrix)<-source_vertex_name[,2]
rownames(source_adjacent_matrix)<-source_vertex_name[,2]
write.xlsx(source_adjacent_matrix,"source_adjacent_matrix.xlsx")


dso<-69; dso_71<-71;
nso_4methodnear71_coa_modified<-71*70*c(1:5)
nso_4methodnear71_69_modified<-seq(floor(71*70/69)*69,floor(71*70*5/69)*69,72*69)
n_nso_4methodnear71_69_modified<-length(nso_4methodnear71_69_modified)
n_nso_4methodnear71_coa_modified<-length(nso_4methodnear71_coa_modified)

#srs
n_mz<-40;n_zs<-300/n_mz;
for(i in 1:n_nso_4methodnear71_69_modified){
  assign(paste("result_so_srs_", i, sep=""), c())
}
time_srs<-matrix(0,nrow=n_zs,ncol=n_nso_4methodnear71_69_modified)
for(ii in 1:n_zs){
  for(jj in 1:n_nso_4methodnear71_69_modified){
    print(paste(paste("round:",ii,sep=""),paste("; m:",nso_4methodnear71_69_modified[jj],sep=""),sep=""))
    time0_srs<-Sys.time()
    t_srs<-strict.so.bx(69,nso_4methodnear71_69_modified[jj],'SRS',val.social,n_mz,source_adjacent_matrix,40)
    t_srs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_so_srs_", jj, sep=""), rbind(get(paste("result_so_srs_", jj, sep="")),t_srs))
    time_srs[ii,jj]<-as.numeric(Sys.time()-time0_srs,units="secs")
    tsave_srs<-get(paste("result_so_srs_", jj, sep=""))
    save(tsave_srs,file=eval(paste(paste("result_so_srs_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_srs[ii,jj],sep=""))
  }
}
result_so_srs<-c()
for(i in 1:n_nso_4methodnear71_69_modified){
  result_so_srs<-rbind(result_so_srs,get(paste("result_so_srs_", i, sep="")))
}


#strrs
for(i in 1:n_nso_4methodnear71_69_modified){
  assign(paste("result_so_strrs_", i, sep=""), c())
}
time_strrs<-matrix(0,nrow=n_zs,ncol=n_nso_4methodnear71_69_modified)
for(ii in 1:n_zs){
  for(jj in 1:n_nso_4methodnear71_69_modified){
    print(paste(paste("round:",ii,sep=""),paste("; m:",nso_4methodnear71_69_modified[jj],sep=""),sep=""))
    time0_strrs<-Sys.time()
    t_strrs<-strict.so.bx(69,nso_4methodnear71_69_modified[jj],'StrRS',val.social,n_mz,source_adjacent_matrix,40)
    t_strrs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_so_strrs_", jj, sep=""), rbind(get(paste("result_so_strrs_", jj, sep="")),t_strrs))
    time_strrs[ii,jj]<-as.numeric(Sys.time()-time0_strrs,units="secs")
    tsave_strrs<-get(paste("result_so_strrs_", jj, sep=""))
    save(tsave_strrs,file=eval(paste(paste("result_so_strrs_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_strrs[ii,jj],sep=""))
  }
}
result_so_strrs<-c()
for(i in 1:n_nso_4methodnear71_69_modified){
  result_so_strrs<-rbind(result_so_strrs,get(paste("result_so_strrs_", i, sep="")))
}

#ls
for(i in 1:n_nso_4methodnear71_69_modified){
  assign(paste("result_so_ls_", i, sep=""), c())
}
time_ls<-matrix(0,nrow=n_zs,ncol=n_nso_4methodnear71_69_modified)
for(ii in 1:n_zs){
  for(jj in 1:n_nso_4methodnear71_69_modified){
    print(paste(paste("round:",ii,sep=""),paste("; m:",nso_4methodnear71_69_modified[jj],sep=""),sep=""))
    time0_ls<-Sys.time()
    t_ls<-strict.so.bx(69,nso_4methodnear71_69_modified[jj],'LS',val.social,n_mz,source_adjacent_matrix,40)
    t_ls[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_so_ls_", jj, sep=""), rbind(get(paste("result_so_ls_", jj, sep="")),t_ls))
    time_ls[ii,jj]<-as.numeric(Sys.time()-time0_ls,units="secs")
    tsave_ls<-get(paste("result_so_ls_", jj, sep=""))
    save(tsave_ls,file=eval(paste(paste("result_so_ls_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_ls[ii,jj],sep=""))
  }
}
result_so_ls<-c()
for(i in 1:n_nso_4methodnear71_69_modified){
  result_so_ls<-rbind(result_so_ls,get(paste("result_so_ls_", i, sep="")))
}

#coa
for(i in 1:n_nso_4methodnear71_coa_modified){
  assign(paste("result_so_coa_", i, sep=""), c())
}
time_coa<-matrix(0,nrow=n_zs,ncol=n_nso_4methodnear71_coa_modified)
for(ii in 1:n_zs){
  for(jj in 1:n_nso_4methodnear71_coa_modified){
    print(paste(paste("round:",ii,sep=""),paste("; m:",nso_4methodnear71_coa_modified[jj],sep=""),sep=""))
    time0_coa<-Sys.time()
    t_coa<-strict.so.coa.bx(dso_71,nso_4methodnear71_coa_modified[jj],'COA',val.social,n_mz,source_adjacent_matrix,40)
    t_coa[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_so_coa_", jj, sep=""), rbind(get(paste("result_so_coa_", jj, sep="")),t_coa))
    time_coa[ii,jj]<-as.numeric(Sys.time()-time0_coa,units="secs")
    tsave_coa<-get(paste("result_so_coa_", jj, sep=""))
    save(tsave_coa,file=eval(paste(paste("result_so_coa_", jj, sep=""),".RData",sep="")))
    print(paste("time(sec):",time_coa[ii,jj],sep=""))
  }
}
result_so_coa<-c()
for(i in 1:n_nso_4methodnear71_coa_modified){
  result_so_coa<-rbind(result_so_coa,get(paste("result_so_coa_", i, sep="")))
}


#pseudo
pseudo_real_sh_911<-rep(0,69)
for(i in 1:n_nso_4methodnear71_69_modified){
  for(j in 1:300){
    pseudo_real_sh_911<-pseudo_real_sh_911+result_so_srs[(i-1)*300+j,3+1:69]*nso_4methodnear71_69_modified[i]
  }
}
pseudo_real_sh_911<-pseudo_real_sh_911/sum(nso_4methodnear71_69_modified)/300
pseudo_real_sh_911

pseudo_top_source<-sort(as.numeric(pseudo_real_sh_911),decreasing = TRUE,index.return=TRUE)$ix
pseudo_top_source
pseudo_top20_shid_so<-data.frame(sort=c(1:20),id=source_vertex_name[pseudo_top_source[1:20],1],name=source_vertex_name[pseudo_top_source[1:20],2],sh=as.numeric(pseudo_real_sh_911[pseudo_top_source[1:20]]))
pseudo_top20_shid_so
pseudo_top25_shid_so<-data.frame(sort=c(1:25),id=source_vertex_name[pseudo_top_source[1:25],1],name=source_vertex_name[pseudo_top_source[1:25],2],sh=as.numeric(pseudo_real_sh_911[pseudo_top_source[1:25]]))
pseudo_top25_shid_so


#sort
#srs
allorder_so_srs<-matrix(0,nrow=nrow(result_so_srs),ncol=dso)
for(i in 1:nrow(allorder_so_srs)){
  allorder_so_srs[i,]<-sort(as.numeric(result_so_srs[i,3+1:dso]),decreasing = TRUE,index.return=TRUE)$ix
}
allorder_so_srs

#strrs
allorder_so_strrs<-matrix(0,nrow=nrow(result_so_strrs),ncol=dso)
for(i in 1:nrow(allorder_so_strrs)){
  allorder_so_strrs[i,]<-sort(as.numeric(result_so_strrs[i,3+1:dso]),decreasing = TRUE,index.return=TRUE)$ix
}
allorder_so_strrs

#ls
allorder_so_ls<-matrix(0,nrow=nrow(result_so_ls),ncol=dso)
for(i in 1:nrow(allorder_so_ls)){
  allorder_so_ls[i,]<-sort(as.numeric(result_so_ls[i,3+1:dso]),decreasing = TRUE,index.return=TRUE)$ix
}
allorder_so_ls

#coa
allorder_so_coa<-matrix(0,nrow=nrow(result_so_coa),ncol=dso)
for(i in 1:nrow(allorder_so_coa)){
  allorder_so_coa[i,]<-sort(as.numeric(result_so_coa[i,3+1:dso]),decreasing = TRUE,index.return=TRUE)$ix
}
allorder_so_coa

#=====================================Figure 3=====================================
social_cols<- c("LS" = "#00CC66", "SRS" = "#FF6666", "StrRS" = "#AB82FF", "COA" = "#6699FF")
social_shape<-c("LS" = 15, "SRS" =16, "StrRS" = 17, "COA" =18)

top_so_srs<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_srs)
top_so_strrs<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_strrs)
top_so_ls<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_ls)
top_so_coa<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_coa)
social_so_set_srs<-data.frame(top5=top_so_srs[,1],top9=top_so_srs[,2],top10=top_so_srs[,3],top12=top_so_srs[,4],top15=top_so_srs[,5],top20=top_so_srs[,6],top24=top_so_srs[,7],size=rep(nso_4methodnear71_69_modified,each=300),method='SRS')
social_so_set_strrs<-data.frame(top5=top_so_strrs[,1],top9=top_so_strrs[,2],top10=top_so_strrs[,3],top12=top_so_strrs[,4],top15=top_so_strrs[,5],top20=top_so_strrs[,6],top24=top_so_strrs[,7],size=rep(nso_4methodnear71_69_modified,each=300),method='StrRS')
social_so_set_ls<-data.frame(top5=top_so_ls[,1],top9=top_so_ls[,2],top10=top_so_ls[,3],top12=top_so_ls[,4],top15=top_so_ls[,5],top20=top_so_ls[,6],top24=top_so_ls[,7],size=rep(nso_4methodnear71_69_modified,each=300),method='LS')
social_so_set_coa<-data.frame(top5=top_so_coa[,1],top9=top_so_coa[,2],top10=top_so_coa[,3],top12=top_so_coa[,4],top15=top_so_coa[,5],top20=top_so_coa[,6],top24=top_so_coa[,7],size=rep(nso_4methodnear71_coa_modified,each=300),method='COA')


social_so_set_1<-rbind(social_so_set_srs,social_so_set_strrs,social_so_set_ls)
social_so_set_1$size<-as.character(social_so_set_1$size)
social_so_set_coa$size<-as.character(social_so_set_coa$size)


social_set<-rbind(social_so_set_1,social_so_set_coa)
social_set$size<-as.character(social_set$size)

social_so_percset_exact_1<-data.frame(wrong.exactperc(social_so_set_1,c(5,9,10,12,15,20,24),300),size=rep(nso_4methodnear71_69_modified,3),method=rep(c('SRS','StrRS','LS'),each=length(nso_4methodnear71_69_modified)))
social_so_percset_exact_coa<-data.frame(wrong.exactperc(social_so_set_coa,c(5,9,10,12,15,20,24),300),size=nso_4methodnear71_coa_modified,method='COA')


social_so_percset_exact<-rbind(social_so_percset_exact_coa,social_so_percset_exact_1)

social_so_percset_exact[social_so_percset_exact$size==4968,]$size<-'4968(4970)'
social_so_percset_exact[social_so_percset_exact$size==4970,]$size<-'4968(4970)'
social_so_percset_exact[social_so_percset_exact$size==c(9936),]$size<-'9936(9940)'
social_so_percset_exact[social_so_percset_exact$size==c(9940),]$size<-'9936(9940)'
social_so_percset_exact[social_so_percset_exact$size==c(14904),]$size<-'14904(14910)'
social_so_percset_exact[social_so_percset_exact$size==c(14910),]$size<-'14904(14910)'
social_so_percset_exact[social_so_percset_exact$size==c(19872),]$size<-'19872(19880)'
social_so_percset_exact[social_so_percset_exact$size==c(19880),]$size<-'19872(19880)'
social_so_percset_exact[social_so_percset_exact$size==c(24840),]$size<-'24840(24850)'
social_so_percset_exact[social_so_percset_exact$size==c(24850),]$size<-'24840(24850)'

social_so_percset_exact$size<-as.factor(social_so_percset_exact$size)
levels(social_so_percset_exact$size)


#top 12
ggplot(data=social_so_percset_exact,aes(x=factor(size,level=size_so), y=top12*100,color=factor(method,levels=c('SRS','StrRS','LS','COA')),group=factor(method,levels=c('SRS','StrRS','LS','COA')),shape=factor(method,levels=c('SRS','StrRS','LS','COA')),linetype=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=size_so)+
  scale_color_manual(values=social_cols)+scale_shape_manual(values=social_shape)+
  labs(x='Sample size',y='Error rate(%)')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.88,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 11))

#top 20
ggplot(data=social_so_percset_exact,aes(x=factor(size,level=size_so), y=top20*100,color=factor(method,levels=c('SRS','StrRS','LS','COA')),group=factor(method,levels=c('SRS','StrRS','LS','COA')),shape=factor(method,levels=c('SRS','StrRS','LS','COA')),linetype=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=size_so)+
  scale_color_manual(values=social_cols)+scale_shape_manual(values=social_shape)+
  labs(x='Sample size',y='Error rate(%)')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.88,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 11))


#=======================================Figure 4(a)===================================
allresult_so<-rbind(result_so_srs,result_so_strrs,result_so_ls,result_so_coa)
allresult_so_top<-allresult_so[,c(3+as.numeric(top_source),1,2,3+69+1)]
allresult_so_n19880<-allresult_so_top[(allresult_so_top$size==19880)|(allresult_so_top$size==19872),]
allresult_so_n19880$method<-as.factor(allresult_so_n19880$method)

allresult_so_n19880_top5<-c(); nc<-ncol(allresult_so_n19880);
for(i in 1:5){
  t<-allresult_so_n19880[,c(i,nc-2,nc)]; names(t)[1]<-'estimate';
  allresult_so_n19880_top5<-rbind(allresult_so_n19880_top5,data.frame(t,variable=source_vertex_name[as.numeric(pseudo_top_source)[i],2]))
}
pseudo_real_top5_so<-as.numeric(pseudo_real_sh_911[pseudo_top_source[1:5]])
allresult_so_n19880_top5<-data.frame(allresult_so_n19880_top5,real=rep(pseudo_real_top5_so,each=300*4))
names_top5_so<-source_vertex_name[as.numeric(pseudo_top_source)[1:5],2]
allresult_so_n19880_top5$variable <- factor(allresult_so_n19880_top5$variable,levels=names_top5_so)
ggplot(data=allresult_so_n19880_top5,aes(x=variable, y=estimate, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=names_top5_so,labels = function(x) str_wrap(x, width = 10)) +
  scale_color_manual(values=social_cols)+scale_fill_manual(values=social_cols)+
  labs(x='Member',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  geom_errorbar(aes(y=real, ymax=real, ymin=real),linetype='dashed',color='black',size=0.7)+
  guides(colour=guide_legend(nrow=2))+
  theme(legend.position=c(0.81,0.83),text=element_text(size = 17),axis.text.x = element_text(size = 11))


#======================================================Figure 4(b)================================================
allmse_so<-c()
for(i in 1:nrow(allresult_so)){
  allmse_so<-c(allmse_so,sum((allresult_so[i,3+c(1:69)]-pseudo_real_sh_911)^2))
}
allmse_so<-data.frame(mse=allmse_so,size=allresult_so$size,rt=allresult_so$rt,method=allresult_so$method)

allmse_so[allmse_so$size==4968,]$size<-'4968(4970)'
allmse_so[allmse_so$size==4970,]$size<-'4968(4970)'
allmse_so[allmse_so$size==c(9936),]$size<-'9936(9940)'
allmse_so[allmse_so$size==c(9940),]$size<-'9936(9940)'
allmse_so[allmse_so$size==c(14904),]$size<-'14904(14910)'
allmse_so[allmse_so$size==c(14910),]$size<-'14904(14910)'
allmse_so[allmse_so$size==c(19872),]$size<-'19872(19880)'
allmse_so[allmse_so$size==c(19880),]$size<-'19872(19880)'
allmse_so[allmse_so$size==c(24840),]$size<-'24840(24850)'
allmse_so[allmse_so$size==c(24850),]$size<-'24840(24850)'

allmse_so$method<-as.factor(allmse_so$method)
allmse_so$size<-as.factor(allmse_so$size)

size_so<-as.factor(c('4968(4970)','9936(9940)','14904(14910)','19872(19880)','24840(24850)'))
ggplot(data=allmse_so,aes(x=factor(size,levels=size_so), y=mse, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_so)+
  scale_color_manual(values=social_cols)+scale_fill_manual(values=social_cols)+
  labs(x='Sample size',y='Squared loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 11))

