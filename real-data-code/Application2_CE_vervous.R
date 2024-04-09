##########################################################################
##########Code for Application(2): C.E nervous system with d=279##########
##########################################################################
library(openxlsx)
library(ggplot2)
library(MASS)
library(foreach)
library(doParallel)
library(stringr)
library(igraph)

####================================Functions========================================
##---------------------0.DFS-----------------------
#judge whether the network is connected using DFS method
DFS.f<-function(visited,v,adjacent){
  visited[v]<-1; num<-length(visited);
  for(i in 1:num){
    if((adjacent[v,i]==1)&&(visited[i]==0)){
      visited<-DFS.f(visited,i,adjacent)
    }
    if(i==num){
      return(visited)
    }
  }
}

##---------------------1.value function-----------------------
val.net<-function(sets,adjacent){
  DFS<-function(visited,v,adjacent){
    visited[v]<-1; num<-length(visited);
    for(i in 1:num){
      if((adjacent[v,i]==1)&&(visited[i]==0)){
        visited<-DFS(visited,i,adjacent)
      }
      if(i==num){
        return(visited)
      }
    }  
  }
  nsets<-length(sets) 
  if(nsets==1) value<-0
  else{
    subadjacent<-adjacent[sets,sets];
    visited<-rep(0,nsets);
    value<-ifelse(sum(DFS(visited,1,subadjacent))==nsets,1,0)
  }
  return(value)
}

##---------------------2.estimate of shapley---------------------
#---------------------------2.1 SRS---------------------------
est.shsrs.net<-function(d,n,val,adjacent){
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

est.shls.net<-function(d,n,val,adjacent){
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


#---------------------------2.3 Summary----------------------------
strict.est.sh.net<-function(d,n,samplemethod,val,adjacent){ 
  if(samplemethod=='SRS'){return(est.shsrs.net(d,n,val,adjacent))}
  else if(samplemethod=='LS'){return(est.shls.net(d,n,val,adjacent))}
  else{stop("Error: please change the method.")}
}

#---------------------------2.4 Parallel computing---------------------
strict.net.bx<-function(d,n,samplemethod,val,rt,adjacent,cores){
  nc<-d+3; nn<-length(n); nall<-nn*rt; nother<-nn;
  resultall<-matrix(0,nall,nc);  
  for(i in 1:nn){
    cores<-min(detectCores(logical=F),cores)
    cl <- makeCluster(cores)
    registerDoParallel(cl,cores=sepnum)
    result<-foreach(j=1:rt, .combine= "rbind",.export =c("strict.est.sh.net","est.shsrs.net","est.shls.net","onels")) %dopar%
      {
        t0<-Sys.time()
        sh<-strict.est.sh.net(d,n[i],samplemethod,val,adjacent); 
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





#--------------------------2.5 top set and error rate --------------------------
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

####================================Data========================================
# The data set 'Connectome.csv' are downloaded from https://github.com/3BIM20162017/CElegansTP/tree/master/data.
net_connection <- read.csv("Connectome.csv", header=TRUE)

#vertex
net_allvertex<-c(levels(net_connection$Neuron),levels(net_connection$Target))
net_vertex<-net_allvertex[!duplicated(net_allvertex)]
net_num_vertex<-length(net_vertex)

#adjacent matrix
net_adjacent<-matrix(0,nrow=net_num_vertex,ncol=net_num_vertex)
for(i in 1:net_num_vertex){
  for(j in 1:nrow(net_connection)){
    if(net_connection[j,2]==net_vertex[i]) {
      tid<-which(net_vertex==net_connection[j,3])
      net_adjacent[i,tid]<-1;
      net_adjacent[tid,i]<-1;
    }
  }
}

num_edge<-sum(net_adjacent==1)
View(net_adjacent)

#choose a connected subgraph
gt<- graph.adjacency(net_adjacent,mode="undirected")
plot(gt,vertex.label=net_vertex)

part1<-which(DFS.f(rep(0,299),1,net_adjacent)==1)
n_part1<-length(part1)
net_adjacent_part1<-net_adjacent[part1,part1]
net_vertex_part1<-net_vertex[part1]



####================================Simulations================================
d_net<-n_part1
n_net<-d_net*c(1,c(1:4)*5)
n_nnet<-length(n_net)

##----------------------1.srs---------------
n_mz<-30;n_zs<-90/n_mz;


for(i in 1:n_nnet){
  assign(paste("result_net_srs_", i, sep=""), c())
}
time_net_srs<-matrix(0,nrow=n_zs,ncol=n_nnet)
for(ii in 1:n_zs){
  for(jj in 1:n_nnet){
    print(paste(paste("round:",ii,sep=""),paste("; m:",n_net[jj],sep=""),sep=""))
    time0_srs<-Sys.time()
    t_srs<-strict.net.bx(d_net,n_net[jj],'SRS',val.net,n_mz,net_adjacent_part1,30)
    t_srs[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_net_srs_", jj, sep=""), rbind(get(paste("result_net_srs_", jj, sep="")),t_srs))
    time_net_srs[ii,jj]<-as.numeric(Sys.time()-time0_srs,units="mins")
    tsave_net_srs<-get(paste("result_net_srs_", jj, sep=""))
    save(tsave_net_srs,file=eval(paste(paste("result_net_srs_", jj, sep=""),".RData",sep="")))
    print(paste("time(mins):",time_net_srs[ii,jj],sep=""))
    print(paste("result:",t_srs[1,],sep=""))
  }
}
result_net_srs<-c()
for(i in 1:n_nnet){
  result_net_srs<-rbind(result_net_srs,get(paste("result_net_srs_", i, sep="")))
}


##----------------------2.ls---------------
n_mz<-30;n_zs<-90/n_mz;


for(i in 1:n_nnet){
  assign(paste("result_net_ls_", i, sep=""), c())
}
time_net_ls<-matrix(0,nrow=n_zs,ncol=n_nnet)
for(ii in 1:n_zs){
  for(jj in 1:n_nnet){
    print(paste(paste("round:",ii,sep=""),paste("; m:",n_net[jj],sep=""),sep=""))
    time0_ls<-Sys.time()
    t_ls<-strict.net.bx(d_net,n_net[jj],'LS',val.net,n_mz,net_adjacent_part1,30)
    t_ls[,2]<-c(((ii-1)*n_mz+1):(ii*n_mz))
    assign(paste("result_net_ls_", jj, sep=""), rbind(get(paste("result_net_ls_", jj, sep="")),t_ls))
    time_net_ls[ii,jj]<-as.numeric(Sys.time()-time0_ls,units="mins")
    tsave_net_ls<-get(paste("result_net_ls_", jj, sep=""))
    save(tsave_net_ls,file=eval(paste(paste("result_net_ls_", jj, sep=""),".RData",sep="")))
    print(paste("time(mins):",time_net_ls[ii,jj],sep=""))
    print(paste("result:",t_ls[1,],sep=""))
  }
}
result_net_ls<-c()
for(i in 1:n_nnet){
  result_net_ls<-rbind(result_net_ls,get(paste("result_net_ls_", i, sep="")))
}

####================================Analysis================================
##----------------------1.analysis------------------------
d_net<-n_part1
allresult_net<-rbind(result_net_srs,result_net_ls)

size_net<-as.factor(n_net)
cols_net<- c("LS" = "#00CC66", "SRS" = "#FF6666")
shape_net<-c("LS" = 15, "SRS" =16)


#----------------------1.1 pseudo real value-----------------------
pseudo_realsh_net<-rep(0,d_net)
for(i in 1:n_nnet){
  for(j in 1:90){
    pseudo_realsh_net<-pseudo_realsh_net+result_net_srs[(i-1)*90+j,3+1:d_net]*n_net[i]
  }
}
pseudo_realsh_net<-pseudo_realsh_net/sum(n_net)/90
pseudo_realsh_net
as.numeric(round(sort(pseudo_realsh_net,decreasing = TRUE),6))

pseudo_realsh_net2<-rep(0,d_net)
for(i in 1:n_nnet){
  for(j in 1:90){
    pseudo_realsh_net2<-pseudo_realsh_net2+result_net_ls[(i-1)*90+j,3+1:d_net]*n_net[i]
  }
}
pseudo_realsh_net2<-pseudo_realsh_net2/sum(n_net)/90
as.numeric(round(sort(pseudo_realsh_net2,decreasing = TRUE),6))


pseudo_top_net<-sort(as.numeric(pseudo_realsh_net),decreasing = TRUE,index.return=TRUE)$ix
pseudo_top_net


#----------------------1.2 sort---------------------------
#srs
allorder_net_srs<-matrix(0,nrow=nrow(result_net_srs),ncol=d_net)
for(i in 1:nrow(allorder_net_srs)){
  allorder_net_srs[i,]<-sort(as.numeric(result_net_srs[i,3+1:d_net]),decreasing = TRUE,index.return=TRUE)$ix
}
allorder_net_srs


#ls
allorder_net_ls<-matrix(0,nrow=nrow(result_net_ls),ncol=d_net)
for(i in 1:nrow(allorder_net_ls)){
  allorder_net_ls[i,]<-sort(as.numeric(result_net_ls[i,3+1:d_net]),decreasing = TRUE,index.return=TRUE)$ix
}
allorder_net_ls

#----------------------1.3 Figure 5(a)---------------------------
top_net_srs<-top.diffset(t=c(2,5,13,20,24),top=pseudo_top_net,etop=allorder_net_srs)
top_net_ls<-top.diffset(t=c(2,5,13,20,24),top=pseudo_top_net,etop=allorder_net_ls)

net_set_srs<-data.frame(top2=top_net_srs[,1],top5=top_net_srs[,2],top13=top_net_srs[,3],top20=top_net_srs[,4],top24=top_net_srs[,5],size=rep(n_net,each=90),method='SRS')
net_set_ls<-data.frame(top2=top_net_ls[,1],top5=top_net_ls[,2],top13=top_net_ls[,3],top20=top_net_ls[,4],top24=top_net_ls[,5],size=rep(n_net,each=90),method='LS')

net_set<-rbind(net_set_srs,net_set_ls)
net_set$size<-as.character(net_set$size)
net_set_percset_exact<-data.frame(wrong.exactperc(net_set,c(2,5,13,20,24),90),size=rep(n_net,2),method=rep(c('SRS','LS'),each=length(n_net)))


#top 5
ggplot(data=net_set_percset_exact,aes(x=factor(size,level=size_net), y=top5*100,color=factor(method,levels=c('SRS','LS')),group=factor(method,levels=c('SRS','LS')),shape=factor(method,levels=c('SRS','LS')),linetype=factor(method,levels=c('SRS','LS'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=size_net)+
  scale_color_manual(values=cols_net)+scale_shape_manual(values=shape_net)+
  labs(x='Sample size',y='Error rate (%)')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.85,0.81),text=element_text(size = 17))

#----------------------1.4 Figure 5(b)---------------------------
allmse_net<-c()
for(i in 1:nrow(allresult_net)){
  allmse_net<-c(allmse_net,sum((allresult_net[i,3+c(1:d_net)]-pseudo_realsh_net)^2))
}
allmse_net<-data.frame(mse=allmse_net,size=allresult_net$size,rt=allresult_net$rt,method=allresult_net$method)


allmse_net$method<-as.factor(allmse_net$method)
allmse_net$size<-as.factor(allmse_net$size)
ggplot(data=allmse_net,aes(x=factor(size,levels=size_net), y=mse, fill=factor(method,levels=c('SRS','LS')),color=factor(method,levels=c('SRS','LS'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_net)+
  scale_color_manual(values=cols_net)+scale_fill_manual(values=cols_net)+
  labs(x='Sample size',y='Squared loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.81),text=element_text(size = 17))

