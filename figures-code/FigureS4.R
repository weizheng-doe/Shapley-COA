library(ggplot2) 
library(patchwork)
library(stringr)

load("simulation-results/simulations_results.RData")

d<-7
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
pa<-ggplot(data=allest_bs_n42,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=c('Pclass','Sex','Age','SibSp','Parch','Fare','Embarked'))+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.13,0.78))
#(b)
pb<-ggplot(data=allresult_bs,aes(x=size, y=d*mse, fill=method,color=method))+
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
pc<-ggplot(data=allest_ces_n42,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=c('Pclass','Sex','Age','SibSp','Parch','Fare','Embarked'))+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.13,0.78))
#(d)
pd<-ggplot(data=allresult_ces,aes(x=size, y=d*mse, fill=method,color=method))+
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
pe<-ggplot(data=allest_cs_n42,aes(x=variable, y=estimate, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels=c('Pclass','Sex','Age','SibSp','Parch','Fare','Embarked'))+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.13,0.78))
#(f)
pf<-ggplot(data=allresult_cs,aes(x=size, y=d*mse, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(labels=n2)+
  labs(x='Sample size',y='Square loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))



ggsave(pa+pb+pc+pd+pe+pf+ plot_layout(ncol = 2), file="figures/FigureS4.pdf",width=14,height=4*3)
