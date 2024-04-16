library(ggplot2) 
library(patchwork)


load("simulation-results/simulations_results.RData")


d<-5
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
p1<-ggplot(data=allest_msmd5_n20,aes(x=variable, y=estimate, fill=method,color=method))+
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

p2<-ggplot(data=allresult_msmd5_partial_clear,aes(x=size, y=d*mse, fill=method,color=method))+
  geom_boxplot()+
  scale_x_discrete(labels=n2)+
  labs(x='Sample size',y='Square loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.8,0.8))

ggsave(p1+p2, file="figures/FigureS3.pdf",width=12,height=4)