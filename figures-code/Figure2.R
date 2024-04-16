library(ggplot2) 
library(patchwork)


load("simulation-results/simulations_results.RData")


allresult_air<-rbind(result_air_srs,result_air_strrs,result_air_ls,result_air_coa)
size_air<-as.factor(c('10100','20200','30300','40400','50500','60600','70700','80800','90900','101000'))
cols_air<- c("LS" = "#00CC66", "SRS" = "#FF6666", "StrRS" = "#AB82FF", "COA" = "#6699FF")


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
allest_air_n10100_partial$variable<-factor(allest_air_n10100_partial$variable,levels = c('x1','x9','x19','x26','x39','x51','x62','x72','x87','x97'))

p1<-ggplot(data=allest_air_n10100_partial,aes(x=variable, y=estimate, fill=factor(method,level=c('SRS','StrRS','LS','COA')),color=factor(method,level=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=c('x1','x9','x19','x26','x39','x51','x62','x72','x87','x97'),labels=c('x1','x9','x19','x26','x39','x51','x62','x72','x87','x97'))+
  scale_color_manual(values=cols_air)+scale_fill_manual(values=cols_air)+
  labs(x='Variable',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.2,0.75),text = element_text(size = 17),axis.text.x = element_text(size = 10))

#---------------------1.2 Figure 2(b): square loss------------------------
allresult_air$size<-as.character(allresult_air$size)
p2<-ggplot(data=allresult_air,aes(x=factor(size,level=size_air), y=d_air*mse, fill=factor(method,level=c('SRS','StrRS','LS','COA')),color=factor(method,level=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_air)+
  scale_color_manual(values=cols_air)+scale_fill_manual(values=cols_air)+
  labs(x='Sample size',y='Squared loss')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.85,0.7),text = element_text(size = 17),axis.text.x = element_text(size = 10))

ggsave(p1+p2, file="figures/Figure2.pdf",width=10,height=4)