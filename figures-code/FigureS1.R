library(ggplot2) 
library(patchwork)
library(stringr)

load("simulation-results/simulations_results.RData")

#---------------------1. summarize------------------------
cols_sspi<- c("SRS" = "#FF6666", "StrRS" = "#AB82FF", "LS" = "#00CC66", "COA" = "#6699FF")
shape_sspi<-c("SRS" =16, "StrRS" = 17, "LS" = 15, "COA" =18)

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
pb<-ggplot(data=allmse_sspi_partial_1_clear,aes(x=factor(size,levels=size_sspi), y=mse, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_sspi,labels = function(x) str_wrap(x, width = 1)) +
  scale_color_manual(values=cols_sspi)+scale_fill_manual(values=cols_sspi)+
  labs(x='Sample size',y='Square loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 12))

#Figure S1.(c)
pc<-ggplot(data=allmse_sspi_partial_2,aes(x=factor(size,levels=size_sspi), y=mse, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
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

sspi_names_xx<-as.vector(sspi_names[,2])
sspi_names_xx[8]<-"DC"

pa<-ggplot(data=allest_sspi_n2805,aes(x=variable, y=estimate, fill=factor(method,level=c('SRS','StrRS','LS','COA')),color=factor(method,level=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=tnames,labels = sspi_names_xx)+
  scale_color_manual(values=cols_sspi)+scale_fill_manual(values=cols_sspi)+
  labs(x='State',y='Estimated Shapley value')+labs(color='Method')+labs(fill='Method')+
  theme(legend.position=c(0.95,0.8),text = element_text(size = 17),axis.text.x = element_text(size = 10))

ggsave(pa/(pb|pc), file="figures/FigureS1.pdf",width=14,height=7)
