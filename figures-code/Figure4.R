library(ggplot2) 
library(patchwork)
library(stringr)

load("simulation-results/simulations_results.RData")
social_cols<- c("LS" = "#00CC66", "SRS" = "#FF6666", "StrRS" = "#AB82FF", "COA" = "#6699FF")


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
p1<-ggplot(data=allresult_so_n19880_top5,aes(x=variable, y=estimate, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
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
p2<-ggplot(data=allmse_so,aes(x=factor(size,levels=size_so), y=mse, fill=factor(method,levels=c('SRS','StrRS','LS','COA')),color=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_so)+
  scale_color_manual(values=social_cols)+scale_fill_manual(values=social_cols)+
  labs(x='Sample size',y='Squared loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 11))

ggsave(p1+p2, file="figures/Figure4.pdf",width=12,height=4)