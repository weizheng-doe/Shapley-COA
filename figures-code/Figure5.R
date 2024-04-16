library(ggplot2) 
library(patchwork)


load("simulation-results/simulations_results.RData")


size_net<-as.factor(n_net)
cols_net<- c("SRS" = "#FF6666","LS" = "#00CC66")
shape_net<-c("SRS" =16,"LS" = 15)
linetype_net<-c("SRS" =1,"LS" = 2)

#----------------------Figure 5(a)---------------------------
top_net_srs<-top.diffset(t=c(2,5,13,20,24),top=pseudo_top_net,etop=allorder_net_srs)
top_net_ls<-top.diffset(t=c(2,5,13,20,24),top=pseudo_top_net,etop=allorder_net_ls)

net_set_srs<-data.frame(top2=top_net_srs[,1],top5=top_net_srs[,2],top13=top_net_srs[,3],top20=top_net_srs[,4],top24=top_net_srs[,5],size=rep(n_net,each=90),method='SRS')
net_set_ls<-data.frame(top2=top_net_ls[,1],top5=top_net_ls[,2],top13=top_net_ls[,3],top20=top_net_ls[,4],top24=top_net_ls[,5],size=rep(n_net,each=90),method='LS')

net_set<-rbind(net_set_srs,net_set_ls)
net_set$size<-as.character(net_set$size)
net_set_percset_exact<-data.frame(wrong.exactperc(net_set,c(2,5,13,20,24),90),size=rep(n_net,2),method=rep(c('SRS','LS'),each=length(n_net)))


#top 5
p1<-ggplot(data=net_set_percset_exact,aes(x=factor(size,level=size_net), y=top5*100,color=factor(method,levels=c('SRS','LS')),group=factor(method,levels=c('SRS','LS')),shape=factor(method,levels=c('SRS','LS')),linetype=factor(method,levels=c('SRS','LS'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=size_net)+
  scale_color_manual(values=cols_net)+scale_shape_manual(values=shape_net)+scale_linetype_manual(values=linetype_net)+
  labs(x='Sample size',y='Error rate (%)')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.85,0.81),text=element_text(size = 17))

#----------------------Figure 5(b)---------------------------
allmse_net<-c()
for(i in 1:nrow(allresult_net)){
  allmse_net<-c(allmse_net,sum((allresult_net[i,3+c(1:d_net)]-pseudo_realsh_net)^2))
}
allmse_net<-data.frame(mse=allmse_net,size=allresult_net$size,rt=allresult_net$rt,method=allresult_net$method)


allmse_net$method<-as.factor(allmse_net$method)
allmse_net$size<-as.factor(allmse_net$size)
p2<-ggplot(data=allmse_net,aes(x=factor(size,levels=size_net), y=mse, fill=factor(method,levels=c('SRS','LS')),color=factor(method,levels=c('SRS','LS'))))+
  geom_boxplot()+
  scale_x_discrete(breaks=size_net)+
  scale_color_manual(values=cols_net)+scale_fill_manual(values=cols_net)+
  labs(x='Sample size',y='Squared loss',fill='Method',color='Method')+
  theme(legend.position=c(0.85,0.81),text=element_text(size = 17))

ggsave(p1+p2, file="figures/Figure5.pdf",width=12,height=4)
