#########################################################
#################Code for Examples 1 and 3###############
#########################################################

library(ggplot2)

var_shape1<-c("SRS" =16,"COA or LS" = 15)
d_ex1<-8;m_ex1<-d_ex1*(d_ex1-1)*c(1:10)
var1_ex1<-(d_ex1-1)/m_ex1/d_ex1^2
var2_ex1<-rep(0,length(var1_ex1))
var1_ex1<-data.frame(var=var1_ex1,m=m_ex1,method='SRS')
var2_ex1<-data.frame(var=var2_ex1,m=m_ex1,method='COA or LS')

m_ex1_g<-seq(min(m_ex1),max(m_ex1),1)
var_ex1<-rbind(var1_ex1,var2_ex1)
var_ex1$m<-as.character(var_ex1$m)

#----------------Figure 1(a)---------------------
ggplot(data=var_ex1,aes(x=factor(m,level=m_ex1_g), y=sqrt(abs(var)),color=factor(method,levels=c('SRS','COA or LS')),group=factor(method,levels=c('SRS','COA or LS')),shape=factor(method,levels=c('SRS','COA or LS')),linetype=factor(method,levels=c('SRS','COA or LS'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=m_ex1,expand = expansion(0.05))+scale_shape_manual(values=var_shape1)+
  labs(x='Sample size',y='Standard deviation of estimates')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.83,0.78),text=element_text(size = 17))
