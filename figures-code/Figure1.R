library(ggplot2)
library(patchwork)

var_shape1<-c("SRS" =16,"COA or LS" = 15)
d_ex1<-8;m_ex1<-d_ex1*(d_ex1-1)*c(1:10)
var1_ex1<-(d_ex1-1)/m_ex1/d_ex1^2
var2_ex1<-rep(0,length(var1_ex1))
var1_ex1<-data.frame(var=var1_ex1,m=m_ex1,method='SRS')
var2_ex1<-data.frame(var=var2_ex1,m=m_ex1,method='COA or LS')

m_ex1_g<-seq(min(m_ex1),max(m_ex1),1)
var_ex1<-rbind(var1_ex1,var2_ex1)
var_ex1$m<-as.character(var_ex1$m)



var_shape2<-c("SRS or LS" =16, "COA" = 15)
beta2<-2; beta1<-1
sigma21<--0.5; sigma2<-0.8
g1<-beta2^2*(sigma2-sigma21/sigma2)
g2<-(sigma21*beta1+sigma2*beta2)^2/sigma2
d_ex2<-11;m_ex2<-d_ex2*(d_ex2-1)*c(1:10)
var1_ex2<-(g1-g2)^2/4/m_ex2
var2_ex2<-rep(0,length(var1_ex2))
var1_ex2<-data.frame(var=var1_ex2,m=m_ex2,method='SRS or LS')
var2_ex2<-data.frame(var=var2_ex2,m=m_ex2,method='COA')

m_ex2_g<-seq(min(m_ex2),max(m_ex2),1)
var_ex2<-rbind(var1_ex2,var2_ex2)
var_ex2$m<-as.character(var_ex2$m)

#----------------Figure 1(a)---------------------
p1<-ggplot(data=var_ex1,aes(x=factor(m,level=m_ex1_g), y=sqrt(abs(var)),color=factor(method,levels=c('SRS','COA or LS')),group=factor(method,levels=c('SRS','COA or LS')),shape=factor(method,levels=c('SRS','COA or LS')),linetype=factor(method,levels=c('SRS','COA or LS'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=m_ex1,expand = expansion(0.05))+scale_shape_manual(values=var_shape1)+
  labs(x='Sample size',y='Standard deviation of estimates')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.83,0.78),text=element_text(size = 17))


#----------------Figure 1(b)---------------------
p2<-ggplot(data=var_ex2,aes(x=factor(m,level=m_ex2_g), y=sqrt(abs(var)),color=factor(method,levels=c('SRS or LS','COA')),group=factor(method,levels=c('SRS or LS','COA')),shape=factor(method,levels=c('SRS or LS','COA')),linetype=factor(method,levels=c('SRS or LS','COA'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=m_ex2,expand = expansion(0.05))+scale_shape_manual(values=var_shape2)+
  labs(x='Sample size',y='Standard deviation of estimates')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.83,0.78),text=element_text(size = 17))

ggsave(p1+p2, file="figures/Figure1.pdf",width=10,height=4)