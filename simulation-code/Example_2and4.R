#########################################################
#################Code for Examples 2 and 4###############
#########################################################

library(ggplot2)

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

#----------------Figure 1(b)---------------------
ggplot(data=var_ex2,aes(x=factor(m,level=m_ex2_g), y=sqrt(abs(var)),color=factor(method,levels=c('SRS or LS','COA')),group=factor(method,levels=c('SRS or LS','COA')),shape=factor(method,levels=c('SRS or LS','COA')),linetype=factor(method,levels=c('SRS or LS','COA'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=m_ex2,expand = expansion(0.05))+scale_shape_manual(values=var_shape2)+
  labs(x='Sample size',y='Standard deviation of estimates')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.83,0.78),text=element_text(size = 17))
