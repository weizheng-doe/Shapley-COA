library(ggplot2) 
library(patchwork)


load("simulation-results/simulations_results.RData")


social_cols<- c("LS" = "#00CC66", "SRS" = "#FF6666", "StrRS" = "#AB82FF", "COA" = "#6699FF")
social_shape<-c("LS" = 15, "SRS" =16, "StrRS" = 17, "COA" =18)
social_linetype<-c("LS" = 2, "SRS" =1, "StrRS" = 3, "COA" =2)

top_so_srs<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_srs)
top_so_strrs<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_strrs)
top_so_ls<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_ls)
top_so_coa<-top.diffset(t=c(5,9,10,12,15,20,24),top=pseudo_top_source,etop=allorder_so_coa)
social_so_set_srs<-data.frame(top5=top_so_srs[,1],top9=top_so_srs[,2],top10=top_so_srs[,3],top12=top_so_srs[,4],top15=top_so_srs[,5],top20=top_so_srs[,6],top24=top_so_srs[,7],size=rep(nso_4methodnear71_69_modified,each=300),method='SRS')
social_so_set_strrs<-data.frame(top5=top_so_strrs[,1],top9=top_so_strrs[,2],top10=top_so_strrs[,3],top12=top_so_strrs[,4],top15=top_so_strrs[,5],top20=top_so_strrs[,6],top24=top_so_strrs[,7],size=rep(nso_4methodnear71_69_modified,each=300),method='StrRS')
social_so_set_ls<-data.frame(top5=top_so_ls[,1],top9=top_so_ls[,2],top10=top_so_ls[,3],top12=top_so_ls[,4],top15=top_so_ls[,5],top20=top_so_ls[,6],top24=top_so_ls[,7],size=rep(nso_4methodnear71_69_modified,each=300),method='LS')
social_so_set_coa<-data.frame(top5=top_so_coa[,1],top9=top_so_coa[,2],top10=top_so_coa[,3],top12=top_so_coa[,4],top15=top_so_coa[,5],top20=top_so_coa[,6],top24=top_so_coa[,7],size=rep(nso_4methodnear71_coa_modified,each=300),method='COA')


social_so_set_1<-rbind(social_so_set_srs,social_so_set_strrs,social_so_set_ls)
social_so_set_1$size<-as.character(social_so_set_1$size)
social_so_set_coa$size<-as.character(social_so_set_coa$size)


social_set<-rbind(social_so_set_1,social_so_set_coa)
social_set$size<-as.character(social_set$size)

social_so_percset_exact_1<-data.frame(wrong.exactperc(social_so_set_1,c(5,9,10,12,15,20,24),300),size=rep(nso_4methodnear71_69_modified,3),method=rep(c('SRS','StrRS','LS'),each=length(nso_4methodnear71_69_modified)))
social_so_percset_exact_coa<-data.frame(wrong.exactperc(social_so_set_coa,c(5,9,10,12,15,20,24),300),size=nso_4methodnear71_coa_modified,method='COA')


social_so_percset_exact<-rbind(social_so_percset_exact_coa,social_so_percset_exact_1)

social_so_percset_exact[social_so_percset_exact$size==4968,]$size<-'4968(4970)'
social_so_percset_exact[social_so_percset_exact$size==4970,]$size<-'4968(4970)'
social_so_percset_exact[social_so_percset_exact$size==c(9936),]$size<-'9936(9940)'
social_so_percset_exact[social_so_percset_exact$size==c(9940),]$size<-'9936(9940)'
social_so_percset_exact[social_so_percset_exact$size==c(14904),]$size<-'14904(14910)'
social_so_percset_exact[social_so_percset_exact$size==c(14910),]$size<-'14904(14910)'
social_so_percset_exact[social_so_percset_exact$size==c(19872),]$size<-'19872(19880)'
social_so_percset_exact[social_so_percset_exact$size==c(19880),]$size<-'19872(19880)'
social_so_percset_exact[social_so_percset_exact$size==c(24840),]$size<-'24840(24850)'
social_so_percset_exact[social_so_percset_exact$size==c(24850),]$size<-'24840(24850)'

social_so_percset_exact$size<-as.factor(social_so_percset_exact$size)
levels(social_so_percset_exact$size)

social_so_percset_exact$method<-factor(social_so_percset_exact$method,levels = 'SRS','StrRS','LS','COA')

#top 12
p1<-ggplot(data=social_so_percset_exact,aes(x=factor(size,level=size_so), y=top12*100,color=factor(method,levels=c('SRS','StrRS','LS','COA')),group=factor(method,levels=c('SRS','StrRS','LS','COA')),shape=factor(method,levels=c('SRS','StrRS','LS','COA')),linetype=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=size_so)+
  scale_color_manual(values=social_cols)+scale_shape_manual(values=social_shape)+scale_linetype_manual(values=social_linetype)+
  labs(x='Sample size',y='Error rate(%)')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.88,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 11))

#top 20
p2<-ggplot(data=social_so_percset_exact,aes(x=factor(size,level=size_so), y=top20*100,color=factor(method,levels=c('SRS','StrRS','LS','COA')),group=factor(method,levels=c('SRS','StrRS','LS','COA')),shape=factor(method,levels=c('SRS','StrRS','LS','COA')),linetype=factor(method,levels=c('SRS','StrRS','LS','COA'))))+
  geom_point()+geom_line()+
  scale_x_discrete(breaks=size_so)+
  scale_color_manual(values=social_cols)+scale_shape_manual(values=social_shape)+scale_linetype_manual(values=social_linetype)+
  labs(x='Sample size',y='Error rate(%)')+labs(color='Method')+labs(shape='Method')+labs(linetype='Method')+
  theme(legend.position=c(0.88,0.75),text=element_text(size = 17),axis.text.x = element_text(size = 11))

ggsave(p1+p2, file="figures/Figure3.pdf",width=12,height=4)