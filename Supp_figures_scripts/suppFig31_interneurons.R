library(ggplot2)
library(patchwork)
library(stringr)
library(ggpubr)
library(rstatix)
input_path='/scratch/users/banwang3/celllines/scripts/RScripts/figscripts/suppfigscripts/alexdata/interneurons/'
output_path='/scratch/users/banwang3/celllines/scripts/RScripts/Rplots/'
alldf =read.csv(paste0(input_path,'All_Interneurons_Pseudo.txt'),sep='\t')
alldf$type='All'
alldf$bc=str_split_fixed(alldf$X,'-',2)[,1]
CGE=read.csv(paste0(input_path,'CGE_Interneurons_Pseudo.txt'),sep='\t')
CGE$type='CGE'
CGE$bc=str_split_fixed(CGE$X,'-',2)[,1]
LGE=read.csv(paste0(input_path,'LGE_Interneurons_Pseudo.txt'),sep='\t')
LGE$type='LGE'
LGE$bc=str_split_fixed(LGE$X,'-',2)[,1] 
MGE=read.csv(paste0(input_path,'MGE_Interneurons_Pseudo.txt'),sep='\t')
MGE$type='MGE'
MGE$bc=str_split_fixed(MGE$X,'-',2)[,1]   
humancell=read.csv(paste0(input_path,'cell_to_species_human.txt'),sep='\t')
chimpcell=read.csv(paste0(input_path,'cell_to_species_chimp.txt'),sep='\t')
celltable=rbind(humancell,chimpcell)
colnames(celltable)=c('bc','Species')
celltable$Species=str_to_title(celltable$Species)
alldf=merge(alldf,celltable,all.x=T)
alldf$Species=factor(alldf$Species,levels=c('Human','Chimp'))
CGE=merge(CGE,celltable,all.x=T)
CGE$Species=factor(CGE$Species,levels=c('Human','Chimp'))
MGE=merge(MGE,celltable,all.x=T)
MGE$Species=factor(MGE$Species,levels=c('Human','Chimp'))
LGE=merge(LGE,celltable,all.x=T)
LGE$Species=factor(LGE$Species,levels=c('Human','Chimp'))
commoncells=intersect(CGE$bc,LGE$bc)
commoncells=intersect(commoncells,MGE$bc)
df = rbind(alldf,CGE,MGE,LGE)
df$progenitor=ifelse(df$bc %in% commoncells,'Progenitor','Non-progenitor')
binpseudo_equalbinsize<-function(df,filter0=T,fixonprogenitor=T){
    if (filter0){  
        df = df[df$GAD1_Normed_Counts>0,]
    }
    types=unique(df$type)
    df$bin=0
    df = df[order(df$Pseudotime),]
    final=c()
    for(type in types){
        subdf = df[df$type==type,]
        if (fixonprogenitor){
            subdf1=subdf[subdf$progenitor=='Progenitor',]
            subdf2=subdf[subdf$progenitor!='Progenitor',]
            steps = floor(nrow(subdf1)/2)
            subdf1[0:steps,'bin']=1
            subdf1[(steps+1):nrow(subdf1),'bin']=2
            steps = floor(nrow(subdf2)/3)
            subdf2[0:steps,'bin']=3
            subdf2[(steps+1):(steps*2),'bin']=4
            subdf2[(2*steps+1):nrow(subdf2),'bin']=5
            final=rbind(final,subdf1,subdf2)
        }else{
            steps = floor(nrow(subdf)/5)
            subdf[0:steps,'bin']=1
            subdf[(steps+1):(steps*2),'bin']=2
            subdf[(2*steps+1):(steps*3),'bin']=3
            subdf[(3*steps+1):(steps*4),'bin']=4
            subdf[(4*steps+1):(nrow(subdf)),'bin']=5
            final=rbind(final,subdf)
        }
    }
    final$bin=factor(final$bin)
    final$type=factor(final$type)
    return(final)
}
fixonprogenitor=T
df = binpseudo_equalbinsize(df,fixonprogenitor=fixonprogenitor) 
stat.test <- df %>% group_by(type,bin) %>% wilcox_test(GAD1_Normed_Counts ~ Species) %>%add_xy_position(x = "bin", dodge = 0.8) 
stat.test$plab<-ifelse(stat.test$p<0.005,'***',ifelse(stat.test$p<0.01,'**',ifelse(stat.test$p<0.05,'*','ns')))

#p<-ggplot(df,aes(x=bin,y=GAD1_Normed_Counts,fill=Species))+
#    geom_boxplot(outlier.size = 0.25,facet.by="type")+
p<-ggboxplot(df,x="bin",y="GAD1_Normed_Counts",fill="Species",outlier.size=0.25,facet.by="type",nrow=4)+    
    ylab('GAD1 Normed Counts')+
    xlab('Pseudotime bin by equal size')+
    theme_bw()+
    scale_fill_manual(values=c("#FF6A6A","#1E90FF"))+
    stat_pvalue_manual(stat.test,label = "p", tip.length = 0.02)+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))+
    theme(legend.position="top",panel.grid.major=element_blank(),panel.grid.minor=element_blank())

pdf(paste0(output_path,'interneuron_binpseudo_equalsizebin_fixonprog_',fixonprogenitor,'.pdf'),width=6,height=12,useDingbats=FALSE)
print(p)
dev.off()
for (type in c('All','lineage')){
#for(type in unique(df$type)){
    for (fixonprogenitor in c(T,F)){
    #p<-binpseudo_fixvalue(df[df$type==type,])
    #pdf(paste0(output_path,type,'_interneuron_binpseudo_fixvalue.pdf'),width=6,height=4,useDingbats=FALSE)
    #print(p)
    #dev.off()
    if (type=='All'){
        p<-binpseudo_equalbinsize(df[df$type=='All',],fixonprogenitor=fixonprogenitor)
        pdf(paste0(output_path,type,'_interneuron_binpseudo_equalsizebin_fixonprog_',fixonprogenitor,'.pdf'),width=6,height=4,useDingbats=FALSE)
    print(p)
    dev.off()
    }else{
        p<-binpseudo_equalbinsize(df[df$type!='All',],fixonprogenitor=fixonprogenitor)
        pdf(paste0(output_path,type,'_interneuron_binpseudo_equalsizebin_fixonprog_',fixonprogenitor,'.pdf'),width=6,heigh=10,useDingbats=FALSE)
    print(p)
    dev.off()
}
        
    
}
    }
