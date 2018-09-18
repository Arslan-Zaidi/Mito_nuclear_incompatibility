library(plyr)
library(ggplot2)
library(ggtern)
library(gridExtra)
library(reshape2)

#read individual id info for admixture results
fam<-read.table('1kg_nam_flipped_pruned.fam',header=F,stringsAsFactors = F)
fam<-fam[,c(1,2)]
colnames(fam)<-c("FID","IID")
#read file linking IDs to pop groups
pop<-read.table('../../Data_tables/pop_1kg.txt',header=T,sep="\t",stringsAsFactors = F)

#link fam file to pop designations
fam<-join(fam,pop,by="IID")

#fill in pop designations for native americans
nat.index<-which(is.na(fam$pop))
fam$pop[nat.index]<-sapply(fam$IID[nat.index],function(x){unlist(strsplit(as.character(x),split="_"))[1]})

#read global ancestry from ADMIXTURE

#write function to read and format *.Q file for different values of K
read.q<-function(k){
  qk<-read.table(paste('1kg_nam_flipped_pruned.',k,'.Q',sep=""),header=F)
  colnames(qk)<-paste("C",seq(1,k),sep="")
  qk$IID<-fam$IID
  qk$pop<-fam$pop
  #continental grouping (African, European, Native American, or admixed) for plotting
  qk$cont_group<-NA
  qk$cont_group[which(qk$pop=="YRI")]<-"AFR"
  qk$cont_group[which(qk$pop=="CEU")]<-"EUR"
  qk$cont_group[which(qk$pop%in%c("MAYAN","AYMARAN","QUECHUAN","NAHUAN"))]<-"NAT"
  qk$cont_group[which(is.na(qk$cont_group)=="TRUE")]<-"Admixed"
  return(qk)
}

q3<-read.q(3)

#read MT haplogroup info
haplo<-read.table('../../Data_tables/AMR_mt.haplo',header=T,sep="\t")
colnames(haplo)[1]<-c("IID")
haplo<-haplo[,c(1,4)]
q3<-join(q3,haplo,by="IID")
q3<-q3[-which(q3$Major_Haplogroup=="M"),]
q3$region_haplogroup<-NA
q3$region_haplogroup[which(q3$Major_Haplogroup%in%c("A","B","C","D"))]<-"NAT"
q3$region_haplogroup[which(q3$Major_Haplogroup%in%c("K","J","U","H","T","V","W","I"))]<-"EUR"
q3$region_haplogroup[which(q3$Major_Haplogroup%in%c("L"))]<-"AFR"


fig_2a<-ggtern()+
  geom_point(data=q3[which(q3$cont_group=="Admixed"),],aes(C1,C2,C3),color="grey",alpha=0.6)+
  geom_point(data=q3[which(q3$cont_group%in%c("NAT","EUR","AFR")),],aes(C1,C2,C3,color=cont_group),size=6,alpha=0.2)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="Continental\ngroup",title="A.")

fig_2b<-ggtern()+
  geom_point(data=q3[which(q3$cont_group=="Admixed"),],aes(C1,C2,C3,color=pop),alpha=0.6)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="Population",x="EUR",y="AFR",z="NAT",title="B.")

fig_2c<-ggtern()+
  geom_point(data=q3,aes(C1,C2,C3,color=Major_Haplogroup),alpha=0.6)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="mtDNA\nhaplogroup",x="EUR",y="AFR",z="NAT",title="C.")

fig_2d<-ggtern()+
  geom_point(data=q3,aes(C1,C2,C3,color=region_haplogroup),alpha=0.6)+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        plot.title=element_text(size=20))+
  labs(color="mtDNA\nhaplogroup",x="EUR",y="AFR",z="NAT",title="D.")

fig_2.combined<-ggtern::grid.arrange(fig_2a,fig_2b,fig_2c,fig_2d,nrow=2,padding=unit(0.01,"line"))

ggsave("ancestry_ternary_plots.pdf",fig_2.combined,height=13,width=13)

##### ADMIXTURE plot for different values of K #######
q2<-read.q(2)
q3<-read.q(3)
q4<-read.q(4)
q5<-read.q(5)

q2$k<-2
q3$k<-3
q4$k<-4
q5$k<-5

#rename columns of q4 and q5 to match the identity of components in q2 and q3
colnames(q4)[c(1:4)]<-c("C3","C4","C2","C1")
colnames(q5)[c(1:5)]<-c("C2","C3","C4","C5","C1")

mq2<-melt(q2,id.vars=c("IID","pop","cont_group","k"),value.name="anc.prop",variable.name="ancestry")
mq3<-melt(q3,id.vars=c("IID","pop","cont_group","k"),value.name="anc.prop",variable.name="ancestry")
mq4<-melt(q4,id.vars=c("IID","pop","cont_group","k"),value.name="anc.prop",variable.name="ancestry")
mq5<-melt(q5,id.vars=c("IID","pop","cont_group","k"),value.name="anc.prop",variable.name="ancestry")


mq.all<-rbind(mq2,mq3,mq4,mq5)

k.plt<-ggplot(mq.all,aes(IID,anc.prop,fill=ancestry))+
  geom_bar(stat="identity",width=2)+
  facet_grid(k~cont_group,scales="free")+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=14))+
  scale_fill_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"))+
  labs(x="Individuals",y="Ancestry fraction",fill="Ancestry component")

ggsave("Admixture_k_plot_09182018.pdf",k.plt,height=7,width=10)
