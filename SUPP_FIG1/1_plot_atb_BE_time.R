library(phyloseq)
library(tidyverse)
library(egg)
library(ggpubr)
library(grid)
load("../DATA_INPUT/EC880F1B.rda")

#MAIN
#rearrange meta data
source("../HARMONIZE_FIGURES/1_plot_utils.R")
met<-data.frame(sample_data(EC880))

#subset
s<-data.frame(met[,c(1:3,8,23:27)])
s$SampleID<-rownames(met)

#Coli, Doxy, Tobr fused to group Other
#remember sample IDs
meC<-which(s$Coli_In!=0)
meD<-which(s$Doxy_Or!=0)
#replace OThers with value
s$Others<-s$Tobr_In
s$Others[meC]<-s$Coli_In[meC]
s$Others[meD]<-s$Doxy_Or[meD]

#Clear names
s$treat<-NA
s$treat[s$Azit_Or!=0]<-"Azithromycin"
s$treat[s$Aztr_In!=0]<-"Aztreonam"
s$treat[s$Others!=0]<-"Others"
s$treat<-factor(s$treat,levels=c("Azithromycin", "Aztreonam", "Others"))
s$BETRU<-ifelse(s$BETRU==1,"Baseline","Exacerbation")

#order ECIds by patientA-K
me<-s[s$ECid=="10",]
me[,"ECid"]<-"3"
s<-s[s$ECid!="10",]
s[s$ECid=="9","ECid"]<-"10"
s[s$ECid=="8","ECid"]<-"9"
s[s$ECid=="7","ECid"]<-"8"
s[s$ECid=="6","ECid"]<-"7"
s[s$ECid=="5","ECid"]<-"6"
s[s$ECid=="4","ECid"]<-"5"
s[s$ECid=="3","ECid"]<-"4"
s<-rbind(s,me)

#group by time series
ss <- s %>% group_by(ECid)
ss$ECid<-factor(ss$ECid, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))

#subset data for common legend
d1<-ss[,c(1,3:5,12)]
d11<-d1[!is.na(d1$treat),]
d2<-ss[,c(1,3:4,7,12)]
d3<-ss[,c(1, 3:4,11,12)]
colnames(d1)[4]<-"Maintenance"
colnames(d11)[4]<-"Maintenance"	
colnames(d2)[4]<-"Maintenance"	
colnames(d3)[4]<-"Maintenance"	

####
# FIGURE SI 1B
####

#common legends
p11<-ggplot(d11, aes(-SampleDay2T, ECid, color=factor(Maintenance), fill=treat, shape=factor(BETRU))) + geom_point(aes(alpha=factor(Maintenance), shape=factor(BETRU))) + theme_bw()+ggtitle("Daily maintenance antibiotic therapy")+scale_color_manual(values=c("#C0C0C0", "#56B4E9"))+scale_alpha_manual(values=c(0.3, 1))+guides(color="none")+guides(alpha="none")+guides(size="none")+guides(shape=guide_legend(title="Clinical state", override.aes = list(size=2)))+ guides(fill=guide_legend(title="Maintenance\nantibiotics",override.aes = list(col = c( "#56B4E9","#CC3366", "#E69F00"), shape=c(4,4,4))))+ylab("")+theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +  theme(axis.ticks.length.x = unit(0, "cm")) + theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) +scale_shape_manual(values=c(4,17))
p1<-ggplot(d1, aes(-SampleDay2T, ECid, color=factor(Maintenance), shape=factor(BETRU))) + geom_point(aes(alpha=factor(Maintenance))) + theme_bw()+ggtitle("Daily maintenance antibiotic therapy")+scale_color_manual(values=c("#C0C0C0", "#56B4E9"))+scale_alpha_manual(values=c(0.3, 1))+ guides(alpha="none")+guides(color="none")+ylab("")+theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +  theme(axis.ticks.length.x = unit(0, "cm")) + theme(plot.margin = unit(c(0,0.5,0,0), "cm"))+scale_shape_manual(values=c(4,17))
p2<-ggplot(d2, aes(-SampleDay2T, ECid, color=factor(Maintenance), shape=factor(BETRU))) + geom_point(aes(alpha=factor(Maintenance))) + theme_bw()+scale_color_manual(values=c("#C0C0C0", "#CC3366"))+scale_alpha_manual(values=c(0.3, 1))+ guides(alpha="none")+guides(color="none")+ylab("")+theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +  theme(axis.ticks.length.x = unit(0, "cm")) + theme(plot.margin = unit(c(0,0.5,0,0), "cm"))+scale_shape_manual(values=c(4,17))
p3<-ggplot(d3, aes(-SampleDay2T, ECid, color=factor(Maintenance), shape=factor(BETRU))) + geom_point(aes(alpha=factor(Maintenance))) + theme_bw()+scale_color_manual(values=c("#C0C0C0", "#E69F00"))+scale_alpha_manual(values=c(0.3, 1),labels=c("free","used"))+ guides(alpha="none")+guides(color="none") +ylab("")+scale_x_continuous(breaks=seq(-60,0, by=10))+xlab("Days to start of PEx therapy")+theme(plot.margin = unit(c(0,0.5,0,0), "cm"))+scale_shape_manual(values=c(4,17))

#legend
f1<-get_legend(p11)
#arrangement
ff<-ggarrange(p1,p2, p3, nrow = 3, legend="none")
ff<-annotate_figure(ff, left = textGrob("PEx time series", rot = 90, vjust = 1.5, gp = gpar(fontsize= 11)))

####
#FIGURE SI 1A
####
#patients vs time lines
#paper color code
coh<-paper_pal("Paper")(11)
names(coh)<-levels(as.factor((as(sample_data(EC880),"data.frame"))$Patient))

#plot
p4<-ggplot(ss, aes(-SampleDay2T, ECid, color=factor(Patient))) + geom_point() + geom_hline(aes(yintercept=3.5))+ geom_hline(aes(yintercept=4.5))+ geom_hline(aes(yintercept=5.5))+ geom_hline(aes(yintercept=6.5))+ geom_hline(aes(yintercept=9.5))+ geom_hline(aes(yintercept=10.5))+ geom_hline(aes(yintercept=11.5))+ geom_hline(aes(yintercept=12.5))+ geom_hline(aes(yintercept=15.5))+ geom_hline(aes(yintercept=16.5))+ggtitle("Sputum samples")+scale_colour_manual(values=coh,labels=c("Subject A","Subject B","Subject C","Subject D","Subject E", "Subject F","Subject G","Subject H","Subject I","Subject J", "Subject K"),name="Cohort")+ theme_bw()+ylab("PEx time series")+xlab("Days to start of PEx therapy")+scale_x_continuous(breaks=seq(-60,0, by=10))

#legend
l1<-get_legend(p4)

####
#arrangement SI Figure1
lg<-plot_grid(l1,f1,axis="t",align="h", ncol=2)
p5<-plot_grid(p4+theme(legend.position="none"), lg,nrow=2, labels=c("A",""))
SIX<-plot_grid(p5,ff+theme(legend.position="none"), rel_widths=c(1,1.14), labels=c("","B"))

save_plot(SIX, file="SI_figure1.pdf", base_height=6.5, base_asp=1.1)
save_plot(SIX, file="SI_figure1.png", base_height=6.5, base_asp=1.1)

