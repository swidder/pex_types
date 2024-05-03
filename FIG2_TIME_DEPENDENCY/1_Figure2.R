library(phyloseq)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(effectsize)
library(robCompositions)
library(dplyr)
library(grid)
source("../HARMONIZE_FIGURES/1_plot_utils.R")


###########
#LOAD DATA
###########
load("../DATA_INPUT/EC880F1B.rda")
dat<-as(sample_data(EC880),"data.frame")
dat<-dat[!is.na(dat$PExClust),]

########
#FIGURE2
########

###########
#DMM TIME EVOLUTION 2A
############
#harmonize color code
dmm<-paper_pal("DMM")(6)
names(dmm)<-c("1","2","3","4","5","6")
clear<-c("1-pathogens", "2-fac. anaerobes", "3-pathogens","4-anaerobes", "5-anaerobes", "6-anaerobes")
names(clear)<-names(dmm)

#keep original patient/time series color code
coh<-paper_pal("Paper")(11)
names(coh)<-levels(as.factor((as(sample_data(EC880),"data.frame"))$Patient))
coh1<-coh[!names(coh)%in%"H"]
ts<-paper_pal("Age")(16)

#time series ids ordered by patients
dat[dat$ECid=="10","ECid"]<-"X"
dat[dat$ECid=="9","ECid"]<-"10"
dat[dat$ECid=="8","ECid"]<-"9"
dat[dat$ECid=="7","ECid"]<-"8"
dat[dat$ECid=="6","ECid"]<-"7"
dat[dat$ECid=="5","ECid"]<-"6"
dat[dat$ECid=="4","ECid"]<-"5"
dat[dat$ECid=="3","ECid"]<-"4"
dat[dat$ECid=="X","ECid"]<-"3"
#order levels numerically
dat$ECid<-factor(dat$ECid, levels=c("1","2","3","4","5","6","7","8","10","11","13","14","15","16","17","18"))
names(ts)<-levels(as.factor(dat$ECid))

#harmonize order of time groups
my_comparisons <- list( c("<24", "24-60") )
dat$D23<-factor(dat$D23, levels=c("<24", "24-60"))
xx<-dat$class
dat$classN<-ifelse(xx=="1-5","-1", ifelse(xx=="6-10", "-6", ifelse(xx=="11-15", "-11", ifelse(xx=="16-20", "-16", ifelse(xx=="21-25", "-21", ifelse(xx=="26-30", "-26", ifelse(xx=="31-35", "-31", ifelse(xx=="36-40", "-36", ifelse(xx=="41-45", "-41", ifelse(xx=="46-50", "-46", ifelse(xx=="51-55", "-51", "-56")))))))))))
dat$classN<-factor(dat$classN,levels=c("-56","-51","-46","-41","-36","-31","-26","-21","-16","-11","-6","-1"))

#plot/calculate by PEx cluster
for (n in 1: 3){
if(n==1){
#cluster 1
me<-dat[dat$PExClust%in%"PAT",]
lab="PAT"
tit<-"Community class evolution"
tit2<-"Samples by subject"
tit3<-"Samples by time series"
tit4<-"Diversity in time groups"
tit5<-"Richness in time groups"
yy<-"#F1BB7B"
}
if(n==2){
#cluster 2
me<-dat[dat$PExClust%in%"AN1",]
lab="AN1"
tit5<-tit4<-tit3<-tit2<-tit<-""
yy<-"#FD6467"
}
if(n==3){
#cluster 3
me<-dat[dat$PExClust%in%"AN2",]
lab="AN2"
tit5<-tit4<-tit3<-tit2<-tit<-""
yy<-"#5B1A18"
}

cc<-me %>% count(classN, DMM)
cc$title<-lab
p<-ggplot(cc,aes(x=as.factor(classN),y=n,fill=DMM)) + geom_bar(position="stack", stat="identity") + meg_stef_theme()+ theme(axis.text.x = element_text(angle = 45, vjust=0.5)) + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_fill_manual(values=dmm,name="Community\nclasses",labels=clear) + ggtitle(tit) + xlab("Days to start of PEx therapy") + ylab("Sample number")+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+theme(legend.key.size = unit(0.5, "cm"))+ylim(c(0,25))

cc<-me %>% count(classN, Patient)
cc$title<-lab
sp<-ggplot(cc,aes(x=as.factor(classN),y=n,fill=Patient)) + geom_bar(position="stack", stat="identity") + meg_stef_theme() + theme(axis.text.x = element_text(angle = 45, vjust=0.5))+ guides(color=guide_legend(override.aes=list(fill=NA))) + scale_fill_manual(values=coh1,name="Cohort") + ggtitle(tit2) + xlab("Days to start of PEx therapy") + ylab("Sample number")+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+theme(legend.key.size = unit(0.5, "cm"))+ylim(c(0,25))

cc<-me %>% count(classN, ECid)
cc$title<-lab
sec<-ggplot(cc,aes(x=as.factor(classN),y=n,fill=ECid)) + geom_bar(position="stack", stat="identity") + meg_stef_theme() + theme(axis.text.x = element_text(angle = 45, vjust=0.5))+ guides(color=guide_legend(override.aes=list(fill=NA))) + scale_fill_manual(values=ts,name="Time series")+ ggtitle(tit3) + xlab("Days to start of PEx therapy") + ylab("Sample number")+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+theme(legend.key.size = unit(0.5, "cm"))+ylim(c(0,25))

cc<-me %>% count(Shannon, D23)
cc$title<-lab
sw<-ggplot(cc,aes(x=D23,y=Shannon,col=D23))+geom_boxplot(outlier.alpha=0., outlier.size=0.5, size=0.3) +stat_compare_means(comparisons = my_comparisons,label = "p.format")+ meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_color_paper(discrete=T,palette="Time",name="Days to\ntreatment")+ ggtitle(tit4) + xlab("Time to treatment") + ylab("Shannon")+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+theme(legend.key.size = unit(0.5, "cm"))+ylim(c(0,4))+geom_jitter(aes(col=D23),shape=16, position=position_jitter(0.2), size=0.7, alpha=0.5) 

cc<-me %>% count(Chao1, D23)
cc$title<-lab
sr<-ggplot(cc,aes(x=D23,y=Chao1,col=D23))+geom_boxplot(outlier.alpha=0., outlier.size=0.5, size=0.3) +stat_compare_means(comparisons = my_comparisons,label = "p.format")+ meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_color_paper(discrete=T,palette="Time",name="Days to\ntreatment")+ ggtitle(tit5) + xlab("Time to treatment") + ylab("Chao1")+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+theme(legend.key.size = unit(0.5, "cm"))+ylim(c(NA,max(cc$Chao1)+20))+geom_jitter(aes(col=D23),shape=16, position=position_jitter(0.2), size=0.7, alpha=0.5) 


if(n==1){
d1<-p
sp1<-sp
sec1<-sec
sw1<-sw
sr1<-sr
}
if(n==2){
d2<-p
sp2<-sp
sec2<-sec
sw2<-sw
sr2<-sr
}
if(n==3){
d3<-p
sp3<-sp
sec3<-sec
sw3<-sw
sr3<-sr
}}

#legend 2 DMMs color profile
cc<-dat %>% count(class, DMM)
dlegend<-get_legend(ggplot(cc,aes(x=as.factor(class),y=n,fill=DMM)) + geom_bar(position="stack", stat="identity") + meg_stef_theme() + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_fill_manual(values=dmm,name="Community\nclasses",labels=clear)) 
#legend supp 2C subjects color profile
cc<-dat %>% count(class, Patient)
sl1<-get_legend(ggplot(cc,aes(x=as.factor(class),y=n,fill=Patient))+ meg_stef_theme()+ geom_bar(position="stack", stat="identity")+scale_fill_discrete(guide = "none")+ scale_fill_manual(values=coh1,labels=c("Subject A","Subject B","Subject C","Subject D","Subject E", "Subject F","Subject G","Subject I","Subject J", "Subject K"),name="Cohort")+guides(fill=guide_legend(ncol=2)))
#legend supp2D time series color profile
cc<-dat %>% count(class, ECid)
sl2<-get_legend(ggplot(cc,aes(x=as.factor(class),y=n,fill=ECid))+ meg_stef_theme()+ geom_bar(position="stack", stat="identity")+scale_fill_discrete(guide = "none")+ scale_fill_manual(values=ts,name="Time series")+guides(fill=guide_legend(ncol=3)))


#############
# SHANNON TIME EVOLUTION
# FIG2B
#############

e1<-ggplot(dat,aes(x=classN,y=Shannon, color=PExClust))+geom_boxplot(outlier.alpha=0.5, outlier.size=0.5, size=0.3)+geom_smooth(method = "lm",formula=y ~ x, se=TRUE,size=0.5, aes(group=PExClust))+xlab("Days to start of PEx therapy")+meg_stef_theme()+ggtitle("Diversity evolution")+theme(axis.text.x = element_text(angle = 45, vjust=0.5))+guides(color=guide_legend(override.aes=list(fill=NA)))+scale_color_paper(discrete = TRUE, palette = "CF",name="PEx types", breaks=c('PAT', 'AN1', 'AN2'))


###########
# CHAO1 TIME EVOLUTION
# FIG2C
############

e2<-ggplot(dat,aes(x=classN,y=Chao1, color=PExClust))+geom_boxplot(outlier.alpha=0.5, outlier.size=0.5, size=0.3)+geom_smooth(method = "lm",formula=y ~ x, se=TRUE,size=0.5, aes(group=PExClust))+xlab("Days to start of PEx therapy")+scale_y_continuous(limit = c(5, 70))+meg_stef_theme() +theme(axis.text.x = element_text(angle = 45, vjust=0.5))+ guides(color=guide_legend(override.aes=list(fill=NA)))+scale_color_paper(discrete = TRUE, palette = "CF",name="PEx types",breaks=c('PAT', 'AN1', 'AN2'))+ggtitle("Richness evolution")

rlegend<-get_legend(e2)

#############
# COMMUNITY TURNOVER/AITCHISON DISTANCE
# FIG2D
#############

##consecutive calculation by cluster
for (n in 1: 3){
if(n==1){
# cluster 1
me<-rownames(dat[dat$PExClust%in%"PAT",])
lab="PAT"
tit<-"Community turnover T"
yy<-"#F1BB7B"
ml<-24
mu<-36
grob1 <- grobTree(textGrob(expression(paste("T"["<24"],"=0.1650")), x=0.1,  y=0.9, hjust=0,
  gp=gpar(col="#7294D4", fontsize=9, fontface="italic")))
grob2 <- grobTree(textGrob(expression(paste("T"["24-60"],"=0.2586")), x=0.1,  y=0.78, hjust=0,
  gp=gpar(col="#02401B", fontsize=9, fontface="italic")))
}
if(n==2){
#cluster 2
me<-rownames(dat[dat$PExClust%in%"AN1",])
lab="AN1"
tit<-""
yy<-"#FD6467"
ml<-26
mu<-38
grob1 <- grobTree(textGrob(expression(paste("T"["<24"],"=0.4096")), x=0.1,  y=0.9, hjust=0,
  gp=gpar(col="#7294D4", fontsize=9, fontface="italic")))
grob2 <- grobTree(textGrob(expression(paste("T"["24-60"],"=0.3363")), x=0.1,  y=0.78, hjust=0,
  gp=gpar(col="#02401B", fontsize=9, fontface="italic")))
}
if(n==3){
#cluster 3
me<-rownames(dat[dat$PExClust%in%"AN2",])
lab="AN2"
tit<-""
yy<-"#5B1A18"
ml<-37
mu<-49
grob1 <- grobTree(textGrob(expression(paste("T"["<24"],"=0.3793")), x=0.1,  y=0.9, hjust=0,
  gp=gpar(col="#7294D4", fontsize=9, fontface="italic")))
grob2 <- grobTree(textGrob(expression(paste("T"["24-60"],"=0.2458")), x=0.1,  y=0.78, hjust=0,
  gp=gpar(col="#02401B", fontsize=9, fontface="italic")))
}

#generate data frame holding sampling time differences
E<-data.frame(sample=as.numeric(me))
rownames(E)<-me
dE<-dist(E)

#Aitchison distance among samples
#get rarified counts
ES<-otu_table(prune_samples(me,EC880))
#remove empty ASVs
ES<-ES[rowSums(ES)!=0,]
#add pseudocount for Aitchison
ES<-ES+runif(length(ES),0.0000001, 0.001)
#aitchison distance
dES<-aDist(t(ES))

#transform df to long format
dd<-sum(dE==1)
dfE<-data.frame(delta=rep(1,dd),AI=dES[dE==1])
for (i in 2:20){
dd<-sum(dE==i)
dfE<-rbind(dfE,(cbind(delta=rep(i,dd),AI=dES[dE==i])))
}

#Generate array with names of 2nd samples IDENTICAL to distance object
me<-labels(dE)
me<-me[-1]
nm<-me
while(length(me)>0){
me<-me[-1]
nm<-c(nm,me)
}

#df long format including SampleDay2T and 2nd sample
mett<-data.frame(dat)
days<-unname(mett[nm,"SampleDay2T"])
dd<-sum(dE==1)
dfET<-data.frame(delta=rep(1,dd),AI=dES[dE==1], SampleDays2T=days[dE==1], SecSample=nm[dE==1])
for (i in 2:20){
dd<-sum(dE==i)
dfET<-rbind(dfET,(cbind(delta=rep(i,dd),AI=dES[dE==i],SampleDays2T=days[dE==i], SecSample=nm[dE==i])))
}

#time groups <24 days and >23 days before treatment
tmp<-dfET
tmp$Time2Treatment<-(ifelse(tmp$SampleDays2T<24,"0-23 days","24-60 days"))
tmp$delta<-factor(tmp$delta, level=c("1","2","3","4","5","6","7","8", "9", "10","11","12","13","14","15","16","17","18", "19", "20"))
print(summary(aov(as.numeric(AI)~as.numeric(delta)+as.factor(Time2Treatment), data=tmp)))

######
##Turnover
tmp$title<-lab
#graph with mean value by sampling distance
p<-ggplot(tmp,aes(x=delta,y=as.numeric(AI), color=Time2Treatment))+stat_summary(geom="point",fun.y="mean", size=1)+geom_smooth(method = "lm",formula=y ~ x, se=TRUE,size=0.5, aes(group=Time2Treatment))+ylab("Aitchison distance")+xlab("Sampling distance (days)")+meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA)))+ggtitle(tit)+scale_color_paper(discrete = TRUE, palette = "Time", name="Days to\ntreatment", labels=c("<24", "24-60")) +scale_x_discrete(breaks=every_nth(n = 2))+ scale_y_continuous(breaks=seq(ml,mu,3))+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+annotation_custom(grob1)+annotation_custom(grob2)

if(n==1){t1<-p}
if(n==2){t2<-p}
if(n==3){t3<-p}

#quantify turnover (regression slopes)
sP<-tmp[tmp$Time2Treatment=="0-23 days",]
sP$AI<-as.numeric(sP$AI)
sP$delta<-as.numeric(sP$delta)
print("Turnover")
print(summary(lm(sP$AI~sP$delta)))

sB<-tmp[tmp$Time2Treatment=="24-60 days",]
sB$AI<-as.numeric(sB$AI)
sB$delta<-as.numeric(sB$delta)
print(summary(lm(sB$AI~sB$delta)))
}

alegend<-get_legend(t1)


#######
#ASSEMBLY FIGURE2
######
#arrange plots
legends1<-plot_grid(rlegend,alegend, labels=c("",""), align = "v", axis="b",ncol=1,label_size = 12, label_fontface = "bold")
legends<-plot_grid(legends1,dlegend, labels=c("",""), align = "h", axis="t",ncol=2,label_size = 12, label_fontface = "bold",rel_widths=c(0.5,0.5))
up<-plot_grid(d1+theme(legend.position = "none"), d2+theme(legend.position = "none"),d3+theme(legend.position = "none"), labels=c("A","",""), align = "h", ncol=3,label_size = 12, label_fontface = "bold")
mid<-plot_grid(e1+theme(legend.position = "none"), e2+theme(legend.position = "none"), labels=c("B","C"), align = "h", ncol=2,label_size = 12, label_fontface = "bold")
low<-plot_grid(t1+theme(legend.position = "none"), t2+theme(legend.position = "none"),t3+theme(legend.position = "none"), labels=c("D","",""), align = "h", ncol=3,label_size = 12, label_fontface = "bold")
mid<-plot_grid(mid, legends,ncol=2,labels=c("","",""), rel_widths=c(0.66,0.33))
figure2<-plot_grid(up,mid,low, labels=c("","",""), align = "hv",axis="r", ncol=1,label_size = 12, label_fontface = "bold")

#save plots
save(figure2, file="figure2.rda")
save_plot("figure2.pdf",figure2, ncol=1, base_height=7,  base_asp=1.1)
save_plot("figure2.png",figure2, ncol=1, base_height=7,  base_asp=1.1)


###########
##ANALYSES
###########

#############
#SHANNON TIME EVOLUTION
#############

m1<-(lmerTest::lmer(Shannon~as.factor(D23)+(1|Agegroup) +(1|Patient)+(1|ECid), data=dat[dat$PExClust=="PAT",]))
m2<-(lmerTest::lmer(Shannon~as.numeric(classN)+(1|Agegroup) +(1|Patient)+(1|ECid)+(1|Agegroup), data=dat[dat$PExClust=="AN1",]))
m3<-(lmerTest::lmer(Shannon~as.factor(D23)+(1|Agegroup) +(1|Patient)+(1|ECid), data=dat[dat$PExClust=="AN2",]))

#calulate effect sizes
am1<-anova(m1)
am1$F
am1$Num
am1$Den
F_to_eta2(am1$F,am1$Num,am1$Den)

am1<-anova(m2)
am1$F
am1$Num
am1$Den
F_to_eta2(am1$F,am1$Num,am1$Den)

am1<-anova(m3)
am1$F
am1$Num
am1$Den
F_to_eta2(am1$F,am1$Num,am1$Den)

#Wilcoxon
wilcox.test(Shannon~D23, data=dat[dat$PExClust=="PAT",])
wilcox.test(Shannon~D23, data=dat[dat$PExClust=="AN1",])
wilcox.test(Shannon~D23, data=dat[dat$PExClust=="AN2",])

###########
#CHAO1 TIME EVOLUTION
############

#Sign CF, F, not N
m1<-(lmerTest::lmer(Chao1~as.factor(D23)+(1|Agegroup) +(1|Patient) +(1|ECid), data=dat[dat$PExClust=="PAT",]))
m2<-(lmerTest::lmer(Chao1~as.factor(D23) +(1|Agegroup)+(1|Patient)+(1|ECid), data=dat[dat$PExClust=="AN1",]))
m3<-(lmerTest::lmer(Chao1~as.factor(D23)+(1|Agegroup) +(1|Patient)+(1|ECid), data=dat[dat$PExClust=="AN2",]))

#models, effect sizes
am1<-anova(m1)
am1
am1$F
am1$Num
am1$Den
F_to_eta2(am1$F,am1$Num,am1$Den)

am1<-anova(m2)
am1
am1$F
am1$Num
am1$Den
F_to_eta2(am1$F,am1$Num,am1$Den)

am1<-anova(m3)
am1
am1$F
am1$Num
am1$Den
F_to_eta2(am1$F,am1$Num,am1$Den)

#Wilcoxon
wilcox.test(Chao1~D23, data=dat[dat$PExClust=="PAT",])
wilcox.test(Chao1~D23, data=dat[dat$PExClust=="AN1",])
wilcox.test(Chao1~D23, data=dat[dat$PExClust=="AN2",])


##################
###SUPP FIGURES
##################

#########
#SUPP3 BC
#########
#arrange plots
sSubj<-plot_grid(sp1+theme(legend.position ="none"),sp2+theme(legend.position ="none"),sp3+theme(legend.position ="none"),sl1, labels=c("C","","",""), align = "hv", ncol=4,label_size = 12, label_fontface = "bold", axis="tblr")
sEC<-plot_grid(sec1+theme(legend.position ="none"),sec2+theme(legend.position ="none"),sec3+theme(legend.position ="none"),sl2, labels=c("D","","",""), align = "hv", ncol=4,label_size = 12, label_fontface = "bold",axis="tblr")
s2CD<-plot_grid(sSubj,sEC,ncol=1)

#save plots
save(s2CD, file="SI_figure3CD.rda")
save_plot("SI_figure3CD.pdf",s2CD, ncol=1, base_height=5)
save_plot("SI_figure3CD.png",s2CD, ncol=1, base_height=5)


#######
#SUPP5
#######
#arrange plots
sl<-get_legend(sw1)
ss<-plot_grid(sw1+theme(legend.position ="none"),sw2+theme(legend.position ="none"),sw3+theme(legend.position ="none"),sl, labels=c("A","","",""), align = "vh", ncol=4,label_size = 12, label_fontface = "bold",axis="tblr")
rr<-plot_grid(sr1+theme(legend.position ="none"),sr2+theme(legend.position ="none"),sr3+theme(legend.position ="none"),sl, labels=c("B","","",""), align = "vh", ncol=4,label_size = 12, label_fontface = "bold",axis="tblr")
s3<-plot_grid(ss,rr,ncol=1)

#save plots
save(s3, file="SI_figure5.rda")
save_plot("SI_figure5.pdf",s3, ncol=1, base_height=5)
save_plot("SI_figure5.png",s3, ncol=1, base_height=5)

