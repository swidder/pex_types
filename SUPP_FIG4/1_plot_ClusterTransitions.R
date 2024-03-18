library(phyloseq)
library(tidyverse)
library(igraph)
library(magick)
library(MESS)
library(cowplot)
library(qgraph)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggsankey)
load("../DATA_INPUT/EC880F1B.rda")
source("../HARMONIZE_FIGURES/1_plot_utils.R")

######AUX
get_transition<-function(tmp){
#weighted adj matrix per patient
m<-matrix(0,nrow=3, ncol=3)
colnames(m)<-c("PAT","AN1", "AN2")
rownames(m)<-c("PAT","AN1", "AN2")
#by time series
nEC<-names(table(tmp$ECid))[table(tmp$ECid)>0]
for (j in 1: length(nEC)){
ttmp<-tmp[tmp$ECid==nEC[j],]
#run window and count
for (k in 1: (nrow(ttmp)-1)){
m[ttmp$K3[k],ttmp$K3[k+1]]<-m[ttmp$K3[k],ttmp$K3[k+1]]+1
}}
#weighted transition matrix row->col
m
}
#############


#####MAIN

#subset PExClust, patient, Time, sample
met<-as(sample_data(EC880),"data.frame")
s<-data.frame(met[,c(2:3,8,20:21)])
s$SampleID<-rownames(met)
s<-s[!is.na(s$PExClust),]

#ordertime series by patients
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
ss$ECid<-factor(ss$ECid, levels=c("1","2","3","4","5","6","7","8","10","11","13","14","15","16","17","18"))
ss$K3<-factor(ss$K3,levels=c("PAT", "AN1", "AN2"))

#Pex types by patient by time lines
# add hline after every patient
p1<-ggplot(ss, aes(-SampleDay2T, ECid, color=factor(K3))) + geom_point() + geom_hline(aes(yintercept=3.5))+ geom_hline(aes(yintercept=4.5))+ geom_hline(aes(yintercept=5.5))+ geom_hline(aes(yintercept=6.5))+ geom_hline(aes(yintercept=8.5))+ geom_hline(aes(yintercept=9.5))+ geom_hline(aes(yintercept=10.5))+ geom_hline(aes(yintercept=13.5))+ geom_hline(aes(yintercept=14.5))+ggtitle("Time series by PEx type")+scale_colour_manual(values=c("#F1BB7B", "#FD6467", "#5B1A18"),name="PEx type")+ meg_stef_theme()+ylab("PEx time series")+xlab("Days to start of PEx therapy")+scale_x_continuous(breaks=seq(-60,0, by=10))+theme(legend.position="bottom")


#########TRANSITIONS
#group by patient
ss <- s %>% group_by(Patient)
ss$ECid<-factor(ss$ECid, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))

#make a neighbour list with weights (transition events) by patient
pl<-list()
pn<-names(table(ss$Patient))
for (i in 1: length(table(ss$Patient))){
tmp<- ss[ss$Patient==pn[i],]
tmp <- tmp %>% group_by(ECid)
pl[[i]]<-get_transition(tmp)
}
names(pl)<-names(table(ss$Patient))

#make graphs
# transitions in %
pp<-lapply(pl, MESS::round_percent)
nm<-toupper(letters[c(1:7,9:11)])

#save individual transition graphs as pdf
for (i in 1: length(nm)){
f<-paste("P_",nm[i],"_g",i,".pdf", sep="")
pdf(f)
qgraph(pp[[i]],edge.labels=TRUE, trans=FALSE, fade=FALSE, label.prop = 0.9, vsize=20, edge.label.cex=3, repulsion=1, mar=c(6,9.5,9,9.5),layout="groups")
dev.off()
}

#read in to arrange with cowplot
ppA <- ggdraw() + draw_image("P_A_g1.pdf", scale=0.95) + theme(plot.background = element_rect(color = "black"))
ppB <- ggdraw() + draw_image("P_B_g2.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppC <- ggdraw() + draw_image("P_C_g3.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppD <- ggdraw() + draw_image("P_D_g4.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppE <- ggdraw() + draw_image("P_E_g5.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppF <- ggdraw() + draw_image("P_F_g6.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppG <- ggdraw() + draw_image("P_G_g7.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppI <- ggdraw() + draw_image("P_I_g8.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppJ <- ggdraw() + draw_image("P_J_g9.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))
ppK <- ggdraw() + draw_image("P_K_g10.pdf", scale=0.95)+ theme(plot.background = element_rect(color = "black"))

SIXBperc<-plot_grid(ppA,ppB,ppC,ppD,ppE,ppF,ppG,ppI,ppJ,ppK,labels=c("Sub A", "Sub B", "Sub C", "Sub D","Sub E", "Sub F", "Sub G", "Sub I", "Sub J", "Sub K"), ncol=2, nrow=5,label_size = 10, label_fontface = "plain")

#fix title margin
SIXBPT<-annotate_figure(SIXBperc, top=arrangeGrob(textGrob("PEx type transitions by subject (%)", hjust=0.75, gp=gpar(fontsize=10)
), zeroGrob(),widths = unit(1, 'npc'),  heights = unit(c(0.45, 0.2), c('cm', 'npc')),  as.table = FALSE      )      )

#SANKEY from total transition data
#sum transitions (EC-wise)
tot<-pl[[1]]
for (i in 2: length(pl)){
tot<-tot+pl[[i]]
}
df<-data.frame()
for (i in 1: nrow(tot)){
for (j in 1: ncol(tot)){
df<-rbind(df,t(replicate(tot[i,j],c(rownames(tot)[i],colnames(tot)[j]))))
}
}
colnames(df)<-c("FROM", "TO")

#formate
dat<-df%>%make_long(FROM, TO)
#order PEx types
dat$node<-factor(dat$node,levels=c("PAT", "AN1", "AN2"))
#add frequency
reagg <- dat%>% dplyr::group_by(node)%>% tally()
dat2 <- merge(dat,reagg,  by.x = "node", by.y = "node", all.x = TRUE)

#calulate % transitions
dat3<-dat%>%count(x, node, next_node)%>%drop_na
g<-c("PAT","AN1","AN2")
perc<-array()
for (i in 1:length(g)){
tmp<-dat3[dat3$node==g[i],]
tmp$next_node
perc<-c(perc,MESS::round_percent(tmp$n))
}
perc<-perc[-1]
dat3$P<-perc

#Sankey plot
SIS<-ggplot(dat2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node,  label = node)) + scale_fill_manual(values=c("#F1BB7B", "#FD6467", "#5B1A18"), name="PEx types") + meg_stef_theme()  + geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE)  + xlab("Transitions") + geom_sankey_label(size = 3, color = "black", fill = "white") + theme(axis.title = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank())+ggtitle("Transitions among PEx types")

#arrange
T<-plot_grid(NULL,SIS, labels=c("A",""),ncol=2, rel_widths=c(0.05,1))
SIA<-plot_grid(T, p1, labels=c("","B"),nrow=2)
B<-plot_grid(NULL,SIXBPT, labels=c("C",""),ncol=2, rel_widths=c(0.1,1))
SIT<-plot_grid(SIA, B, ncol=2, axis="t", align="h")

#save figure
save_plot(SIT, file="SI_figure4.pdf", base_height=8, base_asp=0.95)
save_plot(SIT, file="SI_figure4.png", base_height=8, base_asp=0.95)













