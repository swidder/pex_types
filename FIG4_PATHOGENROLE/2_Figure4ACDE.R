#Author Stefanie Widder 24062023

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(phyloseq)
library(igraph)
library(plyr)
source("../HARMONIZE_FIGURES/1_plot_utils.R")

#########
##AUX
#########

############
# page rank on giant component
############
get.pageRankGC<-function(g){
#giant
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#page rank
b<-page_rank(tmp, directed=F)$vector
#join to vertex array 
d<-setNames(array(rep(NA,vcount(g))),V(g)$name)
d[names(b)]<-b
d
}

############
# extract NW_IDs whith majority in clusterK and return only NWS with <10 iSamples
###########
get_NWs_majority<-function(maP,spl){
#find majority k
tmp<-lapply(maP, function(x,spl){table(spl[x])},spl=spl)
tmp<-lapply(tmp, function(x){x[order(-x)]})
tmp<-lapply(tmp, function(x){names(x)[1]})

#find representative (dummy) samples
rep<-array()
for (j in 1:length(tmp)){
rep[j]<-names(which(spl[maP[[j]]]==tmp[[j]]))[1]
}
names(rep)<-names(tmp)
rep
}


#####
# make GC and remove pathogen, return subsetted graph
#####
remove_pathogensGC<-function(g, TAX, who){
#subset giant component
co <- components(g, mode = 'STRONG')
g<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#find node indices to remove
#if all CF pathogens together
if(who=="S4"){
me<-V(g)$name[TAX[V(g)$name,"S4"]=="C"]
#na from non-pathogens patient-wise composition
me<-na.omit(me)
}else{
#single ASV
me<-who
if(!me%in%V(g)$name){return(NA)}
}
g - me
}


########
#Change of components number (strong components; normalized to GC sizes)
#########
comp_noComp<-function(g,g1){
#subset giant componet of graph g
co <- components(g, mode = 'STRONG')
g<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#graph g1 is already subsetted to giant component

#breaking apart +/-pathogen, normalized to original size 
cc1<-clusters(g1)$no/vcount(g1)
cc2<-clusters(g)$no/vcount(g)
#delta
cc1-cc2
}

########
#Size change of giant component (GC; normalized to original sizes)
########
comp_bigClustSize<-function(g,g1){
#subset giant componet of graph g
co <- components(g, mode = 'STRONG')
g<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#graph g1 is already subsetted to giant component

#biggest cluster +/- pathogen, normalized to original size
cc1<-sort(clusters(g1)$csize, decreasing=T)[1]/vcount(g1)
cc2<-sort(clusters(g)$csize, decreasing=T)[1]/vcount(g)
#Delta
cc1-cc2
}

#######
#Change of walktrap modularity on topology only (unsigned adjacencies, GC)
########
comp_walktrapMod<-function(g,g1){
#subset giant componet of graph g
co <- components(g, mode = 'STRONG')
g<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#unsign weights
E(g)$weight<-abs(E(g)$weight)
#graph g1 is already subsetted to giant component
#unsign weights
E(g1)$weight<-abs(E(g1)$weight)
#Delta walktrap modularity
cc1<-modularity(g1,(walktrap.community(g1))$membership)
cc2<-modularity(g,(walktrap.community(g))$membership)
cc<-cc1-cc2
#fetch error modularity cases
if(if.na(cc1)|if.na(cc2)){cc<-NA}
cc
}

#######
#Change of walktrap modularity on topology only (unsigned adjacencies, GC, absolute modularity values)
########
comp_ABSwalktrapMod<-function(g,g1){
#subset giant componet of graph g
co <- components(g, mode = 'STRONG')
g<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#unsign weights
E(g)$weight<-abs(E(g)$weight)
#graph g1 is already subsetted to giant component
#unsigned weights
E(g1)$weight<-abs(E(g1)$weight)
#Delta walktrap modularity for absolute change
cc1<-modularity(g1,(walktrap.community(g1))$membership)
cc2<-modularity(g,(walktrap.community(g))$membership)
cc<-abs(cc1-cc2)
if(if.na(cc1)|if.na(cc2)){cc<-NA}
cc
}


##########
# Identify dominant pathogen in NW, trajectory wise
########
get_pathogen<-function(maP, OTU, Plist){
pathogen<-array()
for (i in 1: length(maP)){
#subset to available samples in NW
me<-maP[[i]][maP[[i]]!="iSample"]
tmp<-OTU[,colnames(OTU)%in%me]
#subset to pathogens
tmp<-tmp[rownames(tmp)%in%Plist,]

#remember pathogen with most counts
pathogen[i]<-names(which.max(rowSums(tmp)))
#in case no pathogen is available set NA
if(pathogen[i]==0){pathogen[i]<-NA}
}
#array of pathogens by NW
pathogen
}



########
#MAIN
########

load("../DATA_INPUT/EC880F1B.rda")
dat<-as.data.frame(sample_data(EC880))
#subset to NW inclusion criteria
dat<-dat[!is.na(dat$PExClust),]
#taxonomy file
TAX<-as.data.frame(tax_table(EC880))
OTU<-as.data.frame(otu_table(EC880))
Plist<- rownames(TAX)[TAX$S4=="C"]


#################
#FLAT NW FORMAT
###

#find all NW directories
me<-dir("../DATA_INPUT/NWs", pattern="EC", full.names=T)
#rememeber cluster category by sample
spl<-dat$K3
names(spl)<-rownames(dat)

#intitialize dfs
ASV<-setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("patient","EC", "K3","NW_ID","D2T","ASV","S4","pageRank"))
ABX<-setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("patient","EC", "NW_IDs", "sampleID1", "K3","D2T","removedASV","absModularityGCP","modularityGCP","componentsGCP","sizeGCP"))

#for each time series
for (i in 1: length(me)){
d<-strsplit(me[i],"/")[[1]][4]
dd<-unlist(strsplit(d, "_"))
f<-paste("NWs_",dd[1], "_patient",dd[3],".rda", sep="")
targ<-paste(me[i], "/",f, sep="")

#read NW
NW<-get(load(targ))

#remember samples used for NW inference
maP<-NW$samples

#NW classification to PExClust by sample majority
mes<-get_NWs_majority(maP,spl)

#assess disruptivity 
nm<-names(mes)
k3<-dat[unname(mes),"K3"]

#intitialize
c_p<-array()
s_p<-array()
m_p<-array()
am_p<-array()
rem_pat<-array()
pathogen<-get_pathogen(maP, OTU, Plist)

#add ASV features in NW to ASV
#for each NW
for (h in 1: length(mes)){
#get NW
g<-NW$fullNW[[h]]

#ASV importance (page rank)
ASV<-rbind(ASV,data.frame(patient=rep(NW$patient, vcount(g)),EC=rep(dd[1], vcount(g)),K3=rep(unname(unlist(k3[h])), vcount(g)),NW_ID=rep(nm[h], vcount(g)),D2T=rep(as.numeric(NW$days2treat[h]),vcount(g)),ASV=V(g)$name,S4=TAX[V(g)$name,"S4"],pageRank=get.pageRankGC(g)))


#subset to NW with dominant pathogen in giant component
#remove dominant pathogen from giant component to assess disruptivity
rem_pat[h]<-pathogen[h]
c_p[h]<-s_p[h]<-m_p[h]<-am_p[h]<-NA
#now overwrite
if(!is.na(rem_pat[h])&(pathogen[h]%in%V(g)$name)){
g1<-remove_pathogensGC(g,TAX, rem_pat[h])
if(!is.na(g1)){
c_p[h]<-comp_noComp(g,g1)
s_p[h]<-comp_bigClustSize(g,g1)
m_p[h]<-comp_walktrapMod(g,g1)
am_p[h]<-comp_ABSwalktrapMod(g,g1)}
}}

#remember all features, flat file format
ABX<-rbind(ABX,data.frame(patient=rep(NW$patient, length(mes)),EC=rep(dd[1], length(mes)) ,NW_IDs=nm,sampleID1=mes,K3=k3, D2T=as.numeric(NW$days2treat[nm]),removedASV=rem_pat,absModularityGCP=am_p, modularityGCP=m_p,componentsGCP=c_p, sizeGCP=s_p))

}


##########
#PLOTTING
#########

#######################
## TOPOLOGY of AEROBICITY GROUPS Fig4A
###############

#subset core ASVs 
orgs<-unique(ASV$ASV)
col<-rownames(TAX)
morgs<-setdiff( orgs, col)
ASV<-ASV[!ASV$ASV%in%morgs,]
#remove uncultivated organisms
tmp<-ASV[ASV$S4%in%c("N","C","F","A"),]

#set plotting order
tmp$S4<-factor(tmp$S4,levels=c("C","A","F","N"))
tmp$K3<-factor(tmp$K3,levels=c("PAT","AN1", "AN2"))
tmp$patient<-gsub("ient","ient_",tmp$patient)
comp=list(c("A","C"),c("F","C"),c("C","N") )

#fix color code
coh<-paper_pal("Paper")(11)
names(coh)<-levels(as.factor((as(sample_data(EC880),"data.frame"))$Patient))
coh1<-coh[!names(coh)%in%"H"]


#split into facettes
for (n in 1: 3){
if(n==1){
# cluster 1
me<-tmp[tmp$K3%in%"PAT",]
lab="PAT"
tit<-"Community importance of taxonomic groups"
yy<-"#F1BB7B"
tl<-c(0.3,0.37,0.44)

}
if(n==2){
#cluster 2
me<-tmp[tmp$K3%in%"AN1",]
lab="AN1"
tit<-""
yy<-"#FD6467"
tl<-c(0.25,0.3,0.35)

}
if(n==3){
#cluster 3
me<-tmp[tmp$K3%in%"AN2",]
lab="AN2"
tit<-""
yy<-"#5B1A18"
tl<-c(0.1,0.13,0.16)

}

me$title<-lab
p<-ggplot(me,aes(x=as.factor(S4),y=pageRank))+geom_jitter(aes(col=patient),shape=16, position=position_jitter(0.2), size=0.5, alpha=0.5) +scale_color_manual(values=coh1,name="Cohort",labels=c("Subject A","Subject B","Subject C","Subject D","Subject E", "Subject F","Subject G","Subject I","Subject J", "Subject K"))+geom_boxplot(outlier.shape=NA, fill=NA)+stat_compare_means(comparisons=comp,label = "p.signif", size=4.2, method="wilcox.test",label.y=tl)+ylab("Page rank")+xlab("Aerobicity class") + meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA,size=1.5,alpha=1)))+ ggtitle(tit)+facet_grid(. ~ title)+theme(strip.background =element_rect(fill=alpha(yy, 0.5)),strip.text = element_text(colour = 'black', size=10))+theme(legend.key.size = unit(0.5, "cm"))

ifelse(n==1,d4A<-p,ifelse(n==2,d4B<-p,d4C<-p))
}

me<-tmp
l1<-get_legend(ggplot(me,aes(x=as.factor(S4),y=pageRank))+geom_jitter(aes(col=patient),shape=16, position=position_jitter(0.2), size=0.5, alpha=0.5) +scale_color_manual(values=coh1,name="Cohort",labels=c("Subject A","Subject B","Subject C","Subject D","Subject E", "Subject F","Subject G","Subject I","Subject J", "Subject K")) + meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA,size=1.5,alpha=1)))+ theme(legend.key.size = unit(0.95, 'lines')) )



#######
#TREATMENT SIMULATION (giant component)
#######

#check NW numbers with pathogen in giant component
sum(!is.na(tmp$sizeGCP))
#[1] 312


tmp<-ABX
tmp$absSizeGCP<-abs(tmp$sizeGCP)
#fix color code
cfcol<-c("PAT"="#F1BB7B","AN1"="#FD6467","AN2"="#5B1A18")

#Density plots
#https://stackoverflow.com/questions/51385455/geom-density-y-axis-goes-above-1
#https://towardsdatascience.com/how-to-compare-two-or-more-distributions-9b06ee4d30bf

###########
#scaled density plots
###########

p3A<-ggplot(tmp, aes(x=sizeGCP,..scaled..,fill=K3)) + geom_density(alpha=.5)+ meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA)))+ scale_fill_manual(values=cfcol,name="PEx types",breaks=c( 'PAT','AN1', 'AN2'))+ggtitle("Community reduction\nupon treatment")+xlab(expression(Delta*"Nodes (%)"))+ylab("Probability density")

p3B<-ggplot(tmp, aes(x=absModularityGCP,..scaled..,fill=K3)) + geom_density(alpha=.5)+ meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA)))+ scale_fill_manual(values=cfcol,name="PEx types",breaks=c( 'PAT','AN1', 'AN2'))+ggtitle("Modularity change\nupon treatment")+xlab(expression(Delta*"Modularity (%)"))+ylab("Probability density")

p3C<-ggplot(tmp, aes(x=componentsGCP,..scaled..,fill=K3)) + geom_density(alpha=.5)+ meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA)))+ scale_fill_manual(values=cfcol,name="PEx types",breaks=c( 'PAT','AN1', 'AN2'))+ggtitle("Community breakage\nupon treatment")+xlab(expression(Delta*"Components (%)"))+ylab("Probability density")





##############
#FIGURE4
############



f41<-plot_grid(d4A+theme(legend.position = "none"),d4B+theme(legend.position = "none"),d4C+theme(legend.position = "none"), ncol=3,align = "hv",labels=c("A","",""),label_size = 12, label_fontface = "bold")

l3<-get_legend(p3C)
f43<-plot_grid(p3A+theme(legend.position = "none"),p3B+theme(legend.position = "none"),p3C+theme(legend.position = "none"), ncol=3,align = "hv",labels=c("C","D","E"),label_size = 12, label_fontface = "bold")

load("p4B.rda")
l2<-get_legend(p4B)
legs<-plot_grid(l1,l2,l3,ncol=1)

plts<-plot_grid(f41,p4B+theme(legend.position = "none"),f43, ncol=1, labels=c("","B",""),label_size = 12, label_fontface = "bold")

F4<-plot_grid(plts,legs,ncol=2, rel_widths=c(1,0.2))
save(F4,file="figure4.rda")

save_plot("figure4.pdf",F4, base_height=7, base_asp=1)
save_plot("figure4.png",F4, base_height=7, base_asp=1)


############
#ANALYSIS density distribution
############

#########
#kolmogorov smirnov test for similarity of the cumulative distributions
###########

#two sided KS absolute modularity change
#CF vs N
x<-tmp[tmp$K3=="PAT","absModularityGCP"]
y<-tmp[tmp$K3=="AN2","absModularityGCP"]
ks.test(x,y)
#D = 0.44777, p-value = 5.368e-10
#CF vs F
y<-tmp[tmp$K3=="AN1","absModularityGCP"]
ks.test(x,y)
#D = 0.29847, p-value = 0.002745
#N vs F
x<-tmp[tmp$K3=="AN2","absModularityGCP"]
ks.test(x,y)
#D = 0.4277, p-value = 1.677e-06


###########
#two sided KS size reduction of GC
#CF vs N
x<-tmp[tmp$K3=="PAT","sizeGCP"]
y<-tmp[tmp$K3=="AN2","sizeGCP"]
ks.test(x,y)
#D = 0.3772, p-value = 3.831e-08
#CF vs F
y<-tmp[tmp$K3=="AN1","sizeGCP"]
ks.test(x,y)
#D = 0.33273, p-value = 0.0002699
#N vs F
x<-tmp[tmp$K3=="AN2","sizeGCP"]
ks.test(x,y)
#D = 0.27867, p-value = 0.007089


###########
#two sided KS disruptivity GC
#CF vs N
x<-tmp[tmp$K3=="PAT","componentsGCP"]
y<-tmp[tmp$K3=="AN2","componentsGCP"]
ks.test(x,y)
#D = 0.49458, p-value = 1.077e-13
#CF vs F
y<-tmp[tmp$K3=="AN1","componentsGCP"]
ks.test(x,y)
#D = 0.56424, p-value = 1.487e-11
#N vs F
x<-tmp[tmp$K3=="AN2","componentsGCP"]
ks.test(x,y)
#D = 0.26049, p-value = 0.01445
































