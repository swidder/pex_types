#Author Stefanie Widder 24062023

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(igraph)
library(plyr)
library(dplyr)
source("../HARMONIZE_FIGURES/1_plot_utils.R")
library(effectsize)
library(lmerTest)

#########
##AUX
#########


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


############
#Betweenness centrality giant component
############
get.betwGC<-function(g){
#giant
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#centrality
b<-betweenness(tmp,directed=F, normalize=T, weight=abs(E(tmp)$weight))

#join to vertex array 
d<-setNames(array(rep(NA,vcount(g))),V(g)$name)
d[names(b)]<-b
d
}

#########
#Vertex count on giant component
#########
vcountGC<-function(g){
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#ecc
vcount(tmp)
}

##########
# edgecount giant component
##########
ecountGC<-function(g){
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#ecc
ecount(tmp)
}

##############
#clustering coefficient giant component 
##############
clustCoeffGC<-function(g){
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
transitivity(tmp)
}


######### Mean/Median diversity from sample set building graph; iSamples removed
get.diversity<-function(smp,PO, div){
smp<-smp[!smp=="iSample"]
tmp<-estimate_richness(prune_samples(smp,PO), split=T, measures=div)
c(mean(tmp[[1]]), median(tmp[[1]]))
}




########
#MAIN
########

####
# LOAD DATA
####

load("../DATA_INPUT/EC880F1B.rda")
dat<-as.data.frame(sample_data(EC880))
dat<-dat[!is.na(dat$PExClust),]
#taxonomy file
TAX<-as.data.frame(tax_table(EC880))


####
#FLATTEN NW INFORMATION
#####


#find all NW directories
me<-dir("../DATA_INPUT/NWs", pattern="EC", full.names=T)
#rememeber cluster category by sample
spl<-dat$K3
names(spl)<-rownames(dat)

#intitialize dfs
TOP<-setNames(data.frame(matrix(ncol = 12, nrow = 0)), c("patient","EC", "NW_IDs", "sampleID1", "K3","D2T","ecount", "vcount","betweennessC", "Clust", "Shannon", "Chao1"))
#for each time series
for (i in 1: length(me)){
#construct path to file
d<-strsplit(me[i],"/")[[1]][3]
dd<-unlist(strsplit(d, "_"))
f<-paste("NWs_",dd[1], "_patient",dd[3],".rda", sep="")
targ<-paste(me[i], "/",f, sep="")

#read NW
NW<-get(load(targ))

#remember samples used for NW inference
maP<-NW$samples

#NW classification to PExClust by sample majority
mes<-get_NWs_majority(maP,spl)

#collect NW features
nm<-names(mes)
k3<-unname(unlist((dat[unname(mes),"K3"])))

#initialize arrays
edges<-array()
nodes<-array()
#Graph-wise clustering, Shannon mean, Chao1 mean
clc<-array()
#graph wise
S<-array()
R<-array()

#for every NW
for (h in 1: length(mes)){
g<-NW$fullNW[[h]]
edges[h]<-ecountGC(g)
nodes[h]<-vcountGC(g)
#graph-wise
clc[h]<-clustCoeffGC(g)
tmp<-get.diversity(maP[[names(mes[h])]],EC880, "Shannon" )
S[h]<-tmp[1]
tmp<-get.diversity(maP[[names(mes[h])]],EC880, "Chao1" )
R[h]<-tmp[1]
}

names(edges)<-names(mes)
names(nodes)<-names(mes)

#########
# Topological NW features
TOP<-rbind(TOP, data.frame(patient=rep(NW$patient, length(mes)),EC=rep(dd[1], length(mes)) ,NW_IDs=nm,sampleID1=mes,K3=k3, D2T=as.numeric(NW$days2treat[nm]),ecount=edges, vcount=nodes,betweennessC=as.numeric(NW$betweennessC[nm]), Clust=clc,Shannon=S, Chao1=R ))

}

###########
#RESULT LOOKUP
############
#GLMM K3
############

#BetweennessCentrality
m1<-(lmerTest::lmer(betweennessC~as.factor(K3)+Shannon+Chao1+(1|patient), data=TOP))
#get effect size
am1<-anova(m1)
am1$F
am1$Num
am1$Den
es1<-F_to_eta2(am1$F,am1$Num,am1$Den)
Eta2 (partial) |       95% CI
-----------------------------
0.35           | [0.10, 1.00]
0.12           | [0.00, 1.00]
0.05           | [0.00, 1.00]
#K3 **

#Clustering
m2<-(lmerTest::lmer(Clust~as.factor(K3)+Shannon+Chao1+(1|patient), data=TOP))
#get effect size
am2<-anova(m2)
am2$F
am2$Num
am2$Den
es2<-F_to_eta2(am2$F,am2$Num,am2$Den)

Eta2 (partial) |       95% CI
-----------------------------
0.11           | [0.00, 1.00]
0.02           | [0.00, 1.00]
8.88e-03       | [0.00, 1.00]
# K3*

#Size giant component
m3<-(lmerTest::lmer(vcount~as.factor(K3)+Shannon+Chao1+(1|patient), data=TOP))
#get effect size
am3<-anova(m3)
am3$F
am3$Num
am3$Den
es3<-F_to_eta2(am3$F,am3$Num,am3$Den)
Eta2 (partial) |       95% CI
-----------------------------
0.05           | [0.03, 1.00]
5.57e-03       | [0.00, 1.00]
0.02           | [0.01, 1.00]

#K3, Chao1 ***

#Connectivity
m4<-(lmerTest::lmer(ecount~as.factor(K3)+Shannon+Chao1+(1|patient), data=TOP))
#get effect size
am4<-anova(m4)
am4$F
am4$Num
am4$Den
es4<-F_to_eta2(am4$F,am4$Num,am4$Den)
Eta2 (partial) |       95% CI
-----------------------------
0.58           | [0.43, 1.00]
0.17           | [0.06, 1.00]
0.26           | [0.19, 1.00]
#K3, Shannon, Chao ***


############
#PLOTTING EFFECT SIZES AND P VALUES
#######
df<-data.frame(omega2_1= es1$Eta2_partial,omega2_0=0, p=am1[6][[1]])
#same covariate labels
df$covariates<-c("PEx types","Shannon","Chao1")
rownames(df)<-c("PEx types","Shannon","Chao1")
df$Pvalue<-ifelse(df$p>=0.05,"n.s.",ifelse(df$p>=0.01, "*", ifelse(df$p>=0.001, "**", "***")))
p1<-df%>%ggplot(aes( x = omega2_1, y = reorder(covariates,omega2_1),fill=Pvalue)) + geom_bar(stat="identity", position="dodge")+ meg_stef_theme()+xlab(expression(paste("Partial ",eta^2)))+ylab("Covariates")+ggtitle("Effect sizes on betweenness centrality")+scale_fill_paper(discrete = TRUE, palette = "Bi2")

df<-data.frame(omega2_1= es2$Eta2_partial,omega2_0=0, p=am2[6][[1]])
#same covariate labels
df$covariates<-c("PEx types","Shannon","Chao1")
rownames(df)<-c("PEx types","Shannon","Chao1")
df$Pvalue<-ifelse(df$p>=0.05,"n.s.",ifelse(df$p>=0.01, "*", ifelse(df$p>=0.001, "**", "***")))
p2<-df%>%ggplot(aes( x = omega2_1, y = reorder(covariates,omega2_1),fill=Pvalue)) + geom_bar(stat="identity", position="dodge")+ meg_stef_theme()+xlab(expression(paste("Partial ",eta^2)))+ylab("Covariates")+ggtitle("Effect sizes on clustering coefficient")+scale_fill_paper(discrete = TRUE, palette = "Bi2")


df<-data.frame(omega2_1= es3$Eta2_partial,omega2_0=0, p=am3[6][[1]])
#same covariate labels
df$covariates<-c("PEx types","Shannon","Chao1")
rownames(df)<-c("PEx types","Shannon","Chao1")
df$Pvalue<-ifelse(df$p>=0.05,"n.s.",ifelse(df$p>=0.01, "*", ifelse(df$p>=0.001, "**", "***")))
p3<-df%>%ggplot(aes( x = omega2_1, y = reorder(covariates,omega2_1),fill=Pvalue)) + geom_bar(stat="identity", position="dodge")+ meg_stef_theme()+xlab(expression(paste("Partial ",eta^2)))+ylab("Covariates")+ggtitle("Effect sizes on vertex count (size)")+scale_fill_paper(discrete = TRUE, palette = "Bi2")

df<-data.frame(omega2_1= es4$Eta2_partial,omega2_0=0, p=am4[6][[1]])
#same covariate labels
df$covariates<-c("PEx types","Shannon","Chao1")
rownames(df)<-c("PEx types","Shannon","Chao1")
df$Pvalue<-ifelse(df$p>=0.05,"n.s.",ifelse(df$p>=0.01, "*", ifelse(df$p>=0.001, "**", "***")))
p4<-df%>%ggplot(aes( x = omega2_1, y = reorder(covariates,omega2_1),fill=Pvalue)) + geom_bar(stat="identity", position="dodge")+ meg_stef_theme()+xlab(expression(paste("Partial ",eta^2)))+ylab("Covariates")+ggtitle("Effect sizes on edge count (connectivity)")+scale_fill_paper(discrete = TRUE, palette = "Bi2")

#########
# ASSEMBLE FIGURE SI-6
############

up<-plot_grid(p1,p2, ncol=2,align = "hv",labels=c("A","B"),label_size = 12, label_fontface = "bold", axis="bottom")
low<-plot_grid(p3,p4, ncol=2,align = "hv",labels=c("C","D"),label_size = 12, label_fontface = "bold", axis="bottom")
f<-plot_grid(up, low, nrow=2)
save_plot("SI_figure6.pdf",f,base_height=5, base_asp=1.5)
save_plot("SI_figure6.png",f,base_height=5, base_asp=1.5)







