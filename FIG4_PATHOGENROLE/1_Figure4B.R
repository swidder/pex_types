library(phyloseq)
library(seqtime)
library(dplyr)
library(tidyr)
library(ggpubr)
source("../HARMONIZE_FIGURES/1_plot_utils.R")

############
#AUX
#############
get_color_counts<-function(dat, sample, color){
if(length(color)==0){return (0)}
sum(dat[color,sample])
}

####calc relative abundance of noise types for plotting
relab<-function(table){
all<-sum(length(table$white),length(table$pink),length(table$brown), length(table$black))
c(length(table$white)/all, length(table$pink)/all,length(table$brown)/all,length(table$black)/all) 
}

###################
returns df with counts and real abundance of colors by organims group (C, N, F, A)
get_org_color_df<-function(dfC, S4S, NTT, tax, group){
dfC$whiteC<-NA
dfC$pinkC<-NA
dfC$brownC<-NA
dfC$blackC<-NA
dfC$white<-NA
dfC$pink<-NA
dfC$brown<-NA
dfC$black<-NA

#by EC
uu<-length(unique(dfC$ECid))
for (i in 1: uu){
#by sample
for (j in 1: ncol(S4S[[i]])){
#Counts
tot<-0
tot<-tot+ (dfC[colnames(S4S[[i]])[j],"whiteC"]<-get_color_counts(S4S[[i]],colnames(S4S[[i]])[j],NTT[[i]]$white[tax[NTT[[i]]$white,"S4"]%in%group]))
tot<-tot+(dfC[colnames(S4S[[i]])[j],"pinkC"]<-get_color_counts(S4S[[i]],colnames(S4S[[i]])[j],NTT[[i]]$pink[tax[NTT[[i]]$pink,"S4"]%in%group]))
tot<-tot+(dfC[colnames(S4S[[i]])[j],"brownC"]<-get_color_counts(S4S[[i]],colnames(S4S[[i]])[j],NTT[[i]]$brown[tax[NTT[[i]]$brown,"S4"]%in%group]))
tot<-tot+(dfC[colnames(S4S[[i]])[j],"blackC"]<-get_color_counts(S4S[[i]],colnames(S4S[[i]])[j],NTT[[i]]$black[tax[NTT[[i]]$black,"S4"]%in%group]))
#relAb
if(tot>0){
dfC[colnames(S4S[[i]])[j],"white"]<-dfC[colnames(S4S[[i]])[j],"whiteC"]/tot
dfC[colnames(S4S[[i]])[j],"pink"]<-dfC[colnames(S4S[[i]])[j],"pinkC"]/tot
dfC[colnames(S4S[[i]])[j],"brown"]<-dfC[colnames(S4S[[i]])[j],"brownC"]/tot
dfC[colnames(S4S[[i]])[j],"black"]<-dfC[colnames(S4S[[i]])[j],"blackC"]/tot
}else{
dfC[colnames(S4S[[i]])[j],"white"]<-dfC[colnames(S4S[[i]])[j],"pink"]<-dfC[colnames(S4S[[i]])[j],"brown"]<-dfC[colnames(S4S[[i]])[j],"black"]<-0
}}}
dfC<-na.omit(dfC)
dfC
}


#############
#MAIN
############
########
#LOAD DATA
#######
YY<-get(load("../DATA_INPUT/EC880F1B.rda"))
#remove excluded time series
me<-sample_names(EC880)[!is.na(sample_data(EC880)$PExClust)]
YY<-prune_samples(me,EC880)

#split to sub files
otu<-otu_table(YY)
met<-sample_data(YY)
tax<-as.data.frame(tax_table(YY))
met$SampleID<-rownames(met)
met$day<-NA

#harmonize time series names
ats<-c(1, 10, 11, 13, 14, 15, 16, 17, 18,  2,  3,  4,  5,  6,  7,  9 ) 
M<-list()
S<-list()
for (i in 1:length(ats)){
M[[i]]<-met[met$ECid==ats[i],c(22,2,23,3,8)]
S[[i]]<-otu[,met$ECid==ats[i]]
}
names(S)<-names(M)<-paste("EC",ats[1:16],sep="_")

#remove empty ASVs
for (i in 1:length(ats)){
S[[i]]<-S[[i]][rowSums(S[[i]])>0,]
}

#########
#PREDICT NOISE COLORS by ASV
#########

#interpolate timepoints
SS<-list()
for (i in 1:length(ats)){
SS[[i]]<-interpolate(S[[i]],time.vector=rev(M[[i]]$SampleDay2T)-(rev(M[[i]]$SampleDay2T)[1]-1))
}

#round to integer, remove <0 from interpolation
SS<-lapply(SS,function(x)replace(x,x<0,0))
SS<-lapply(SS,function(x)round(x,digit=0))

#calculate noise colors
NT<-list()
for(i in 1:length(ats)){
NT[[i]]<-identifyNoisetypes(SS[[i]], smooth=TRUE)
}

#######
#PLOT 
#######

#clearID
tax$clearID<-gsub("_s__","",gsub("g__","",paste( tax[,6], tax[,7], 1:nrow(tax), sep="_")))

#prepare lists/dfs
NTT<-NT
for (i in 1: length(ats)){
NTT[[i]]$white<-rownames(SS[[i]])[NTT[[i]]$white]
NTT[[i]]$pink<-rownames(SS[[i]])[NTT[[i]]$pink]
NTT[[i]]$brown<-rownames(SS[[i]])[NTT[[i]]$brown]
NTT[[i]]$black<-rownames(SS[[i]])[NTT[[i]]$black]
}

#merge noise colors and ASV abundance
dfC<-data.frame(met[,c(21,3,19,29,8)])
dfC<-get_org_color_df(dfC, S, NTT, tax, "C")
dfN<-data.frame(met[,c(21,3,19,29,8)])
dfN<-get_org_color_df(dfN, S, NTT, tax, "N")

#calculate abundance ratio
dfC$white2N<-dfC$whiteC/dfN$whiteC
dfC$pink2N<-dfC$pinkC/dfN$pinkC

tt<-data.frame(SampleID=unlist(dfC$SampleID),ExacerbationRegime=unlist(dfC$PExClust),white=unlist(dfC$white2N),pink=unlist(dfC$pink2N),D23=unlist(dfC$D23))

#remove empty samples
tt[tt==0]<-NA
tt[tt==Inf]<-NA

# transform to wider format
ha<-tt%>%pivot_longer(cols=c(3,4),names_to="color")
ha<-na.omit(ha)
#remove outliers
ha<-ha[ha$value<150,]

######
#PLOT
######

cfcol<-c("PAT"="#F1BB7B","AN1"="#FD6467","AN2"="#5B1A18")
ha$ExacerbationRegime<-factor(ha$ExacerbationRegime, levels=c("PAT","AN1", "AN2"))
ha$color<-factor(ha$color, levels=c("white", "pink"))

p4B<-ggplot(ha,aes(x=ExacerbationRegime, y=log2(value), color=as.factor(color)))+geom_hline(yintercept=0,col="#02401B",linetype="dashed")+geom_boxplot(outlier.shape = NA)+geom_point(shape=16, position=position_jitterdodge(jitter.width=NULL, jitter.height=0,dodge.width=0.75), size=0.7, alpha=0.5)+ylab(expression(Log[2](pat/an)))+xlab("PEx types")+ meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA)))+ scale_color_paper(discrete=T,palette="Noise",name="Noise\ncolor",breaks=c('pink','white'))+ggtitle("Pathogen behavior in distinct community backgrounds")

#save plots
save(p4B,file="p4B.rda")


