#Author Stefanie Widder 24062023

library(phyloseq)
library(microbiome)
library(DirichletMultinomial)
library(lattice)
library(xtable)
library(parallel)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
source("../HARMONIZE_FIGURES/1_plot_utils.R")



######ROW MAX for MATRICES
rowMax<-function(m){
sol<-array()
for (i in 1: nrow(m)){
sol[i]<-which.max(m[i,])
}
names(sol)<-rownames(m)
sol
}

#########JACCARD SIMILARITY
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

#############
rarefy<- function(x,min = 0){
 if(min < 0){
    stop("Min should be either 0 or positive.")
 }
 if(min == 0){
    min=min(colsums=apply(x,2,sum))
    print(paste("Rarefy to minimum count",min))
 }else{
 colsums=apply(x,2,sum)
# there are columns below the minimum
if(min(colsums) < min){
keep=c()
# loop column sums
for(j in 1:ncol(x)){
if(colsums[j] >= min){
keep=c(keep,j)
}
}
print(paste("Number of columns",ncol(x)))
print(paste("Keeping ",length(keep)," columns with column sums equal or above",min))
x=x[,keep]
}
 }
 rar=t(rrarefy(t(x),min))
rar
}

#################
aggregate_sets<-function(pseq, nm){
OTU<-otu_table(pseq)
META<-sample_data(pseq)
TAX<-tax_table(pseq)
NO<-OTU
#get classes
S<-as.vector(unique(unname(TAX[,nm])))
for (i in 1: length(S)){
#get asvs
me<-rownames(TAX[TAX[,nm]==S[i],])
#sum up by sample
NO[i,]<-colSums(OTU[me,])
rownames(NO)[i]<-paste(nm,"__",S[i], sep="")
}
#subset to new size
NO<-NO[1:length(S),]
#susbst taxonomy
TAX<-TAX[1:length(S),nm]
TAX[,nm]<-S
rownames(TAX)<-rownames(NO)
phyloseq(otu_table(NO, taxa_are_rows = TRUE), TAX = tax_table(TAX), META = sample_data(META))
}

########

make_p<-function(jj, S4, meta){
z1<-as.data.frame(otu_table(prune_samples(meta$DMM6==jj,S4)))
z1$taxa<-rownames(z1)
#plotting means
z1<-rowMeans(z1[,-ncol(z1)])
z1<-z1/sum(z1)
data.frame(AE=names(z1),rA=as.numeric(unname(z1)),class=rep(paste0("dmm_",jj),length(z1)))
#for plotting distributions
#z1<-z1%>%pivot_longer(1:(ncol(z1)-1),names_to="names",values_to="relab")
#ggboxplot(z1,x="taxa",y="relab", color="taxa")+theme(legend.position="none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


#MAIN

EC<-get(load("EC60.rda"))

#aggregate at S4 level
S4 = aggregate_sets(EC, "S4")



#rarify
tt<-otu_table(S4)
rar<-rarefy(tt)
rar<-rar[rowSums(rar)>0,]

# dim(rar)
#[1] 5 967


#make exacerbation cycle drawing scheme
meta<-sample_data(S4)

sets<-list()
pat<-unique(meta$patient)
for (i in 1:10000){
subs<-array()
for (k in 1: length(pat)){
tmp<-meta[meta$patient==pat[k],]
who<-unique(tmp$EC60)
#draw from patient, add to subset
idx<-sample(1:length(who),1)
subs<-c(subs,who[idx])
}
subs<-subs[-1]
sets[[i]]<-subs
}
length(sets)

#length(sets[!duplicated(sets)])
#36
sets<-sets[!duplicated(sets)]

otu<-rar
counts<-list()
for (k in 1: length(sets)){
ec<-sets[[k]]
me<-rownames(meta[meta$EC60%in%ec,])
zz<-otu[,colnames(otu)%in%me]
#remove empty orgs
zz<-zz[rowSums(zz)>0,]
counts[[k]]<-t(zz)
}

save(counts,file="counts.rda")

fits<-list()
for (k in 1: length(counts)){
fits[[k]] <- mclapply(1:25, dmn, count=counts[[k]], verbose=TRUE)
}
save(fits,file="fitsS4.rda")



#######
#OPTIMAL NUMBER OF DMMs
name<-1:25
my_comparisons<-list(c("1","2"),c("2","3"),c("3","4"),c("4","5"),c("5","6"),c("6","7"),c("7","8"),c("8","9"),c("9","10"),c("10","11"),c("11","12"),c("12","13"),c("13","14"),c("14","15"),c("15","16"),c("16","17"),c("17","18"),c("18","19"),c("19","20"),c("20","21"),c("21","22"),c("22","23"),c("23","24"),c("24","25"))

#LAPLACE
laplace <- sapply(fits,function(x)sapply(x,laplace))
laplace[is.infinite(laplace)]<- NA
dat<-as.data.frame(laplace)
dat<-cbind(dat,k=factor(rownames(dat), levels=name))
test<-dat%>%pivot_longer(cols=1:36,names_to="run",values_to="laplace")

matplot(1:nrow(laplace),cbind(rowMeans(laplace, na.rm=TRUE), laplace), type = "l", lwd= c(2, rep(0.1, ncol(laplace))), col = c(2, rep(1, ncol(laplace))), lty = 1, xlab = "K Dirichlet components", ylab = "Laplace")

##PLOT SFIG1
my_comp<-list(c("1","2"),c("3","4"),c("5","6"),c("6","7"),c("8","9"),c("10","11"),c("12","13"),c("14","15"),c("16","17"),c("18","19"),c("20","21"),c("22","23"),c("24","25"))
s1B<-ggboxplot(test, x="k", y="laplace",col="k", add="jitter", alpha=0.1)+ylab("Laplace")+xlab("K Dirichlet components")+stat_summary(fun=mean, geom="point", shape=23, size=2)+meg_stef_theme() + stat_compare_means(comparisons = my_comp,label = "p.signif")+ggtitle("Laplace approximation for K components")+scale_x_discrete(breaks=seq(from=2,to=25,by=2))+ylim(NA,16200)

#AIC
aic <- sapply(fits,function(x)sapply(x,AIC))
dat<-as.data.frame(aic)
dat<-cbind(dat,k=factor(rownames(dat), levels=name))
test<-dat%>%pivot_longer(cols=1:36,names_to="run",values_to="aic")
pdf("AIC_DMM.pdf")
matplot(1:nrow(aic),cbind(rowMeans(aic), aic), type = "l", lwd= c(2, rep(0.1, ncol(aic))), col = c(2, rep(1, ncol(aic))), lty = 1, xlab = "K Dirichlet components", ylab = "AIC")
ggboxplot(test, x="k", y="aic",col="k", add="jitter", alpha=0.5)+ylab("AIC")+xlab("k")+stat_summary(fun=mean, geom="point", shape=23, size=4)+theme_bw() + stat_compare_means(comparisons = my_comparisons)
dev.off()

#BIC
bic <- sapply(fits,function(x)sapply(x,BIC))
dat<-as.data.frame(bic)
dat<-cbind(dat,k=factor(rownames(dat), levels=name))
test<-dat%>%pivot_longer(cols=1:36,names_to="run",values_to="bic")
pdf("BIC_DMM.pdf")
matplot(1:nrow(bic),cbind(rowMeans(bic), bic), type = "l", lwd= c(2, rep(0.1, ncol(bic))), col = c(2, rep(1, ncol(bic))), lty = 1, xlab = "K Dirichlet components", ylab = "BIC")
ggboxplot(test, x="k", y="bic",col="k", add="jitter", alpha=0.5)+ylab("BIC")+xlab("k")+stat_summary(fun=mean, geom="point", shape=23, size=4)+theme_bw() + stat_compare_means(comparisons = my_comparisons)
dev.off()


#############
#RESULT K=6
dmm<-sapply(fits,function(x)x[[6]])


#########
#IDENTIFY SAMPLE-K ASSOCIATION across all 36 models
###########

#STEP1: UNIQUE LABELS FOR K across all DMM

#get best k for each patient-corrected instant (36 x samples)
kd<-sapply(dmm,function(x)rowMax(mixture(x)))

#make unique labels (36 x samples) & remember dictionary (samples)
kduni<-kd
MK<-max(unlist(kd))

dict<-1:MK
names(dict)<-1:MK
ind<-MK

for (k in 2: length(kduni)){
kduni[[k]]<-kduni[[k]]+ind
dict<-c(dict,1:MK)
names(dict)<- c(names(dict)[names(dict)!=""], ((1:MK)+ind))
ind<-ind+MK
}

#get all sample names
samp<-unique(unlist(sapply(kd,names)))

#prepare result table
allk<-matrix(NA,nrow=length(samp),ncol=36)
rownames(allk)<-samp


#fill table with k 
for (k in 1: length(samp)){
me<-samp[k]
for (l in 1: 36){
tmp<-kduni[[l]]
if(sum(names(tmp)%in%me)>0){
allk[k,l]<-tmp[names(tmp)%in%me]
}else{
allk[k,l]<-NA
}
}
} 

#pivot longer representation 
ak<-as.data.frame(allk)
ak$sample<-rownames(ak)
test<-ak%>%pivot_longer(1:36, names_to = "repetition", values_to = "k")
#remove degk==NA
test<-test[!is.na(test$k),]

#add degenerate k labels 1:MK
degk<-test$k
test$degk<-dict[degk]


# colnames(test)
# "sample"     "repetition" "k"          "degk"      
# dim(test)
# 23064     4



#STEP2: ASSESS JACCARD SIMILARITY  for sample composition of unique Ks across each repetitions (36)
#jaccard similarity: intersect(a,b)/union(a,b) (a=samples per unique k))

#collect samples per k
me<-unique(test$k)
tl<-list()
for (k in 1: length(me)){
tl[[k]]<-unlist(unname(c(test[test$k==me[k],"sample"])))
}
names(tl)<-paste("k_",me, sep="")

#Jaccard for each k aginst each k
J<-matrix(NA,ncol=36*MK,nrow=36*MK)
colnames(J)<-paste("k_",me, sep="")
rownames(J)<-paste("k_",me, sep="")
for (i in 1: length(tl)){
for (j in 1: length(tl)){
J[i,j]<-jaccard(tl[[i]],tl[[j]])
}
}

#STEP3: KMEANS CLUSTERING of individula Ks using Jaccard similarity among them (for K=10 according to laplace)

#kmeans
set.seed(05012022)
kk<-kmeans(J, centers=MK, nstart=10000)
#get k centers
kc<-fitted(kk,method="centers")
#get k assignment
ka<-fitted(kk,method="classes")
#reads: in 36 trials due to patient stratification i find cluster k x times; 
table(ka)

#STEP4: PROBABILITY of K ASSOCIATION per SAMPLE
#normalize against weight(= de facto usage of sample per configuration)
#number of sample x observed per K/(Sum of sample x observations in total)
#over Ks sums to 1 for each sample

#how often was sample x observed=>weight
weight<-table(unlist(tl))

#normalize every sample in K-stratified LIST (K1-10)
who<-names(ka)[ka==1]
tmp<-table(unlist(unname(tl[who])))
res<-tmp/weight[names(tmp)]
rl<-data.frame(sample=names(res),prob=as.numeric(unname(res)),kmeans=rep(1,length(res)))
for (k in 2: MK){
who<-names(ka)[ka==k]
tmp<-table(unlist(unname(tl[who])))
res<-tmp/weight[names(tmp)]
rl<-rbind(rl,data.frame(sample=names(res),prob=as.numeric(unname(res)),kmeans=rep(k,length(res))))
}

#ASSOCIATION TABLE P sample vs K (==mixture())
R<-rl%>%pivot_wider(names_prefix = "kmeans_",names_from=kmeans, values_from=prob, values_fill=0)
##ASSIGNMENT ARRAY best K for sample (==mixture(assig=T))
kass<-rowMax(R[,2:(MK+1)])
names(kass)<-R$sample

##add kmeans association to test
test$kmeans<-kass[test$sample]

#ADD CLASSIFICATION TO EC880
EC<-get(load("../DATA_INPUT/EC880F1B.rda"))
meta<-sample_data(EC)
dms<-kass[rownames(meta)]

#add DMM information to meta information
meta$DMM6<-dms




#####
#SUPPORTING FIGURE 1
#####

##SUPP1A
#Mean DMM TAXON DISTRIBUTION from original data 
S4 = aggregate_sets(EC, "S4")

p<-data.frame(rbind(make_p(1, S4, meta),make_p(2, S4, meta) ))
for (i in 3:MK){
p<-rbind(p, make_p(i, S4, meta))
}
colnames(p)<-c("Aerobicity","relAb","cclass")

p$Aerobicity<-factor(p$Aerobicity,levels=c("S4__C","S4__A", "S4__F", "S4__N", "S4__U"))
p$cclass<-factor(p$cclass,levels=c("dmm_1","dmm_2", "dmm_3", "dmm_4", "dmm_5", "dmm_6"), labels=c("1","2", "3", "4", "5","6"))

s1A<-ggplot(p,aes(x=cclass,y=relAb,fill=Aerobicity))+geom_bar(stat="identity")+scale_fill_manual(values=c("#F1BB7B","#cccccc", "#965A3F","#A4B6E9","#798E87"),labels=c("S4__C"="CF pathogens","S4__A"="Str. aerobes", "S4__F"="Fac. anaerobes", "S4__N"="Str. anaerobes", "S4__U"="Unknown"))+ylab("Mean rel. abundance")+xlab("Community class (DMM component)")+meg_stef_theme()+ggtitle("Mean composition of DMM comunity classes")+theme(legend.key.size=unit(16, 'points')) 





########
#SUPP 1CD
#######



dmm<-paper_pal("DMM")(6)
names(dmm)<-c("1","2","3","4","5","6")
clear<-c("1-pathogens", "2-fac. anaerobes", "3-pathogens","4-anaerobes", "5-anaerobes", "6-anaerobes")
names(clear)<-names(dmm)
#keep original patient/time series color code
coh<-paper_pal("Paper")(11)
names(coh)<-levels(as.factor((as(sample_data(EC880),"data.frame"))$Patient))

cc<-as(sample_data(EC),"data.frame") %>% count(Patient, DMM)
tit<-"Community states by subject"
s1C<-ggplot(cc,aes(x=as.factor(Patient),y=n,fill=DMM)) + geom_bar(position="stack", stat="identity") + meg_stef_theme()+scale_x_discrete(labels=c('Subj. A', 'Subj. B', 'Subj. C', 'Subj. D', 'Subj. E', 'Subj. F','Subj. G', 'Subj. H', 'Subj. I', 'Subj. J', 'Subj. K')) + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0)) + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_fill_manual(values=dmm,name="Community\nstates",labels=clear) + ggtitle(tit) + xlab("Cohort") + ylab("Sample number")+theme(legend.key.size=unit(16, 'points'))

cc<-as(sample_data(EC),"data.frame") %>% count(DMM, Patient)
tit<-"Subjects by community state"

s1D<-ggplot(cc,aes(x=as.factor(DMM),y=n,fill=Patient)) + geom_bar(position="stack", stat="identity") + meg_stef_theme() + guides(color=guide_legend(override.aes=list(fill=NA))) + scale_fill_manual(values=coh,labels=c("Subject A","Subject B","Subject C","Subject D","Subject E", "Subject F","Subject G","Subject H","Subject I","Subject J", "Subject K"),name="Cohort") + ggtitle(tit) + xlab("Community states") + ylab("Sample number")+theme(legend.key.size=unit(16, 'points'))



######
# SUPP1
#####

s1AB<-plot_grid(s1A,s1B+theme(legend.position="none"), labels=c("A","B"),align = "hv", ncol=2, rel_widths=c(1.2,1),label_size = 12, label_fontface = "bold")
s1CD<-plot_grid(s1C,s1D, labels=c("C","D"), align = "hv", ncol=2,label_size = 12, label_fontface = "bold")
sup1<-plot_grid(s1AB,s1CD,nrow=2,align = "v", rel_heights=c(1,1.2))


save_plot("supp_figure1.pdf",sup1, ncol=1, base_height=6.5,  base_asp=1.4)
save_plot("supp_figure1.png",sup1, ncol=1, base_height=6.5,  base_asp=1.4)
