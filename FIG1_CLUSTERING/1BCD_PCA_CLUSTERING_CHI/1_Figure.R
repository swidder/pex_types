#Author Stefanie Widder 24062023

library(phyloseq)
library(ggbiplot)
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(corrplot)
library(xlsx)
source("../../HARMONIZE_FIGURES/1_plot_utils.R")
load("../../DATA_INPUT/EC880F1B.rda")


####
#FIX color code
####

#keep patient/age group color code
coh<-paper_pal("Paper")(11)
names(coh)<-levels(as.factor((as(sample_data(EC880),"data.frame"))$Patient))
ag<-paper_pal("Age")(3)
names(ag)<-levels(as.factor((as(sample_data(EC880),"data.frame"))$Agegroup))
ag<-ag[c(3,1,2)]


######
#FIG1B
######

meta<-data.frame(sample_data(EC880))
dat<-data.frame(DMM=as.numeric(meta$DMM),CF2N=meta$CF2N,CFdomRelAb=meta$CFdomRelAb,Shannon=meta$Shannon, Chao1=meta$Chao1)
rownames(dat)<-rownames(meta)


#PCA
p<-prcomp(dat, scale=T, center=T)

#plot
p2<-ggbiplot(p,obs.scale = 1, var.scale = 1,groups = as.matrix(meta)[,"Patient"])+meg_stef_theme()+ guides(color=guide_legend(override.aes=list(fill=NA)))+ggtitle("Ordination non-standard descriptors")+scale_color_paper(discrete = TRUE, palette = "Paper")+geom_point(aes(colour = as.matrix(meta)[,"Patient"]),size=0.5)+xlim(c(-5,3))+ylim(c(-1,3))+theme(legend.position="none")+xlab("PC1(66.4%)")+ylab("PC2(17.9%)")

seg <- which(sapply(p2$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
p2$layers[[seg]]$aes_params$colour <- '#333333'
p2$layers <- c(p2$layers, p2$layers[[seg]])
p2$layers[[seg]]<-NULL
txt <- which(sapply(p2$layers, function(x) class(x$geom)[1] == 'GeomText'))
#label = rownames(p$rotation)
p2$layers[[txt]] <- geom_text(aes(x = xvar, y = yvar,label = c("DMM","pat/an","pat","Shannon", "Chao1"),angle = angle, hjust = hjust),size=3,check_overlap = TRUE,nudge_y=0.01,data = p2$layers[[txt]]$data)
p2$layers <- c(p2$layers, p2$layers[[txt]])
p2$layers[[txt]]<-NULL

write.xlsx2(dat,"source_data_figure1.xlsx",sheetName = "figure1_B", append=TRUE)


#########
#LEGEND Cohort
###########

l2<-get_legend(ggbiplot(p,obs.scale = 1, var.scale = 1,groups = as.matrix(meta)[,"Patient"])+ meg_stef_theme()+ggtitle("Ordination non-standard descriptors")+scale_color_paper(discrete = TRUE, palette = "Paper")+geom_point(aes(colour = as.matrix(meta)[,"Patient"]),size=0.5)+xlim(c(-5,3))+ylim(c(-1,3))+xlab("PC1(66.4%)")+ylab("PC2(17.9%)")+ theme(legend.spacing.y = unit(0.01, 'lines'))+scale_color_discrete(guide = "none")+ scale_colour_manual(values=coh,labels=c("Subject A","Subject B","Subject C","Subject D","Subject E", "Subject F","Subject G","Subject H","Subject I","Subject J", "Subject K"),name="Cohort")+guides(color=guide_legend(ncol=3)))



##########
#FIG1CD
##########

#########
#FIG 1C
##############
an2<-data.frame(Age=meta$Agegroup,Patient=meta$Patient,EC=meta$EC)
rownames(an2)<-rownames(meta)

Cop=paper_pal("Paper")(length(unique(an2$Pat)))
names(Cop)<-unique(an2$Pat)
Aop=paper_pal("Age")(length(unique(an2$Age)))
names(Aop)<-unique(an2$Age)
Eop=paper_pal("CF")(length(unique(an2$EC)))
names(Eop)<-unique(an2$EC)

my_col=list(
Patient=Cop,
Age=Aop
)
res<-pheatmap(p$x[,c(1:3)], scale = "row", clustering_distance_rows = "correlation", annotation_row=an2[,1:2], treeheight_row = 10,cluster_col = F,show_rownames=F,show_colnames=F,annotation_names_row=FALSE,legend=T, annotation_legend=FALSE,annotation_color=my_col,color=colorRampPalette(c("#A597B4", "white",  "#FB6C69"))(500))

p3<-as.ggplot(res)
p3<-p3+heat_theme()+ggtitle("Sample clustering")+xlab("PCs")

write.xlsx2(p$x[,c(1:3)],"source_data_figure1.xlsx",sheetName = "figure1_C", append=TRUE)


#######
#FIG 1D
########
an3<-cbind(an2,k2=cutree(res$tree_row,k=2),k3=cutree(res$tree_row,k=3),k4=cutree(res$tree_row,k=4),k5=cutree(res$tree_row,k=5),k6=cutree(res$tree_row,k=6),k7=cutree(res$tree_row,k=7),k8=cutree(res$tree_row,k=8),k9=cutree(res$tree_row,k=9),k10=cutree(res$tree_row,k=10),k11=cutree(res$tree_row,k=11),k12=cutree(res$tree_row,k=12),k13=cutree(res$tree_row,k=13),k14=cutree(res$tree_row,k=14), k15=cutree(res$tree_row,k=15))


xp<-array()
for (i in 4:17){
for (j in 1:3){
tB<-table(an3[,c(j,i)])
chiB<-chisq.test(tB)
if(j==2){xp[i-3]<-chiB$statistic}

#SUPP FIGURE
if (i==5){
nm<-paste("SI_figure3AB",colnames(as.data.frame(tB))[1],"_",colnames(as.data.frame(tB))[2],".pdf",sep="")
if(j==3){pdf(nm, height=35, width=10)}
if(j==2){pdf(nm, height=20, width=10)}
if(j==1){pdf(nm, height=7, width=10)}
corrplot(chiB$stdres, is.cor = FALSE)
dev.off()
}}}


#RESULTS
#Age groups
        Pearson's Chi-squared test
data:  tB
X-squared = 42.514, df = 4, p-value = 1.305e-08
#Patients(2)
X-squared = 799.85, df = 20, p-value < 2.2e-16
#ECs
X-squared = 839.77, df = 34, p-value < 2.2e-16



dp<-data.frame(k=2:15,e=xp)
p4<-ggplot(dp, aes(x=k,y=e))+geom_line()+geom_point(size=2, shape=21, fill='white')+meg_stef_theme()+ylab(expression(chi^{2}))+xlab("k clusters")+ggtitle("Exacerbation grouping")c
write.xlsx2(an3,"source_data_figure1.xlsx",sheetName = "figure1_D", append=TRUE)


########
##LEGEND AGE
#########

an2$Age<-factor(an2$Age, levels=c("<31","31-37","38-52"))
lA<-get_legend(ggplot(an2,aes(x=1:nrow(an2),y=as.numeric(EC),col=Age))+geom_point()+meg_stef_theme()+scale_color_paper(discrete = TRUE, palette = "Age")+ theme(legend.spacing.y = unit(0.01, 'lines'))+scale_color_discrete(guide = "none")+ scale_colour_manual(values=ag,name="Age group (yrs)"))


#######
#FIG1A
#######
p1<-get(load("../1A_EFFECT_SIZES/p1A.rda"))

#######
#ASSEMBLE FIGURE1
#######
up<-plot_grid(p1, p3,p4, labels=c("A","C","D"), align = "h", ncol=3,label_size = 12, label_fontface = "bold")
leg<-plot_grid( l2, lA, labels=c("",""), ncol=1,label_size = 12, label_fontface = "bold", align="v")
bot<- plot_grid(p2, leg, labels=c("B",""), ncol=2,label_size = 12, label_fontface = "bold")
figure1<-plot_grid(up,bot, labels=c("",""), align = "h", ncol=1,label_size = 12, label_fontface = "bold")

save(figure1, file="figure1.rda")
save_plot("figure1.pdf",figure1, ncol=1, base_height=5)
save_plot("figure1.png",figure1, ncol=1, base_height=5)


