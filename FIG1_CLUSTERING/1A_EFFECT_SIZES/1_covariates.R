library(phyloseq)
library(vegan)
library(MicEco)
library(tidyverse)
library(cowplot)
library(ggalt)
source("../../HARMONIZE_FIGURES/1_plot_utils.R")

##########
#DATA
##########
load("../../DATA_INPUT/EC880F1B.rda")
meta<-data.frame(sample_data(EC880))
#subset meta data
df1<-data.frame(CF2N=meta$CF2N,CFdomRelAb=meta$CFdomRelAb,DMM=as.numeric(meta$DMM),Shannon=meta$Shannon,Chao1=meta$Chao1)

##########
#Distances
##########
#distances ASV
dO<-vegdist(t(otu_table(EC880)), method="bray")
#distances non-standard
d1<-vegdist(scale(df1), method="euc")

#############################
#PERMANOVA with marignalised effect size
##############
set.seed(123)
mmO<-adonis2(dO ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup +Patient,data=meta, by="margin")
mm1<-adonis2(d1 ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup +Patient,data=meta, by="margin")
em1<-adonis_OmegaSq(mm1, partial = TRUE)
emO<-adonis_OmegaSq(mmO, partial = TRUE)
#overal model qualities
totO<-adonis2(dO ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup +Patient,data=meta, by=NULL)
tot1<-adonis2(d1 ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup +Patient,data=meta, by=NULL)

#########
#DUMBBELL PLOT 
########
df<-data.frame(omega2_1=em1$parOmegaSq, p_1=em1$P,omega2_0=emO$parOmegaSq, p_0=emO$P )
row.names(em1)[c(1,3,4,5,6)]<-c("State(B|E)","Age","F508 Zygosity","CFTR","Subject")

df$omega2_1[is.na(df$omega2_1)]<-0
df$omega2_0[is.na(df$omega2_0)]<-0
df$covariates<-row.names(em1)
rownames(df)<-row.names(em1)
df[df<0]<-0
df<-df[-nrow(df),]
df<-df[-nrow(df),]

p1A<-df%>%ggplot(aes(x =omega2_0 , xend = omega2_1, y = reorder(covariates,omega2_0 ))) +
  geom_dumbbell(size = 3.5, color = "#dddddd",  colour_x =mycols("Bpink"), colour_xend =mycols("Tgreen"),dot_guide = TRUE, dot_guide_size = 0.25, show.legend=TRUE) +
  meg_stef_theme()+xlab(expression(omega^2))+theme(axis.title.x = element_text(size = 12))+ylab("Covariates")+ggtitle("Effect sizes")+annotate("text", x = 0.43, y = "Subject", label = "ASV", size = 3, color =mycols("Bpink") )+annotate("text", x = 0.14, y = "Subject", label = "non-std", size = 3, color = mycols("Tgreen"))

save(p1A, file="p1A.rda")
save_plot("figure1A.pdf",p1A, ncol=1, base_height=2, base_asp=1.49)







