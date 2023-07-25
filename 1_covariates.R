#Author Stefanie Widder 24062023

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
save_plot("fig1A.pdf",p1A, ncol=1, base_height=2, base_asp=1.49)


########
##RESULTS LOOKUP


#totO
#adonis2(formula = dO ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup + Patient, data = meta, by = NULL)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     13  190.760 0.77864 234.32  0.001 ***
#Residual 866   54.231 0.22136                  
#Total    879  244.992 1.00000                  

# tot1
#adonis2(formula = d1 ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup + Patient, data = meta, by = NULL)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     13   2710.8 0.61678 107.22  0.001 ***
#Residual 866   1684.2 0.38322                  
#Total    879   4395.0 1.00000                  


mmO
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = dO ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup + Patient, data = meta, by = "margin")
              Df SumOfSqs      R2        F Pr(>F)    
BETRU          1    0.245 0.00100   3.9116  0.003 ** 
Sex            0    0.000 0.00000      NaN           
Agegroup       2    1.397 0.00570  11.1522  0.001 ***
F508Zygosity   0    0.000 0.00000      Inf           
CFTRgroup      0    0.000 0.00000      Inf           
Patient        5   57.813 0.23598 184.6385  0.001 ***
Residual     866   54.231 0.22136                    
Total        879  244.992 1.00000
mm1
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = d1 ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup + Patient, data = meta, by = "margin")
              Df SumOfSqs      R2       F Pr(>F)    
BETRU          1      3.8 0.00087  1.9769  0.127    
Sex            0      0.0 0.00000    -Inf           
Agegroup       2     14.9 0.00339  3.8323  0.013 *  
F508Zygosity   0      0.0 0.00000    -Inf           
CFTRgroup      0      0.0 0.00000    -Inf           
Patient        5    607.2 0.13816 62.4456  0.001 ***
Residual     866   1684.2 0.38322                   
Total        879   4395.0 1.00000
#####
 emO
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = dO ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup + Patient, data = meta, by = "margin")
              Df SumOfSqs        F parOmegaSq Pr(>F)    
BETRU          1    0.245   3.9116    0.00330  0.003 ** 
Sex            0    0.000      NaN        NaN           
Agegroup       2    1.397  11.1522    0.02255  0.001 ***
F508Zygosity   0    0.000      Inf        NaN           
CFTRgroup      0    0.000      Inf        NaN           
Patient        5   57.813 184.6385    0.51062  0.001 ***
Residual     866   54.231                               
Total        879  244.992                               
---
em1
Permutation test for adonis under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

adonis2(formula = d1 ~ BETRU + Sex + Agegroup + F508Zygosity + CFTRgroup + Patient, data = meta, by = "margin")
               Df SumOfSqs       F parOmegaSq Pr(>F)    
State(B|E)      1      3.8  1.9769   0.001109  0.127    
Sex             0      0.0    -Inf        NaN           
Age             2     14.9  3.8323   0.006396  0.013 *  
F508 Zygosity   0      0.0    -Inf        NaN           
CFTR            0      0.0    -Inf        NaN           
Patient         5    607.2 62.4456   0.258778  0.001 ***
Residual      866   1684.2                              
Total         879   4395.0                              
---









