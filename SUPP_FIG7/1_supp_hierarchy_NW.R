library(tidyverse)
library(phyloseq)
library(fields)
library(igraph)


##AUX
############
# Degree distribution on giant component
############
get.degreeGC<-function(g){
#giant
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#centrality
b<-degree(tmp)
#join to vertex array 
d<-setNames(array(rep(NA,vcount(g))),V(g)$name)
d[names(b)]<-b
d
}

#############
#Clustering coefficient on giant component
#############
get.transitivityGC<-function(g){
#giant
co <- components(g, mode = 'STRONG')
tmp<-induced.subgraph(g, which(co$membership == which.max(co$csize)))
#centrality
b<-transitivity(tmp, type="local")
names(b)<-V(tmp)$name
#join to vertex array 
d<-setNames(array(rep(NA,vcount(g))),V(g)$name)
d[names(b)]<-b
d
}



#MAIN
#data
load("../DATA_INPUT/EC880F1B.rda")
TAX<-as.data.frame(tax_table(EC880))


#representative NWs
me<-c(10,9,38)
K3<-c("PAT","AN1","AN2")
targ<-c("../DATA_INPUT/NWs/EC16_PATIENT_J_SS/NWs_EC16_patientJ.rda","../DATA_INPUT/NWs/EC15_PATIENT_I_SS/NWs_EC15_patientI.rda","../DATA_INPUT/NWs/EC15_PATIENT_I_SS/NWs_EC15_patientI.rda" )
m<-c("ff9d93d7b7e46787568f2d241caeaf3b","06f825b512d903b9230e1a55d87359ee","32f8fd11d2bee278d609a1d4ab767554")

#read in selected NWs
g<-list()
for (i in 1:length(targ)){
#get nws by EC
#read NW
NW<-get(load(targ[i]))
g[[i]]<-NW$fullNW[[me[i]]]
names(g)[[i]]<-names(NW$fullNW[me[i]])
}


#plot NWs hierarchical (Sugiyama) 
set.seed(123)
lab<-"LAB"
for (i in 1:length(g)){
g1<-g[[i]]

#fix graph attributes for plotting
g1<-delete_edge_attr(g1,"color")
E(g1)$weight<-abs(E(g1)$weight)

#select largest component
co <- components(g1, mode = 'STRONG')
lg<-induced.subgraph(g1, which(co$membership == which.max(co$csize)))

#get vertex clustering coefficient and degree
cc<-get.transitivityGC(lg)
cc[is.na(cc)]<-0
lg<-set_vertex_attr(lg,"cc",index = V(lg),value=cc)
lg<-set_vertex_attr(lg,"deg",index = V(lg),value=get.degreeGC(lg))

#clear names by ASV groups
nms1<-TAX[V(lg)$name,"S4"]
nms1[is.na(nms1)]<-"U"
#pathogens
me<-which(nms1=="C")
#hierarchy winners
win<-which(V(lg)$name %in%m[i])
#replace pathogen by name
nms2<-unname(sapply(TAX[V(lg)$name[me],"Genus"],function(x)strsplit(x,"__")[[1]][2]))
#replace winners by name
nms3<-unname(sapply(TAX[V(lg)$name[win],"Genus"],function(x)strsplit(x,"__")[[1]][2]))
if(lab=="LAB"){V(lg)$name<-nms1}
if(lab!="LAB"){
V(lg)$name<-""
if(length(nms3)>0){V(lg)$name[win]<-nms3}
if(length(nms2)>0){V(lg)$name[me]<-nms2}
}

#make layout
llg<-layout_with_sugiyama(lg,attributes="all",hgap=10,vgap=10)
#node size
size<-V(lg)$deg
size[is.na(size)]<-1
size<-size*3
ff <- factor(c(min(V(lg)$deg),max(V(lg)$deg)))

#node color clustering
hm<-names(table(cc))
col<-heat.colors(length(hm)-1)
col<-c("grey",col)
names(col)<-hm
#node color ASV groups
colA=c( "#E69F00" ,"purple2","#CC3366","#56B4E9","#C0C0C0")
names(colA)<-c("C","A","F","N","U")

if(lab!="LAB"){
#make file name 
f=paste(K3[i],"noLab.pdf",sep="")
pdf(f)
#label offset
if(i==1)plot(llg$extd_graph, vertex.size=size,vertex.color=col[as.character(V(lg)$cc)], vertex.label.cex=1,vertex.label.color="black")
if(i>1)plot(llg$extd_graph, vertex.size=size,vertex.color=col[as.character(V(lg)$cc)], vertex.label.color="black",vertex.label.cex=1,vertex.label.dist=as.numeric(V(lg)$name=="Pseudomonas")*1.5,vertex.label.degree=pi/2)
legend("topleft", legend = levels(ff), pch = 1, pt.cex=c(1,3), bty = "n", title="Degree")
image.plot(legend.only=T, zlim=range(as.numeric(names(col))), col=col, legend.args = list( text = "Clustering\nCoefficient", cex = 1, side = 3,line = .5))
dev.off()
}

if(lab=="LAB"){
f=paste(K3[i],"Lab.pdf",sep="")
pdf(f)
#label offset
plot(llg$extd_graph, vertex.label.cex=1,vertex.label.color="black",vertex.color=colA[V(lg)$name])
dev.off()
}

#png
if(lab!="LAB"){
f=paste(K3[i],"noLab.png",sep="")
png(f, res = 150, height = 800, width = 815)
#label offset
if(i==1)plot(llg$extd_graph, vertex.size=size,vertex.color=col[as.character(V(lg)$cc)], vertex.label.cex=1,vertex.label.color="black") 
if(i>1)plot(llg$extd_graph, vertex.size=size,vertex.color=col[as.character(V(lg)$cc)], vertex.label.cex=1,vertex.label.dist=as.numeric(V(lg)$name=="Pseudomonas")*1.5,vertex.label.degree=pi/2,vertex.label.color="black")
legend("topleft", legend = levels(ff), pch = 1, pt.cex=c(1,3), bty = "n", title="Degree")
image.plot(legend.only=T, zlim=range(as.numeric(names(col))), col=col, legend.args = list( text = "Clustering\nCoefficient", cex = 1, side = 3,line = .5))
dev.off()
}

if(lab=="LAB"){
f=paste(K3[i],"Lab.png",sep="")
png(f, res = 150, height = 800, width = 850)
#label offset
plot(llg$extd_graph, vertex.label.cex=1,vertex.label.color="black",vertex.color=colA[V(lg)$name] ) 
dev.off()
}

}

