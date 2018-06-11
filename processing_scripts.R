clusters<-read.csv('~/Desktop/clusters.csv')
markers<-read.delim('~/Desktop/miRNA/POPULATIONS_Populist_4261_Ideal_Vector_correlation.new.names.txt')
geneClusters<-clusters[-seq(514),]
markers<-markers[which(markers$X%in%geneClusters$Nodes),]
markersNames<-markers$X
markers<-markers[,-1]
hist(as.vector(as.matrix(markers)),breaks = 100, xlab = 'Gene Enrichment Score',main = 'Distribution of Cell-Type-Specificty Scores of Genes')
treshold = mean(as.matrix(markers))+3*sqrt(var(as.vector(as.matrix(markers))))
abline(v = treshold,col = 'red')
dep_treshold = mean(as.matrix(markers))-3*sqrt(var(as.vector(as.matrix(markers))))
abline(v = dep_treshold, col = 'blue')
text(x = 2*treshold, y= 15000,'3*stdev')
text(x = 2.5*treshold, y= 10000,'Enrichment', col = 'red')
text(x = -2*treshold, y= 15000,'3*stdev')
text(x = -2.5*treshold, y= 10000,'Depletion', col = 'blue')

markers<-markers[,which(colSums(markers>treshold)>0)]
enrich<-list()
L.enr<-list()
for( i in seq(22)){
  for (j in seq(dim(markers)[2])){
    enrGene<-markersNames[which(markers[,j] > treshold)]
    L.enr[[j]]<-enrGene[which(enrGene%in%geneClusters$Nodes[which(geneClusters$GW.clusterID==i)])]

  }
  enrich[[i]]<-L.enr
}

enrMat<-matrix(0,22,dim(markers)[2])
p.enr<-enrMat
Odds<-p.enr
for ( i in seq(22)){
  for (j in seq(dim(markers)[2])){
    enrMat[i,j]=length(unlist(enrich[[i]][j]))
    list1<-sum(geneClusters$GW.clusterID==i)
    list2<-sum(markers[,j] > treshold)
    #print(list2)
    #PopSize<-3439
    PopSize<-12979
    #PopSize<-length(markersNames)
    p.enr[i,j]=phyper(enrMat[i,j]-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
    Odds[i,j]=log2(enrMat[i,j]/(list2-enrMat[i,j]))
  }
}

colnames(p.enr)=colnames(markers)
rownames(p.enr)=labels2colors(seq(22))


par(cex.main = 0.8)
require(gplots)
require(WGCNA)
my_palette<-colorRampPalette(c('white','red'))
txt<-ifelse(p.enr<0.1,yes = ifelse(p.enr<0.05,yes = ifelse(p.enr<0.01,yes = '***',no = '**'),no = '*'),no = '')
heatmap.2(-log10(p.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', key.xlab="-log10(p-value)",cellnote = txt)

heatmap.2(p.enr,trace = 'none',cexCol  = 0.8,cexRow = 0.8,main = 'cell-type-enrichment',srtCol = 45, key.xlab="p value")

q.enr=p.enr
b = dim(p.enr)[1]
for(i in seq(b)){
  q.enr[,i]=p.adjust(p.enr[,i],method = 'bonferroni')
}

txt<-ifelse(q.enr<0.1,yes = ifelse(q.enr<0.05,yes = ifelse(q.enr<0.01,yes = '***',no = '**'),no = '*'),no = '')
#heatmap.2(-log10(q.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
#          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', key.xlab="-log10(Bonferroni corrected p-value )",cellnote = txt)
#heatmap.2(-log10(q.enr),trace = 'none',cexCol  = 0.6,cexRow = 0.6,col = my_palette,
#          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', key.xlab="-log10(Bonferroni corrected p-value )",colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),sepcolor="gray",sepwidth=c(0.001,0.001))
#heatmap.2(-log10(q.enr)[,-grep(colnames(q.enr),pattern = 'Unk')],trace = 'none',cexCol  = 0.6,cexRow = 0.6,col = my_palette,
#          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', key.xlab="-log10(Bonferroni corrected p-value )",colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),sepcolor="gray",sepwidth=c(0.001,0.001))
my_palette<-colorRampPalette(c('white','red','purple'))
colnames(q.enr)
myOrder<-c('turquoise','blue','lightgreen','grey60','greenyellow','tan','brown',
           'black','pink','cyan','salmon','lightyellow','green','purple','yellow',
           'red','royalblue','lightcyan','magenta','darkred','darkgreen','midnightblue')
heatmap.2(-log10(q.enr)[myOrder,as.vector(tmp)],trace = 'none',cexCol  = 0.9,cexRow = 0.8,col = my_palette,
          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', 
          key.xlab="-log10(Bonferroni corrected p-value )",colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),
          sepcolor="gray",sepwidth=c(0.001,0.001),dendrogram = 'none', Rowv = FALSE, Colv=FALSE)

#__________________________________________________________________________
#### Enrichment in Sousa AMM, et al. Science modules.
data<-read.csv('~/Desktop/results/Ago2/science_modules.csv')

#head(data)
colnames(data)
idx<-which(!is.na(data$Modules))
for(i in 1:(length(idx)-1)){
  data$Modules[idx[i]:(idx[i+1]-1)]=data$Modules[idx[i]]
}
data$Modules[which(is.na(data$Modules))]=max(data$Modules,na.rm = TRUE)


a = max(geneClusters$GW.clusterID)
b = max(data$Modules)

for( i in seq(a)){
  for (j in seq(b)){
    enrGene<-data$Gene.names[which(data$Modules==j)]
    L.enr[[j]]<-enrGene[which(enrGene%in%geneClusters$Nodes[which(geneClusters$GW.clusterID==i)])]
    
  }
  enrich[[i]]<-L.enr
}

enrMat<-matrix(0,a,b)
p.enr<-enrMat
Odds<-p.enr
for ( i in seq(a)){
  for (j in seq(b)){
    enrMat[i,j]=length(unlist(enrich[[i]][j]))
    list1<-sum(geneClusters$GW.clusterID==i)
    list2<-sum(data$Modules==j)
    #print(list2)
    #PopSize<-3439
    PopSize<-12979
    p.enr[i,j]=phyper(enrMat[i,j]-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
    Odds[i,j]=log2(enrMat[i,j]/(list2-enrMat[i,j]))
  }
}

colnames(p.enr)=1:max(data$Modules)
rownames(p.enr)=labels2colors(seq(a))


par(cex.main = 0.8)
require(gplots)
require(WGCNA)
my_palette<-colorRampPalette(c('white','red'))
txt<-ifelse(p.enr<0.1,yes = ifelse(p.enr<0.05,yes = ifelse(p.enr<0.01,yes = '***',no = '**'),no = '*'),no = '')
heatmap.2(-log10(p.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', key.xlab="-log10(p-value)",cellnote = txt)
heatmap.2(-log10(p.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'cell-type-enrichment',srtCol = 45,notecol = 'black', key.xlab="-log10(p-value)")
heatmap.2(p.enr,trace = 'none',cexCol  = 0.8,cexRow = 0.8,main = 'cell-type-enrichment',srtCol = 45, key.xlab="p value")


q.enr=p.enr
KEEP=seq(b)
for(i in seq(b)){
  q.enr[,i]=p.adjust(p.enr[,i],method = 'fdr')
  KEEP[i]<-sum(min(q.enr[,i])<0.1)
}

txt<-ifelse(q.enr<0.1,yes = ifelse(q.enr<0.05,yes = ifelse(q.enr<0.01,yes = '***',no = '**'),no = '*'),no = '')
heatmap.2(-log10(q.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",cellnote = txt)

FDR.sig<-q.enr[,which(KEEP==1)]
my_palette<-colorRampPalette(c('white','red','purple'))
txt<-ifelse(FDR.sig<0.1,yes = ifelse(FDR.sig<0.05,yes = ifelse(FDR.sig<0.01,yes = '***',no = '**'),no = '*'),no = '')
heatmap.2(-log10(FDR.sig),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,colsep=1:ncol(FDR.sig),rowsep=1:nrow(FDR.sig),
          sepcolor="gray",sepwidth=c(0.001,0.001),
          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",cellnote = txt)

Data<-read.csv('Desktop/results/Ago2/science_modules.csv')
Data$Cell.type.enrichment..Single.cell.RNA.seq.[which(Data$Modules %in% as.numeric(colnames(FDR.sig)))]

### Autism Enrichment
data<-read.delim('~/Desktop/ASD.txt',header = FALSE)
ASD<-data$V1
p.ATS<-seq(22)
for(i in seq(22)){
  overlap<-sum(ASD%in%clusters$Nodes[clusters$GW.clusterID==i])
  list1<-sum(clusters$GW.clusterID==i)
  list2<-length(ASD)
  PopSize<-12979
  p.ATS[i]=phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
}

barplot(ylim = c(0,1.8),-log10(p.adjust(p.ATS)),names.arg = labels2colors(seq(22)),las = 2,cex.names = 0.8, ylab = '-log10(FDR)', main = 'Enrichment in Autism')
text(10.2,y =1.5,'**',cex = 1.5)
abline(h = -log10(0.05),col = 'red', lty = 5)

data<-read.delim('Desktop/ASD.txt',header = FALSE)
ASD<-data$V1[c(which(data$V2==2),which(data$V2==1))]
p.ATS<-seq(22)
for(i in seq(22)){
  overlap<-sum(ASD%in%clusters$Nodes[clusters$GW.clusterID==i])
  list1<-sum(clusters$GW.clusterID==i)
  list2<-length(ASD)
  PopSize<-12979
  p.ATS[i]=phyper(overlap-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
}

barplot(ylim = c(0,1.3),-log10(p.adjust(p.ATS)),names.arg = labels2colors(seq(22)),las = 2,cex.names = 0.8, ylab = '-log10(FDR)', main = 'Enrichment in Autism')
text(10.2,y =1.1,'*',cex = 1.5)
#abline(h = -log10(0.1),col = 'red', lty = 5)



#__________________________________________________________________________
####  Bipartite module detection on Sousa AMM, et al. Science.
#myData<-read.csv('~/Desktop/results/Ago2/Nenads science paper/interactions.csv')
myData<-read.csv('~/Desktop/results/adult_HISTCLIP.csv')
clusters<-read.csv('~/Desktop/clusters.csv')
#myData=myData[-c(1,2),c(3,2)]
myData<-myData[,-3]
myData$val<-TRUE
require(reshape2)
colnames(myData)<-c('miR','geneID','val')
myData$geneID<-gsub(".*\\|","",myData$geneID)
web <- reshape2::dcast(data=myData, miR ~ geneID, value.var="val")
web[is.na(web)]<-FALSE
web.mat<-as.matrix(sapply(web, as.numeric))  
web.mat<-web.mat[,-1]
row.names(web.mat)=web[,1]
web.mat2<-web.mat
web.mat2[web.mat > 0] <- 1

require(lpbrim)
N<-web.mat2
res = list()
for(i in seq(100)){
  res[[i]]<-findModules(N,50,sparse = FALSE)
}
a<-dim(N)[1]+dim(N)[2]
myMat = matrix(data = 0,nrow = a,ncol = a)
for(i in seq(a)){
  for(j in seq(a)){
    for(k in seq(100)){
      cls<-as.numeric(which(res[[k]]$S[i,]==1))
      myMat[i,j]<-myMat[i,j]+res[[k]]$S[j,cls]
    }
  }  
}
myMat = myMat/100
rownames(myMat)<-rownames(res[[1]]$S)
colnames(myMat)<-rownames(res[[1]]$S)

save.image()
heatmap(myMat)

dM<-as.matrix(dist(myMat,method = "euclidean"))
hc <- hclust(d = as.dist(dM), method = "ave")

minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = hc, distM = dM,
                            deepSplit = 2
                            ,method = 'hybrid', pamRespectsDendro= FALSE,
                            minClusterSize = minModuleSize); 
dynamicColors=labels2colors(dynamicMods)

plotDendroAndColors(hc, dynamicColors, "Dynamic Tree Cut", hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = " modules of Sousa AMM, et al. Science network")
df<-data.frame(colnames(myMat),dynamicColors)
colnames(df)<-c('Gene.names','Modules')

write.csv(df,file = '~/Desktop/results/Ago2/Nenads science paper/Sousa_modules.csv')
data<-df
data$Modules<-as.numeric(dynamicMods)

#data<-data[-grep(pattern = 'hsa',data$Gene.names),]
geneClusters<-clusters
#geneClusters<-clusters[-seq(514),]
a = max(geneClusters$GW.clusterID)
b = max(data$Modules)

for( i in seq(a)){
  for (j in seq(b)){
    enrGene<-data$Gene.names[which(data$Modules==j)]
    L.enr[[j]]<-enrGene[which(enrGene%in%geneClusters$Nodes[which(geneClusters$GW.clusterID==i)])]
    
  }
  enrich[[i]]<-L.enr
}

enrMat<-matrix(0,a,b)
p.enr<-enrMat
Odds<-p.enr
for ( i in seq(a)){
  for (j in seq(b)){
    enrMat[i,j]=length(unlist(enrich[[i]][j]))
    list1<-sum(geneClusters$GW.clusterID==i)
    list2<-sum(data$Modules==j)
    #print(list2)
    #PopSize<-3439
    PopSize<-12979
    p.enr[i,j]=phyper(enrMat[i,j]-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
    Odds[i,j]=log2(enrMat[i,j]/(list2-enrMat[i,j]))
  }
}

colnames(p.enr)=1:max(data$Modules)
rownames(p.enr)=labels2colors(seq(a))

q.enr=p.enr
KEEP=seq(b)
for(i in seq(b)){
  q.enr[,i]=p.adjust(p.enr[,i],method = 'fdr')
  KEEP[i]<-sum(min(q.enr[,i])<0.1)
}
FDR.sig<-q.enr[,which(KEEP==1)]
my_palette<-colorRampPalette(c('white','red','purple'))
txt<-ifelse(FDR.sig<0.1,yes = ifelse(FDR.sig<0.05,yes = ifelse(FDR.sig<0.01,yes = '***',no = '**'),no = '*'),no = '')
heatmap.2(-log10(FDR.sig),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",cellnote = txt,colsep=1:ncol(FDR.sig),rowsep=1:nrow(FDR.sig),
          sepcolor="gray",sepwidth=c(0.001,0.001))
txt<-ifelse(q.enr<0.1,yes = ifelse(q.enr<0.05,yes = ifelse(q.sig<0.01,yes = '***',no = '**'),no = '*'),no = '')
heatmap.2(-log10(q.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",cellnote = txt,colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),sepcolor="gray",sepwidth=c(0.001,0.001),dendrogram='none',     
          Rowv=FALSE,
          Colv=FALSE)



#__________________________________________________________________________
### miRNA DE calling ###
require(plyr)
clusters<-read.csv('~/Desktop/results/Tom/clusters.csv')
tmp.data<-read.csv('~/Desktop/results/Tom/microRNA_assay_names.csv',header = FALSE)
data<-read.delim('~/Desktop/results/Tom/repeat.mir.annotated.txt')
colnames(data)<-gsub("\\.", "-",colnames(data))
colnames(data)<-gsub("_", "-",colnames(data))
tmp.data[,1]<-gsub("_", "-",as.character(tmp.data[,1]))
tmp.data[,2]<-gsub("_", "-",as.character(tmp.data[,2]))
colnames(data)<-mapvalues(colnames(data),as.character(tmp.data[,1]),as.character(tmp.data[,2]))
tmp<-read.delim('~/Desktop/results/Tom/gw18_clusteridentity.txt')
data$Cluster<-tmp[,2]
write.csv(file = '~/Desktop/results/Tom/repeat.mir.annotated.updated.csv',x = data)

types<-read.csv('~/Desktop/results/Tom/revised_cluster_interpretation.csv',header = FALSE)$V3
g<-as.numeric(t(data$Cluster))
b<-seq(11)
b<-b[-grep(types,pattern = 'low')]
data<-data[-which(data$Cluster%in%grep(types,pattern = 'low')),1:87] # remove low quality cells and filter for miRNAs
g<-as.numeric(t(data$Cluster))
g<-mapvalues(g,b,seq(8))
table(data$Cluster)
colnames(data)
M<-t(data[,3:dim(data)[2]])
#mat<-t(t(M)/colSums(M)) #normalize
mat<-M
rownames(mat)<-colnames(data[,3:dim(data)[2]])
colnames(mat)<-rownames(data[,3:dim(data)[2]])
n<-max(g)
m<-dim(mat)[1]
p<-matrix(0,m,n)
for(i in seq(n)){
  idx<-which(ifelse(g==i,TRUE,FALSE))
  for(j in seq(m)){
    p[j,i]<-(wilcox.test(mat[j,idx],mat[j,-idx], alternative = "greater"))$p.value
  }
  p[,i]<-p.adjust(p[,i],method = 'bonferroni')
}

rownames(p)<-row.names(mat)
colnames(p)<-types[-grep(types,pattern = 'low')]
q.enr<-p
#txt<-ifelse(q.enr<0.01,yes = ifelse(q.enr<0.001,yes = ifelse(q.enr<0.0001,yes = '***',no = '**'),no = '*'),no = '')
#heatmap.2(-log10(q.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
#          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",dendrogram='none',cellnote = txt)
#num_mrk<-5
#top.markers<-matrix(0,num_mrk,n)
#for(i in seq(n)){top.markers[,i]<-order(q.enr[,i])[1:num_mrk]}
#sub.idx<-as.vector(top.markers)
#top.q<-q.enr[sub.idx,]
#txt<-ifelse(top.q<0.05,yes = ifelse(top.q<0.01,yes = ifelse(top.q<0.001,yes = '***',no = '**'),no = '*'),no = '')
#my_palette<-colorRampPalette(c('white','orange','red','purple'))
#heatmap.2(-log10(top.q),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
#          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",dendrogram='none',Rowv=FALSE,
#          Colv=FALSE,cellnote = txt,colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),sepcolor="gray",sepwidth=c(0.001,0.001))

l<-which(as.vector(apply(q.enr,MARGIN = 1,FUN=function(x){sum(x<0.1)>0}))) #filter for significant p values
top.q<-q.enr[l,]
txt<-ifelse(top.q<0.05,yes = ifelse(top.q<0.01,yes = ifelse(top.q<0.001,yes = '***',no = '**'),no = '*'),no = '')
my_palette<-colorRampPalette(c('white','orange','red','purple'))
par(cex.main=0.8)
require(gplots)
heatmap.2(-log10(top.q),trace = 'none',cexCol  = 0.8,cexRow = 0.45,col = my_palette,
          main = 'cell-type-enrichment analysis',srtCol = 45,notecol = 'black', key.xlab="-log10(Bonferroni corrected p value)",dendrogram='none',Rowv=FALSE,
          Colv=FALSE,cellnote = txt,colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),sepcolor="gray",sepwidth=c(0.001,0.001),lhei=c(0.2,1), lwid=c(2,4.5), keysize=0.75, key.par = list(cex=0.5))
heatmap.2(-log10(top.q),trace = 'none',cexCol  = 0.8,cexRow = 0.45,col = my_palette,
          main = 'cell-type-enrichment analysis',srtCol = 45,notecol = 'black', key.xlab="-log10(Bonferroni corrected p value)", cellnote = txt,colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),
          sepcolor="gray",sepwidth=c(0.001,0.001),lhei=c(0.2,1), lwid=c(2,4.5), keysize=0.75, key.par = list(cex=0.5))

#__________________________________________________________________________
### WGCNA overlap with bipartite ###
data<-read.delim('~/Desktop/results/module.targets.txt') #read WGNCA modules
colnames(data)<-c('Modules','Gene.names')
#data$Gene.names<-gsub("\\.", "-",data$Gene.names)
#data$Gene.names<-gsub("_", "-",data$Gene.names)
data$Modules<-tolower(data$Modules)
clusters<-read.csv('~/Desktop/clusters.csv')
#clusters$Nodes<-gsub("_", "-",clusters$Nodes)
#clusters<-clusters[which(clusters$Nodes%in%data$Gene.names),]

tmp<-seq(dim(table(data$Modules)))-1
data$Modules<-as.numeric(mapvalues(data$Modules,labels2colors(tmp),tmp))
data$Modules=data$Modules+1
a = max(geneClusters$GW.clusterID)
b = max(data$Modules)
for( i in seq(a)){
  for (j in seq(b)){
    enrGene<-data$Gene.names[which(data$Modules==j)]
    L.enr[[j]]<-enrGene[which(enrGene%in%geneClusters$Nodes[which(geneClusters$GW.clusterID==i)])]
  }
  enrich[[i]]<-L.enr
}

enrMat<-matrix(0,a,b)
p.enr<-enrMat
Odds<-p.enr
for ( i in seq(a)){
  for (j in seq(b)){
    enrMat[i,j]=length(unlist(enrich[[i]][j]))
    list1<-sum(geneClusters$GW.clusterID==i)
    list2<-sum(data$Modules==j)
    #print(list2)
    PopSize<-3439
    #PopSize<-12979
    p.enr[i,j]=phyper(enrMat[i,j]-1,list1,PopSize-list1,list2,lower.tail = FALSE, log.p = FALSE)
    Odds[i,j]=log2(enrMat[i,j]/(list2-enrMat[i,j]))
  }
}

colnames(p.enr)=1:max(data$Modules)
rownames(p.enr)=labels2colors(seq(a))

q.enr=p.enr
#KEEP=seq(b)
for(i in seq(b)){
  q.enr[,i]=p.adjust(p.enr[,i],method = 'bonferroni')
  #KEEP[i]<-sum(min(q.enr[,i])<0.1)
}
#FDR.sig<-q.enr[,which(KEEP==1)]
my_palette<-colorRampPalette(c('white','orange','red','purple'))
#txt<-ifelse(FDR.sig<0.1,yes = ifelse(FDR.sig<0.05,yes = ifelse(FDR.sig<0.01,yes = '***',no = '**'),no = '*'),no = '')
#heatmap.2(-log10(FDR.sig),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
#          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(FDR)",cellnote = txt)
txt<-ifelse(q.enr<0.05,yes = ifelse(q.enr<0.01,yes = ifelse(q.enr<0.001,yes = '***',no = '**'),no = '*'),no = '')
colnames(q.enr)<-labels2colors(as.numeric(colnames(q.enr))-1)
heatmap.2(-log10(q.enr),trace = 'none',cexCol  = 0.8,cexRow = 0.8,col = my_palette,
          main = 'module preservation',srtCol = 45,notecol = 'black', key.xlab="-log10(Bonferroni corrected p value)",cellnote = txt,colsep=1:ncol(q.enr),rowsep=1:nrow(q.enr),sepcolor="gray",sepwidth=c(0.001,0.001),dendrogram='none',     
          Rowv=FALSE,
          Colv=FALSE)


#__________________________________________________________________________
### CV relation ###
data<-read.delim('~/Desktop/results/repeat.mir.annotated.txt')
#data<-read.delim('~/Desktop/miRNA/qPCR.data.txt',stringsAsFactors=FALSE,colClasses = NA)
miR<-data[,grep('hsa',colnames(data))]
mRNA<-data[,-grep('hsa',colnames(data))]
mean.mRNA<-apply(sapply(mRNA[,-c(1,2)],as.numeric),MARGIN = 2,FUN = mean)
var.mRNA<-apply(sapply(mRNA[,-c(1,2)],as.numeric),MARGIN = 2,FUN = var)
CV.mRNA<-sqrt(var.mRNA)/mean.mRNA

miR.exp<-read.csv('~/Desktop/miRNA/miR_mRNA.csv')
miR.exp$GW15[which(duplicated(miR.exp$ClusterID))]<-0
miR.exp$GW19[which(duplicated(miR.exp$ClusterID))]<-0
miR.exp$GW15<-miR.exp$GW15/(sum(miR.exp$GW15))*1e6
miR.exp$GW19<-miR.exp$GW19/(sum(miR.exp$GW19))*1e6
miR.exp<-cbind(miR.exp[,c(1,2)],(miR.exp$GW15+miR.exp$GW19)*0.5,miR.exp$ClusterID)
expr=seq(length(labels(CV.mRNA)))
for(i in 1:length(labels(CV.mRNA))){
  expr[i]=sum(miR.exp$`(miR.exp$GW15 + miR.exp$GW19) * 0.5`[which(miR.exp$geneID%in%labels(CV.mRNA)[i])])
  #  duplicates<-sum(duplicated(miR.exp$ClusterID[which(miR.exp$geneID%in%labels(CV.mRNA)[i])]))
  #  expr[i]=ifelse(duplicates == 0, expr[i], exprs[i]/duplicates)
}
plot(expr,CV.mRNA,xlab = 'regulation level (CPM) \n (total mRNA-miRNA interactions)',
     ylab = 'target gene CV',pch = 19, axes = FALSE)
lines(expr[which(expr==0)],CV.mRNA[which(expr==0)],type = 'p',col = 'red',pch = 19)
legend(2000,8,c('Spearman cor = -0.46','p value = 1.507e-11'),box.col = 'white')

axis(side = 2, at = c(1,9,18))
axis(side = 1, at = c(0,2300,4600))


#### clustering 

# Randomized PCA
doFastPCA <- function(X,k_num,center.use = T, scale.use = F,iterative = F) {
  
  
  #First center the data
  X.norm <- scale(X, center = center.use,scale = scale.use)
  cen <- attr(X.norm, "scaled:center")
  sc <- attr(X.norm, "scaled:scale")
  
  #Do the Randomized SVD
  svd.results <- doFastSVD(X.norm,k = k_num,iterative)
  
  #Return the results in the same format as prcomp()
  dimnames(svd.results$v) <- list(colnames(X.norm), paste0("PC", seq_len(ncol(svd.results$v))))
  r <- list(sdev = svd.results$d/sqrt(max(1, nrow(X.norm) - 1)), rotation = svd.results$v , center = if (is.null(cen)) FALSE else cen, 
            scale = if (is.null(sc)) FALSE else sc)
  r$x <- X.norm %*% svd.results$v
  class(r) <- "fastprcomp"
  return(r)
}

doFastSVD <- function(X, k, iterative = F) {
  
  
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  Omega <- matrix(rnorm((2 * k) * m), ncol=2 * k)
  
  
  #random projection and orthonormal decomposition
  if (iterative == F) {
    Y <- X %*% Omega
    Q <- qr.Q(qr(Y))
  }
  else if (iterative == T) {
    Q <- RandomizedSubspaceInteration(X,Omega)
  }
  
  #Projection onto this subspace
  B <- t(Q) %*% X
  # decomposing B gives us singular values and right vectors for A  
  s <- svd(B)
  U <- Q %*% s$u
  
  #Now threshold with k  
  s$v <- s$v[,1:k]
  s$d <- s$d[1:k]
  
  #Return in the same format as svd()
  return (list(u=U, v=s$v, d=s$d))
}

RandomizedSubspaceInteration <- function(A, Omega, q = 2) {
  #Finds Q based on an iterative approach
  #Increases accuracy
  
  Y <- A %*% Omega
  Q <- qr.Q(qr(Y))
  for (j in 1:q) {
    Yhat <- t(A) %*% Q
    Qhat <- qr.Q(qr(Yhat))
    Y <- A %*% Qhat
    Q <- qr.Q(qr(Y))
  }
  return(Q)
}

library(Seurat)
library(gmodels)
library(dplyr)
library(Matrix)
library(RANN)
library(igraph)
setwd('~/Desktop/results/Tom/')
#source("doFastPCA.R")
full.data<-read.csv('~/Desktop/results/Tom/read.csv')
data<-full.data[-c(1),-c(1,2,3)]
#colnames(data)<-seq(dim(data)[2]-1)
#load(file="/diazlab/abhaduri/GW18.RData")
GW18<-data;
gw18 <- new("seurat", raw.data = GW18)
gw18 <- Setup(gw18, min.cells = 1, min.genes = 1, do.logNormalize = F, do.scale=F, project = "10X_v2")
#mito.genes <- grep("^MT-", rownames(gw18@data), value = T)
#percent.mito <- colSums(expm1(gw18@data[mito.genes,]))/colSums(expm1(gw18@data))
#gw18 <- AddMetaData(gw18, percent.mito, "percent.mito")
#gw18 <- SubsetData(gw18, subset.name = "percent.mito", accept.high = 0.10)
#ribo.genes <- grep("^RP-", rownames(gw18@data), value = T)
#percent.ribo <- colSums(expm1(gw18@data[ribo.genes,]))/colSums(expm1(gw18@data))
#gw18 <- AddMetaData(gw18, percent.ribo, "percent.ribo")
#gw18 <- SubsetData(gw18, subset.name = "percent.ribo", accept.high = 0.20)
#gw18 <- SubsetData(gw18, subset.name = "nUMI", accept.high = 10000)
matrix <- data
matrix_scaled <- scale(matrix, center=TRUE)
#matrix_pca <- doFastPCA(t(matrix_scaled),50,center.use = F, scale.use = F,iterative = F)
matrix_pca<-fast.prcomp(t(matrix_scaled));
pca.result<-matrix_pca;
write.table(x = matrix_pca$rotation,file = "pca_GW18.txt", sep="\t", quote=FALSE, col.names=NA)
ev <- matrix_pca$sdev^2
u <- max(which(pca.result$sdev^2>(sqrt(length(row.names(data))/length(colnames(data))) + 1)^2))
sig_PCAs <- matrix_pca$rotation[,1:u]
cells_projected_sig_PCAs <- t(matrix_scaled) %*% sig_PCAs
nearest <- nn2(cells_projected_sig_PCAs, k=10)
#save(gw18, cells_projected_sig_PCAs, matrix_pca, ev, matrix_scaled, file="/diazlab/abhaduri/GW18.RData")
rownames(nearest$nn.idx) <- rownames(cells_projected_sig_PCAs)
write.table(nearest$nn.idx, file="nn2_output_neighbors_sigPCAs_gw18.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(nearest$nn.dists, file="nn2_output_distance_sigPCAs_gw18.txt", sep="\t", quote=FALSE, col.names=NA)
system("perl jaccard_toedges2.pl nn2_output_neighbors_sigPCAs_gw18.txt nn2_output_distance_sigPCAs_gw18.txt > jaccard_weighted_edgesgw18.txt")
edgedata <- read.table("jaccard_weighted_edgesgw18.txt", sep="\t", header=T)
edges <- as.data.frame(edgedata)
input <- igraph::graph_from_data_frame(edges, directed=FALSE)
clusters <- igraph::cluster_louvain(input)
clustered <- igraph::membership(clusters)
write(clustered, file="names_incol1gw18.txt", ncol=1)
names <- edgedata[,1]
uniquenames <- unique(names)
meta <- read.table("names_incol1gw18.txt", sep="\t")
rownames(meta) <- uniquenames

cluster <- new("seurat", raw.data = t(cells_projected_sig_PCAs))
cluster <- Setup(cluster, min.cells = 0, min.genes = 0, do.logNormalize = F, project = "10Xdeep")
cluster <- RunTSNE(cluster, dims.use = 1:u, do.fast = T, check_duplicates=FALSE)
gw18 <- AddMetaData(gw18, meta, "groups")
gw18 <- SetAllIdent(gw18, id="V1")
cluster@ident <- gw18@ident
pdf("gw18_tsne.pdf")
TSNEPlot(cluster, do.label=TRUE)
dev.off()
write.table(gw18@ident, "gw18_clusteridentity.txt", sep="\t", quote=F, col.names=NA)
cluster@data <- gw18@data
Tsne<-cluster@tsne.rot
cluster.markers <- FindAllMarkers(cluster, min.pct = 0.25, thresh.use = 0.25, only.pos=TRUE, max.cells.per.ident=10)
write.table(cluster.markers, "gw18_clustermarkers.txt", sep="\t", quote=F, col.names=T)
#save(cluster, gw18, Tsne, cells_projected_sig_PCAs, matrix_pca, ev, matrix_scaled, file="/diazlab/abhaduri/GW18.RData")
write.csv(cluster@tsne.rot, file = paste("gw18_tSNECoordinates.csv",sep=""))



b='HES1';
tiff(paste(b,".tiff"));print(qplot(Tsne[,1],Tsne[,2],color=data[,a],main=b) + scale_colour_gradient(low="gray", high="blue", space="Lab")+theme_bw(base_size = 10));dev.off();


for(i in 1:length(colnames(data))){pdf(paste(colnames(data)[i],".pdf"));print(qplot(Tsne[,2],Tsne[,3],color=data[,colnames(data)[i]]) + scale_colour_gradient(low="gray", high="blue", space="Lab")+theme_bw(base_size = 10));dev.off();}






