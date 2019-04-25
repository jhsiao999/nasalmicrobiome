library(matrixStats)
library(dplyr)
library(ggplot2)
library(vegan)
library(rafalib)
library(reshape2)
library(devtools)
library(metagenomeSeq)
obj = readRDS("data/nasal_GOM.rds")

# obj=MRobj[,order(colnames(MRobj))]
# lvl='Genus'
# n=20
# plot.factor_levels = labs

plotBar.aggregate <- function(obj, lvl, plot.factor_levels=NULL,
                              cl=colnames(obj), n=20, norm=FALSE, log=FALSE, ord=FALSE,
                              orderby='Other',...){
  if (class(obj) == "MRexperiment") {
    mat = MRcounts(obj, norm = norm, log = log)
    if(length(lvl)==1){
      lvl = fData(obj)[,lvl]
    }
  }
  prop <- cbind(rowSums(mat)/sum(mat),rowSums(mat)/sum(mat))
  aggProp = aggregateByTaxonomy(prop,lvl,out='matrix')
  ordIndex  = order(rowSums(aggProp),decreasing=TRUE)[1:n]
  Other = 1-colSums(aggProp[ordIndex,])
  aggPropSub = cbind(t(aggProp[ordIndex,]),Other)

  dd <- data.frame(all=aggPropSub[1,], variable=colnames(aggPropSub))

  var_factor <- factor(dd$variable, levels=dd$variable[order(dd$all,
                                                             decreasing = T)])

  dd$variable <- var_factor

  if (!is.null(plot.factor_levels)) {
  levels(dd$variable) <- plot.factor_levels
  }
  p <- ggplot(dd,aes(1, all, fill=variable)) +
          geom_bar(position="fill", stat="identity") +
      #    theme(legend.position = "none") +
          ylab("Proportion") +
          labs(fill='Taxa')
  return(list(p=p,dd=dd))
}

#lvl1 <- levels(dd$variable)


# obj = readRDS("data/nasal_GOM.rds")
# obj <- obj[,order(colnames(obj))]
# plot.factor_levels <- lvl1

plotBar <- function(obj, lvl, plot.factor_levels=NULL,
                    cl=colnames(obj), n=10, norm=FALSE, log=FALSE, ord=FALSE,
                    orderby='Other',...){
  if (class(obj) == "MRexperiment") {
    mat = MRcounts(obj, norm = norm, log = log)
    if(length(lvl)==1){
      lvl = fData(obj)[,lvl]
    }
  }
  prop = prop.table(mat,2)
  aggProp = aggregateByTaxonomy(prop,lvl,out='matrix')
  ordIndex  = order(rowSums(aggProp),decreasing=TRUE)[1:n]
  Other = 1-colSums(aggProp[ordIndex,])
  aggPropSub = cbind(t(aggProp[ordIndex,]),Other)
  if(length(unique(cl))!=nrow(aggPropSub)){
    cl = factor(cl)
    aggPropSub = by(aggPropSub,cl,colMeans)
    aggPropSub = Reduce("rbind",aggPropSub)
    rownames(aggPropSub) = levels(cl)
    cl = levels(cl)
  }
  propDF = data.frame(Group=cl,aggPropSub)
  if(ord){
    rord = order(aggPropSub[,orderby])
    propDF=propDF[rord,]
    cl = reorder(cl,rord); propDF$Group = cl # ???
  }
  dd = melt(propDF,id.vars=c("Group"),measure.vars=colnames(propDF)[-1])
  labs <- as.character(dd$variable)
  if (is.null(plot.factor_levels)) {
    dd$variable <- factor(dd$variable)
  } else {
    var_refactor <- factor(labs, levels=plot.factor_levels)
    dd$variable <- var_refactor
  }
  p=ggplot(dd,aes(Group,value,fill=variable)) +
    geom_bar(position="fill", stat='identity') +
    ylab("Proportion") +
    scale_x_discrete(labels = cl,limits=cl)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6)) +
    labs(fill='Taxa')
  return(list(p=p,dd=dd))
}


obj = readRDS("data/nasal_GOM.rds")

dd_agg <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$dd
lvl1 <- levels(dd_agg$variable)
#labs <- as.character(pp$dd$variable)
# bb <- unique(as.character(dd_agg$variable))
# bb <- bb[order(bb)]
# ll <- unique(labs)
# ll <- ll[order(ll)]


#head(cbind(bb,ll))
obj = readRDS("data/nasal_GOM.rds")
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

b <- plotBar(obj[,order(colnames(obj))],lvl='Genus',n=20,
        plot.factor_levels = levels(dd_agg$variable))$p +
  scale_fill_manual(values = getPalette(21)) + xlab('Samples') +
  theme(legend.position = "none")



plotOrd2 = function(obj,tran=TRUE,comp=1:2,norm=TRUE,log=TRUE,usePCA=TRUE,useDist=FALSE,distfun=stats::dist,dist.method="euclidian",n=NULL,...){
  mat = returnAppropriateObj(obj,norm,log)
  if(useDist==FALSE & usePCA==FALSE) stop("Classical MDS requires distances")
  if(is.null(n)) n = min(nrow(mat),1000)

  otusToKeep <- which(rowSums(mat)>0)
  otuVars<-rowSds(mat[otusToKeep,])
  otuIndices<-otusToKeep[order(otuVars,decreasing=TRUE)[seq_len(n)]]
  mat <- mat[otuIndices,]

  if(tran==TRUE){
    mat = t(mat)
  }
  if(useDist==TRUE){
    d <- distfun(mat,method=dist.method)
  } else{ d = mat }

  if(usePCA==FALSE){
    ord = cmdscale(d,k = max(comp))
    xl = paste("MDS component:",comp[1])
    yl = paste("MDS component:",comp[2])
  } else{
    pcaRes <- prcomp(d)
    ord <- pcaRes$x
    vars <- pcaRes$sdev^2
    vars <- round(vars/sum(vars),5)*100
    xl <- sprintf("PCA %s: %.2f%% variance",colnames(ord)[comp[1]], vars[comp[1]])
    yl <- sprintf("PCA %s: %.2f%% variance",colnames(ord)[comp[2]], vars[comp[2]])
  }

  plot(ord[,comp],ylab=yl,xlab=xl,...)
  invisible(ord[,comp])
}


pcaRES <- function(obj,tran=TRUE,comp=1:2,norm=TRUE,log=TRUE,usePCA=TRUE,useDist=FALSE,distfun=stats::dist,dist.method="euclidian",n=NULL,...){
  mat = returnAppropriateObj(obj,norm,log)
  if(useDist==FALSE & usePCA==FALSE) stop("Classical MDS requires distances")
  if(is.null(n)) n = min(nrow(mat),1000)

  otusToKeep <- which(rowSums(mat)>0)
  otuVars<-rowSds(mat[otusToKeep,])
  otuIndices<-otusToKeep[order(otuVars,decreasing=TRUE)[seq_len(n)]]
  mat <- mat[otuIndices,]

  if(tran==TRUE){
    mat = t(mat)
  }
  if(useDist==TRUE){
    d <- distfun(mat,method=dist.method)
  } else{ d = mat }

  if(usePCA==FALSE){
    ord = cmdscale(d,k = max(comp))
    xl = paste("MDS component:",comp[1])
    yl = paste("MDS component:",comp[2])
  } else{
    pcaRes <- prcomp(d)
    ord <- pcaRes$x
    vars <- pcaRes$sdev^2
    vars <- round(vars/sum(vars),5)*100
    xl <- sprintf("PCA %s: %.2f%% variance",colnames(ord)[comp[1]], vars[comp[1]])
    yl <- sprintf("PCA %s: %.2f%% variance",colnames(ord)[comp[2]], vars[comp[2]])
  }
  return(pcaRes)
}


# taxa proportions across samples -----------------------------------------------------
obj = readRDS("data/nasal_GOM.rds")

#obj <- obj[,order(colnames(obj))]; lvl="Genus"; n=20

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/proportion_sample_allaggregate_nolegend.pdf", width=4,height=6)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  cols <- getPalette(21)
  plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$p +
    scale_fill_manual(values = cols) + xlab('Samples') +
    theme(legend.position = "none")
dev.off()

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/proportion_sample_allaggregate.pdf",
    width=12,height=6)
plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$p +
  scale_fill_manual(values = getPalette(21)) + xlab('Samples')
dev.off()


# get table
obj = readRDS("data/nasal_GOM.rds")
obj <- obj[,order(pData(obj)$GOM)]

obj <- obj[,order(colnames(obj))]
lvl <- 'Genus'
n <- 20
norm=FALSE
log=FALSE
ord=FALSE

mat = MRcounts(obj, norm = norm, log = log)
lvl = fData(obj)[,lvl]
prop <- cbind(rowSums(mat)/sum(mat),rowSums(mat)/sum(mat))
aggProp = aggregateByTaxonomy(prop,lvl,out='matrix')
ordIndex  = order(rowSums(aggProp),decreasing=TRUE)[1:n]
Other = 1-colSums(aggProp[ordIndex,])
aggPropSub = cbind(t(aggProp[ordIndex,]),Other)
dd <- data.frame(all=aggPropSub[1,], variable=colnames(aggPropSub))

write.table(dd,
            file= "output-manuscript/proportion_sample_allaggregate.csv",
            sep="\t", quote=F, row.names=F, col.names=F)



# taxa proportions at the sample level ----------------------------------------------------
obj = readRDS("data/nasal_GOM.rds")
dd_agg <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$dd
labs <- as.character(pp$dd$variable)
var_factor <- levels(dd_agg$variable)


pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/proportion_sample.pdf",
    width=15,height=7)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  plotBar(obj[,order(colnames(obj))],lvl='Genus',n=20, plot.factor_levels = var_factor)$p +
    scale_fill_manual(values = getPalette(21)) + xlab('Samples')
dev.off()

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/proportion_sample_nolegend.pdf",
    width=15,height=5)
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))

  plotBar(obj[,order(colnames(obj))],lvl='Genus',n=20, plot.factor_levels = var_factor)$p +
    scale_fill_manual(values = getPalette(21)) + xlab('Samples') +
    theme(legend.position = "none")
dev.off()


##------ sample proportions by GOM membership
pdf("proportion_sample_gom_nolegend.pdf",width=12,height=8)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plotBar(obj,lvl='Genus',cl=pData(obj)$GOM,n=20)$p +
  scale_fill_manual(values = getPalette(21))
#dev.off()
plotBar(obj,lvl='Genus',cl=pData(obj)$GOM,n=10)$p +
  scale_fill_manual(values = getPalette(11))
dev.off()

pdf("proportion_sample_gom.pdf",width=8,height=8)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plotBar(obj,lvl='Genus',cl=pData(obj)$GOM,n=20)$p +
  scale_fill_manual(values = getPalette(21)) +
  theme(legend.position = "none")

# plotBar(obj,lvl='Genus',cl=pData(obj)$GOM,n=10)$p +
#   scale_fill_manual(values = getPalette(11)) +
#   theme(legend.position = "none")
dev.off()

# DA-Testing
mod = model.matrix(~pData(obj)$GOM)
fit = fitFeatureModel(obj,mod)
res = MRfulltable(fit,number = nrow(obj))
res = cbind(fData(obj[rownames(res),])[,-7],res)
write.table(res,file="~/Desktop/gom_comparisons2.tsv",sep="\t",quote=FALSE)

###---PCA
pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/pca.pdf",width=6,heigh=6)
  pcares <- pcaRES(obj,pch=21,bg=pData(obj)$GOM,main="")
  varprop <- (pcares$sdev^2)/sum((pcares$sdev^2))
#  plotOrd(obj,pch=21,bg=pData(obj)$GOM,main="Principal Component Analysis")
  k = plotOrd2(obj,pch=21,bg=pData(obj)$GOM,comp = 1:4,
               main="Principal Component Analysis")
  pairs(k,pch=21,bg=pData(obj)$GOM,main="Principal Component Analysis",
        labels = c("PC1: 9.45%", "PC2: 8.62%", "PC3: 4.40%", "PC4: 4.28%"))
dev.off()

# check outliers on PC2
pc_outlier_indices <- pcares$x[,1] > 20 & pcares$x[,2] < -20
pc_outlier_samples <- rownames(pcares$x)[which(pcares$x[,1] > 20 & pcares$x[,2] < -20)]
all.equal(colnames(obj), rownames(pcares$x))
table(pc_outlier_indices, pData(obj)$anyinf30)
table(pc_outlier_indices, pData(obj)$anyinf6m)

fits <- readRDS(file = "data/gomfits.rds")
fits_clust2_prob <- fits[[1]]$omega[match(rownames(pcares$x),
                                                    rownames(fits[[1]]$omega)),]
all.equal(rownames(fits_clust2_prob), rownames(pcares$x))
fits_clust2_prob[pc_outlier_indices,]
colSums(obj)[pc_outlier_indices]


# check outliers on PC2
pc_clust2_idx <- colnames(obj)[pData(obj)$GOM == 2]
fits <- readRDS(file = "data/gomfits.rds")
fits_clust2_prob <- fits[[1]]$omega[match(rownames(pcares$x),
                                          rownames(fits[[1]]$omega)),]
all.equal(rownames(fits_clust2_prob), rownames(pcares$x))
fits_clust2_prob[which(rownames(fits_clust2_prob) %in% pc_clust2_idx),]




pc_clust2_idx <- pc_clust2_idx[!(pc_clust2_idx %in% pc_outlier_samples)]
# pc_clust2_idx_samples <- rownames(pcares$x)[which(pc_clust2_idx)]
# all.equal(colnames(obj), rownames(pcares$x))
# table(pc_outlier_indices, pData(obj)$anyinf30)
# table(pc_outlier_indices, pData(obj)$anyinf6m)

#colSums(obj)[pc_outlier_indices]


###--- beta-Diversity
pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/beta_diversity.pdf",
    width=4,heigh=3)
pcares <- pcaRES(obj,pch=21,bg=pData(obj)$GOM,main="")
H = diversity(t(MRcounts(obj)))
gom = pData(obj)$GOM
levels(gom) <- c("Group 1", "Group 2")
library(ggplot2)
ggplot(data.frame(H=H,gom=gom), aes(x=gom, y = H, fill = gom)) +
  geom_violin(col = "gray50") +
  scale_fill_manual(values = c(brewer.pal(12,"Set3")[1],
                               brewer.pal(12,"Set3")[2])) +
  geom_boxplot(width=.25) + ylab("Shannon diversity") +
  scale_x_discrete(name = "Cluster membership",
                   labels = c("Group 1", "Group 2")) +
  guides(fill=guide_legend(title="")) +
  ggtitle("Alpha diversity")
dev.off()


###--- alpha-Diversity
pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/shannon.pdf",width=4,height=3)
  H = diversity(t(MRcounts(obj)))
  gom = pData(obj)$GOM
  levels(gom) <- c("Group 1", "Group 2")
  library(ggplot2)
  ggplot(data.frame(H=H,gom=gom), aes(x=gom, y = H, fill = gom)) +
    geom_violin(col = "gray50") +
    scale_fill_manual(values = c(brewer.pal(12,"Set3")[1],
                                   brewer.pal(12,"Set3")[2])) +
    geom_boxplot(width=.25) + ylab("Shannon diversity") +
    scale_x_discrete(name = "Cluster membership",
                     labels = c("Group 1", "Group 2")) +
    guides(fill=guide_legend(title="")) +
    ggtitle("Alpha diversity")
dev.off()
