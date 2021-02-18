# Date: 02-17-2020
#
# Document and reproduce figures/data for the paper
#

# load packages -----------------------------------------
library(matrixStats)
library(dplyr)
library(ggplot2)
library(vegan)
library(reshape2)
library(devtools)
library(metagenomeSeq)
library(tidyr)


# data -----------------------------------------------------
# batch-correctec counts with grades of membership labels
obj = readRDS("data/nasal_GOM.rds")


# Figure 1A ------------------------------------------------
# aggregated count bar chart

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


pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/fig1a_proportion_sample_allaggregate.pdf",
    width=5,height=6)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(21)
fig1_a <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)
labs = as.vector(sapply(levels(fig1_a$p$data$variable), function(x) {
  #strsplit(x, split = "_")[[1]]
  tmp = strsplit(x, split = c("_"))[[1]]
  if (length(tmp) == 1) out = "Other"
  if (length(tmp) == 3) out = tmp[3]
  if (length(tmp) > 3) {
    tmp2 = paste0("Uncl_", tmp[length(tmp)-1])
    out = strsplit(tmp2, split = ";")[[1]][1]
  }
  return(out)
  }))

fig1_a$p$data$variable = factor(fig1_a$p$data$variable,
                                levels = levels(fig1_a$p$data$variable)[c(1,3:21,2)],
                                labels = labs[c(1,3:21,2)])

print(
  fig1_a$p + scale_fill_manual(values = cols) + xlab('Samples') +
  ylab("Relative Abundance") +
  xlab("All Samples") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(#axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14),
        legend.key.size = unit(.5, 'cm'),
        legend.text = element_text(size=6)) +
  guides(fill = guide_legend(title = "Taxa",
                             title.position = "top", ncol=1)) +
    ggtitle("Figure 1A. Relative abundance of top 20 taxa \n to the genus level (aggregate data")
)
dev.off()




# Figure 1A corresponding data table ----------------
obj = readRDS("data/nasal_GOM.rds")
obj <- obj[,order(pData(obj)$GOM)]

# obj <- obj[,order(colnames(obj))]
# lvl <- 'Genus'
# n <- 20
# norm=FALSE
# log=FALSE
# ord=FALSE

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



# Figure 1B --------------------------------
# taxa proportions at the sample level
# obj=obj[,order(colnames(obj))]
# lvl='Genus'
# n=20
# plot.factor_levels = var_factor
#pp <- plotBar(obj[,order(colnames(obj))],lvl='Genus',cl=membership$clust_2,n=20)


plotBar <- function(obj, lvl, plot.factor_levels=NULL,
                    cl=colnames(obj), n=10, norm=FALSE, log=FALSE, ord=FALSE,
                    orderby='Other',...){
  if (class(obj) == "MRexperiment") {
    mat = MRcounts(obj, norm = norm, log = log)
    if(length(lvl)==1){
      lvl = fData(obj)[,lvl]
    }
  }
  # get top 20 order
  prop <- cbind(rowSums(mat)/sum(mat),rowSums(mat)/sum(mat))
  aggProp = aggregateByTaxonomy(prop,lvl,out='matrix')
  ordIndex  = order(rowSums(aggProp),decreasing=TRUE)[1:n]
  features_to_plot = rownames(aggProp)[ordIndex]

  prop = prop.table(mat,2)
  aggProp = aggregateByTaxonomy(prop,lvl,out='matrix')
  #ordIndex  = order(rowSums(aggProp),decreasing=TRUE)[1:n]
  aggProp_contain_features <- aggProp[which(rownames(aggProp)%in%features_to_plot), ]
  aggProp_contain_features <- aggProp_contain_features[match(features_to_plot,
                                                       rownames(aggProp_contain_features)),]
  Other = 1-colSums(aggProp_contain_features)
  aggPropSub = cbind(t(aggProp[which(rownames(aggProp)%in%features_to_plot), ]),Other)
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
  dd <- tidyr::gather(propDF, variable, value, colnames(propDF)[2]:colnames(propDF)[22],
               factor_key=F)
  dd$variable <- gsub(".", ";", dd$variable, fixed=T)
#  dd = melt(propDF,id.vars=c("Group"),measure.vars=colnames(propDF)[-1])
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
fig_a_dd <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$dd
labs <- as.character(fig_a_dd$variable)
var_factor <- levels(fig_a_dd$variable)

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/fig1b_sample_proportion_individual.pdf",
    width=12,height=7)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(21)
plotBar(obj[,order(colnames(obj))],lvl='Genus',n=20, plot.factor_levels = var_factor)$p +
  scale_fill_manual(values = getPalette(21)) + xlab('Samples') +
  theme_classic() + coord_fixed(ratio = 60) + theme(legend.position = "none") +
  ylab("Relative Abundance") + xlab("Individual Samples") +
  theme(axis.text.x = element_text(size = 4, angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14)) +
  ggtitle("Figure 1B. Relative abundance of to 20 taxa to the genus level (aggregate data)")
#  guides(fill = guide_legend(title = "", title.position = "left", ncol=1))
dev.off()





#--- Figure 2A
# Shannon diversity between groups
library(data.table)
shannon_values <- fread("~/Dropbox/current-projects/nasalMicrobiomeManuscript/draft-results/shannon_values.tsv")
obj <- readRDS("data/nasal_GOM.rds")

all.equal(shannon_values$V1, colnames(obj))
shannon_values$cluster <- obj$GOM
wilcox.test(shannon_values$V2 ~ shannon_values$cluster)

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/fig2a_alpha_diversity.pdf",
    width=4,heigh=3)
#pcares <- pcaRES(obj,pch=21,bg=pData(obj)$GOM,main="")
H = diversity(t(MRcounts(obj)))
gom = pData(obj)$GOM
levels(gom) <- c("Group 1", "Group 2")
ggplot(data.frame(H=H,gom=gom), aes(x=gom, y = H, fill = gom)) +
  geom_violin(col = "gray50") +
  scale_fill_manual(values = c(brewer.pal(12,"Set3")[1],
                               brewer.pal(12,"Set3")[2]),
                    labels = c("Group 1", "Group 2")) +
  geom_boxplot(width=.25) + ylab("Shannon diversity") +
  scale_x_discrete(name = "Cluster membership",
                   labels = c("Group 1", "Group 2")) +
  guides(fill=guide_legend(title=""), title.position = "top") +
  ggtitle("Figure 2A. Microbiome diversity metrics \n by cluster class: Alpha (Shannon) diversity") +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  geom_text(
    label="P <.001",
    x=1.5,
    y=3.6,
#    label.padding = unit(0.55, "lines"), # Rectangle size around label
    size = 3
  )
dev.off()





#--- Figure 2B
# Bray-Curtis similarity
MRobj <- readRDS("data/nasal_filtered_normed_batchcorrected.rds")
counts <- MRcounts(MRobj,norm=T,log=T)
colnames(counts) <- colnames(MRobj)

obj <- readRDS("data/nasal_GOM.rds")
pdata <- pData(obj)

all.equal(rownames(pdata),colnames(counts))

library(vegan)
beta_bray <- vegdist(t(counts), method = "bray", upper=T)
beta_bray <- as.matrix(beta_bray)

group1 <- rownames(pdata)[which(pdata$GOM==1)]
group2 <- rownames(pdata)[which(pdata$GOM==2)]

beta_group1 <- lapply(1:length(group1), function(n) {
  beta_bray[which(rownames(beta_bray)==group1[n]),
            which(colnames(beta_bray) %in% group1)]
})
beta_group2 <- lapply(1:length(group2), function(n) {
  beta_bray[which(rownames(beta_bray)==group2[n]),
            which(colnames(beta_bray) %in% group2)]
})


wilcox.test(df$beta~df$clust)

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/fig2b_beta_diversity_mean.pdf",
    width=4,heigh=3)
beta_group1_mn <- sapply(1:length(group1), function(n) {
  mean(beta_bray[which(rownames(beta_bray)==group1[n]),
                 which(colnames(beta_bray) %in% group1)])
})
beta_group2_mn <- sapply(1:length(group2), function(n) {
  mean(beta_bray[which(rownames(beta_bray)==group2[n]),
                 which(colnames(beta_bray) %in% group2)])
})

df <- data.frame(clust=rep(c(1,2), times=c(167,30)),
                 beta=c(beta_group1_mn,beta_group2_mn),
                 row.names = c(group1, group2))
df$clust <- factor(df$clust)
levels(df$clust) <- c("Group 1", "Group 2")
library(ggplot2)
ggplot(df, aes(x=clust, y = beta, fill = clust)) +
  geom_violin(col = "gray50") +
  scale_fill_manual(values = c(brewer.pal(12,"Set3")[1],
                               brewer.pal(12,"Set3")[2])) +
  geom_boxplot(width=.25) + ylab("Mean Bray-Curtis across samples") +
  scale_x_discrete(name = "Cluster membership",
                   labels = c("Group 1", "Group 2")) +
  guides(fill=guide_legend(title="")) +
  ggtitle("Figure 2B. Microbiome diversity metrics \n by cluster class: Beta (Bray-Curtis) diversity") +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8)) +
  geom_text(
    label="P <1E-15",
    x=1.5,
    y=3.6,
    #    label.padding = unit(0.55, "lines"), # Rectangle size around label
    size = 3
  )
dev.off()




#---Figure 2C: PCA
pcaRES <- function(obj,tran=TRUE,comp=1:2,norm=TRUE,log=TRUE,usePCA=TRUE,
                   useDist=FALSE,distfun=stats::dist,dist.method="euclidian",n=NULL,...){
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


pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/fig2c_pca.pdf",
    width=8,height=8)
pcares <- pcaRES(obj,pch=21,bg=pData(obj)$GOM,main="")
varprop <- (pcares$sdev^2)/sum((pcares$sdev^2))
#  plotOrd(obj,pch=21,bg=pData(obj)$GOM,main="Principal Component Analysis")
k = plotOrd2(obj,pch=21,bg=pData(obj)$GOM,comp = 1:4,
             main="Principal Component Analysis")
pairs(k,pch=21,bg=pData(obj)$GOM,
      main="",
      labels = c("PC1: 9.45%", "PC2: 8.62%", "PC3: 4.40%", "PC4: 4.28%"))
mtext(side=3, line=1.7, at=0.05, adj=0, cex=1.3, font = 2,
      "Figure 2C. Microbiome diversity metrics \n by cluster class: Principal Component Analysis")
dev.off()







# Supp? Figure 2 --------------------------------
# sample proportions by GOM membership
pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/proportion_sample_gom.pdf",
    width=12,height=8)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
plotBar(obj,lvl='Genus',cl=pData(obj)$GOM,n=20, plot.factor_levels = var_factor)$p +
  scale_fill_manual(values = getPalette(21)) +
  guides(fill = guide_legend(title = "", title.position = "left", ncol=1)) +
  xlab("Cluster Membership") + theme_light()
dev.off()
