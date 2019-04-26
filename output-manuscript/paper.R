# Date: 04-26-2019
#
# Author: Joyce Hsiao
#
# Document and reproduce figures/data for the paper
#

# load packages -----------------------------------------
library(matrixStats)
library(dplyr)
library(ggplot2)
library(vegan)
#library(rafalib)
library(reshape2)
library(devtools)
library(metagenomeSeq)
library(tidyr)


# data -----------------------------------------------------
# batch-correctec counts with grades of membership labels
obj = readRDS("data/nasal_GOM.rds")


# Figure 1A ------------------------------------------------
# aggregated count bar chart
obj=obj[,order(colnames(obj))]
lvl='Genus'
n=20

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
    dd$variable <- factor(dd$variable, levels = plot.factor_levels)
  }
  p <- ggplot(dd,aes(1, all, fill=variable)) +
    geom_bar(position="fill", stat="identity") +
    #theme(legend.position = "none") +
    ylab("Proportion") +
    labs(fill='Taxa')
  return(list(p=p,dd=dd))
}


pdf("output-manuscript/fig1a_sample_proportion_aggregated.pdf",
    width=12,height=6)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(21)
fig1_a <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)
write.table(fig1_a$dd,
            file= "output-manuscript/fig1a_sample_proportion_aggregated.csv",
            sep="\t", quote=F, row.names=F, col.names=F)

var_levels <- levels(fig1_a$dd$variable)
var_levels <- c(levels(fig1_a$dd$variable)[-2], "Other")
print(
  plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20,
                    plot.factor_levels=var_levels)$p +
    scale_fill_manual(values = cols) + xlab('Samples') +
    ylab("Genus Proportion") +
    xlab("All Samples") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    guides(fill = guide_legend(title = "", title.position = "left", ncol=1))
)
dev.off()


# Figure 1B --------------------------------
# taxa proportions at the sample level
# obj=obj[,order(colnames(obj))]
# lvl='Genus'
# n=20
# plot.factor_levels = var_factor

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
fig_1_dd <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$dd
#labs <- levels(fig_1_dd$variable)
var_factor <- c(levels(fig_1_dd$variable)[-2], "Other")

pdf("output-manuscript/fig1b_sample_proportion_individual.pdf",
    width=15,height=7)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(21)
fig1_b <- plotBar(obj[,order(colnames(obj))],lvl='Genus',n=20,
                  plot.factor_levels = var_factor)

write.table(fig1_b$dd,
            file= "output-manuscript/fig1b_sample_proportion_individual.csv",
            sep="\t", quote=F, row.names=F, col.names=F)

fig1_b$p +
  scale_fill_manual(values = getPalette(21)) + xlab('Samples') +
  theme(legend.position = "none")
#  guides(fill = guide_legend(title = "", title.position = "left", ncol=1))
dev.off()



# Figure 2 --------------------------------
# sample proportions by GOM membership
obj = readRDS("data/nasal_GOM.rds")
fig_1_dd <- plotBar.aggregate(obj[,order(colnames(obj))],lvl='Genus',n=20)$dd
#labs <- levels(fig_1_dd$variable)
var_factor <- c(levels(fig_1_dd$variable)[-2], "Other")

pdf("output-manuscript/fig2_sample_proportion_clusters.pdf",
    width=12,height=8)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fig2 <- plotBar(obj,lvl='Genus',cl=pData(obj)$GOM,n=20, plot.factor_levels = var_factor)

write.table(fig2$dd,
            file= "output-manuscript/fig2_sample_proportion_clusters.csv",
            sep="\t", quote=F, row.names=F, col.names=F)

fig2$p +
  scale_fill_manual(values = getPalette(21)) +
  guides(fill = guide_legend(title = "", title.position = "left", ncol=1)) +
  xlab("Cluster Membership") + theme_light()

dev.off()




# Figure 3A ------------------------------
# Shannon diversity between groups
library(data.table)
shannon_values <- fread("~/Dropbox/current-projects/nasalMicrobiomeManuscript/draft-results/shannon_values.tsv")
obj <- readRDS("data/nasal_GOM.rds")

all.equal(shannon_values$V1, colnames(obj))
shannon_values$cluster <- obj$GOM
wilcox.test(shannon_values$V2 ~ shannon_values$cluster)

write.table(shannon_values,
            file= "output-manuscript/fig3a_alpha_diversity.csv",
            sep="\t", quote=F, row.names=F, col.names=F)

pdf("output-manuscript/fig3a_alpha_diversity.pdf",
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
  ggtitle("Alpha (Shannon) diversity")
dev.off()





# Figure 3B ------
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



pdf("output-manuscript/fig3b_beta_diversity.pdf",
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

write.table(df,
            file= "output-manuscript/fig3b_beta_diversity.csv",
            sep="\t", quote=F, row.names=F, col.names=F)

library(ggplot2)
ggplot(df, aes(x=clust, y = beta, fill = clust)) +
  geom_violin(col = "gray50") +
  scale_fill_manual(values = c(brewer.pal(12,"Set3")[1],
                               brewer.pal(12,"Set3")[2])) +
  geom_boxplot(width=.25) + ylab("Mean Bray-Curtis across samples") +
  scale_x_discrete(name = "Cluster membership",
                   labels = c("Group 1", "Group 2")) +
  guides(fill=guide_legend(title="")) +
  ggtitle("Beta (Bray-Curtis) diversity")
wilcox.test(df$beta~df$clust)
dev.off()
