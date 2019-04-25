# make figures/tables for the manuscript


###--- Compute variance explained under different Ks
# see the link here for original analysis
# https://jhsiao999.github.io/nasalmicrobiome/cluster-variance-explained.html

MRobj <- readRDS("data/nasal_filtered_normed_batchcorrected.rds")
counts <- MRcounts(MRobj,norm=FALSE,log=FALSE)
colnames(counts) <- colnames(MRobj)

fits <- readRDS(file = "data/gomfits.rds")

membership <- lapply(fits, function(xx) {
  apply(xx$omega, 1, which.max)
})

varprop_species <- lapply(1:length(membership), function(index) {
  which_sample <- which(colnames(counts) %in% names(membership[[index]]))
  counts_which <- counts[, which_sample]
  var <- sapply(1:nrow(counts_which), function(index2) {
    res <- lm(log2(counts_which[index2,]+2) ~ membership[[index]])
    sumOfsquares <- anova(res)[[2]]
    sumOfsquares[1]/sum(sumOfsquares) })
  return(var)
})
names(varprop_species) <- names(membership)
saveRDS(varprop_species, file = "data/varprop.rds")


varprop <- readRDS("data/varprop.rds")
# take results up to 5 clusters
varprop <- varprop[2:5]
varprop_long <- data.frame(
  prop = do.call(c, varprop),
  cluster = as.factor(rep(2:5, each = length(varprop[[1]]))))

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/variance-explained.pdf",
    width=6,height=4)
library(ggplot2)
ggplot(varprop_long, aes(x = cluster, y = prop,
                         col = cluster), size = .2) +
  geom_violin(trim = FALSE, size=.5) +
  geom_boxplot(width=0.1, size = .6) +
  labs(x = "Number of clusters", y = "Proportion variance explained")
dev.off()



###--- make flow diagram
raw <- readRDS("data/nasal_raw.rds")
dim(raw)

filt <- readRDS("data/nasal_filtered.rds")
dim(filt)

batch <- readRDS("data/nasal_filtered_normed_batchcorrected.rds")
dim(batch)

batch <- readRDS("data/nasal_filtered_normed_batchcorrected.rds")
dim(batch)

###--- make a table for all taxa identified at species level
obj <- readRDS("data/nasal_GOM.rds")
prop <- sweep(MRcounts(obj), 2, STAT=colSums(obj), "/")
prop <- prop[order(rowSums(prop)),]
write.csv(prop, file = "output-manuscript/taxa-proportions.csv",
            row.names = T, quote=F)


# make heatmap for 66 features with FDR < .05 and abs(logFC) > 1
library(data.table)
gom_comparisons <- fread("output-manuscript/gom_comparisons.tsv")
obj <- readRDS("data/nasal_GOM.rds")

# match gom_comparisons with obj
ord <- match(fData(obj)$Species, gom_comparisons$Species)
gom_comparisons <- gom_comparisons[ord,]
all.equal(gom_comparisons$Species, fData(obj)$Species)

species <- gom_comparisons$Species
fdr <- gom_comparisons$adjPvalues
logFC <- gom_comparisons$logFC
tsamples <- rowSums(gom_comparisons[,2:3])
treads <- rowSums(gom_comparisons[,4:5])

obj_sig <- obj[fdr < .05 & abs(logFC) > 1,]
obj_sig <- obj_sig[,order(pData(obj_sig)$GOM)]
fdr_sig <- fdr[fdr < .05 & abs(logFC) > 1]

# phylum <- sapply(strsplit(rownames(obj_sig),';c__'), function(i) gsub("k__Bacteria;p__","",i)[1])
# genus <- gsub("g__","",sapply(strsplit(rownames(obj_sig),';'), function(i) i[6]))
# genus2 <- genus
# genus2[genus=="" | is.na(genus)] <- paste("Genus", c(1:sum(genus=="" | is.na(genus)) ))

# function to plot heatmap
#plotMRheatmap2<-function(obj, n, norm = TRUE, log = TRUE, fun = sd, ...){
obj = obj_sig; n=nrow(obj_sig); norm=TRUE; log=TRUE  #; fun=sd
mat = returnAppropriateObj(obj, norm, log)
mat2 <- mat[order(fdr_sig),]

ann_col <- data.frame(Group=factor(paste0("Group.",
                          as.integer(factor(pData(obj_sig)$GOM)))),
                        row.names = colnames(mat2))
ann_colors = list(
  Cluster=c(Group.1=brewer.pal(12,"Set3")[1],
            Group.2=brewer.pal(12,"Set3")[2]))
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)

pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/heatmap-fdr05-logfc1.pdf",
    width=10,height=10)
  library(pheatmap)
  pheatmap(mat2,
           main="OTU Abundance Heatmap",
           fontsize_row = 3,
           fontsize_col = 3,
           border_color = NA,
           annotation_col = ann_col,
           cluster_rows = F, cluster_cols = F,
           col = heatmapCols,
           annotation_colors = ann_colors)
dev.off()



pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/heatmap-fdr05-logfc1-top4taxa.pdf",
    width=8,height=4)
library(data.table)
library(metagenomeSeq)
gom_comparisons <- fread("output-manuscript/gom_comparisons.tsv")
obj <- readRDS("data/nasal_GOM.rds")

# match gom_comparisons with obj
ord <- match(fData(obj)$Species, gom_comparisons$Species)
gom_comparisons <- gom_comparisons[ord,]
all.equal(gom_comparisons$Species, fData(obj)$Species)

species <- gom_comparisons$Species
fdr <- gom_comparisons$adjPvalues
logFC <- gom_comparisons$logFC
tsamples <- rowSums(gom_comparisons[,2:3])
treads <- rowSums(gom_comparisons[,4:5])

norm=TRUE; log=TRUE  #; fun=sd
mat = returnAppropriateObj(obj, norm, log)

library(data.table)
ref <- fread("~/Dropbox/current-projects/nasalMicrobiomeManuscript/draft-results/individual_level_results_20180513.csv")
top4 <- ref$top_taxa[ref$adj_pvalue<.05]
mat_plot <- mat[which(rownames(mat) %in% top4),]
mat_plot <- mat_plot[match(top4,rownames(mat_plot)),]


# order by group then infectious status --------------------------------------
ann_col <- data.frame(Group=factor(paste0("Group.",as.integer(factor(pData(obj)$GOM)))),
                      Status=factor(pData(obj)$anyinf6m,
                                          levels=c(0,1),
                                          labels=c("Not_infected", "Infected")),
                      row.names = colnames(mat_plot))
ann_colors = list(
  Group=c(Group.1=brewer.pal(12,"Set3")[1],
            Group.2=brewer.pal(12,"Set3")[2]),
  Status=c(Not_infected=brewer.pal(12,"Paired")[9],
           Infected=brewer.pal(12,"Paired")[10]))

heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)

mat_plot2 <- mat_plot[,order(pData(obj)$GOM, pData(obj)$anyinf6m)]
library(pheatmap)
pheatmap(mat_plot2,
         main="OTU Abundance Heatmap",
         fontsize_row = 3,
         fontsize_col = 3,
         border_color = NA,
         annotation_col = ann_col,
         cluster_rows = F, cluster_cols = F,
         col = heatmapCols,
         annotation_colors = ann_colors)


# order by group  -------------------------------------------------------------------
ann_col <- data.frame(Group=factor(paste0("Group.",as.integer(factor(pData(obj)$GOM)))),
                      # Status=factor(pData(obj)$anyinf6m,
                      #               levels=c(0,1),
                      #               labels=c("Not_infected", "Infected")),
                      row.names = colnames(mat_plot))
ann_colors = list(
  Group=c(Group.1=brewer.pal(12,"Set3")[1],
          Group.2=brewer.pal(12,"Set3")[2]))
mat_plot2 <- mat_plot[,order(pData(obj)$GOM)]
pheatmap(mat_plot2,
         main="OTU Abundance Heatmap",
         fontsize_row = 3,
         fontsize_col = 3,
         border_color = NA,
         annotation_col = ann_col,
         cluster_rows = F, cluster_cols = F,
         col = heatmapCols,
         annotation_colors = ann_colors)



# order by infectious status --------------------------------------
ann_col <- data.frame(Status=factor(pData(obj)$anyinf6m,
                                    levels=c(0,1),
                                    labels=c("Not_infected", "Infected")),
                      row.names = colnames(mat_plot))
ann_colors = list(
  Status=c(Not_infected=brewer.pal(12,"Paired")[9],
           Infected=brewer.pal(12,"Paired")[10]))
mat_plot2 <- mat_plot[,order(pData(obj)$anyinf6m)]
library(pheatmap)
pheatmap(mat_plot2,
         main="OTU Abundance Heatmap",
         fontsize_row = 3,
         fontsize_col = 3,
         border_color = NA,
         annotation_col = ann_col,
         cluster_rows = F, cluster_cols = F,
         col = heatmapCols,
         annotation_colors = ann_colors)
dev.off()






# Compare shannon diversity between groups ---------------------------------------------
library(data.table)
shannon_values <- fread("~/Dropbox/current-projects/nasalMicrobiomeManuscript/draft-results/shannon_values.tsv")
obj <- readRDS("data/nasal_GOM.rds")

all.equal(shannon_values$V1, colnames(obj))
shannon_values$cluster <- obj$GOM
wilcox.test(shannon_values$V2 ~ shannon_values$cluster)


## Compute Bradley Curtis index??
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



pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/beta_diversity_mean.pdf",
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
  geom_boxplot(width=.25) + ylab("Bray-Curtis similarity average across samples") +
  scale_x_discrete(name = "Cluster membership",
                   labels = c("Group 1", "Group 2")) +
  guides(fill=guide_legend(title="")) +
  ggtitle("Beta diversity")
wilcox.test(df$beta~df$clust)
dev.off()



pdf("~/Dropbox/GitHub/nasalmicrobiome/output-manuscript/beta_diversity_raw.pdf",
    width=4,heigh=3)
beta_group1 <- do.call(rbind, lapply(1:length(group1), function(n) {
  beta_bray[which(rownames(beta_bray)==group1[n]),
                 which(colnames(beta_bray) %in% group1)]
}))
beta_group1_raw <- beta_group1[upper.tri(beta_group1)]

beta_group2 <- do.call(rbind, lapply(1:length(group2), function(n) {
  beta_bray[which(rownames(beta_bray)==group2[n]),
                 which(colnames(beta_bray) %in% group2)]
}))
beta_group2_raw <- beta_group2[upper.tri(beta_group2)]

df <- data.frame(clust=rep(c(1,2), times=c(length(beta_group1_raw),
                                           length(beta_group2_raw))),
                 beta=c(beta_group1_raw,beta_group2_raw))
df$clust <- factor(df$clust)
levels(df$clust) <- c("Group 1", "Group 2")
library(ggplot2)
ggplot(df, aes(x=clust, y = beta, fill = clust)) +
  geom_violin(col = "gray50") +
  scale_fill_manual(values = c(brewer.pal(12,"Set3")[1],
                               brewer.pal(12,"Set3")[2])) +
  geom_boxplot(width=.25) + ylab("Bray-Curtis similarity average across samples") +
  scale_x_discrete(name = "Cluster membership",
                   labels = c("Group 1", "Group 2")) +
  guides(fill=guide_legend(title="")) +
  ggtitle("Beta diversity")
wilcox.test(df$beta~df$clust)
dev.off()


# beta_group12 <- do.call(rbind, lapply(1:length(group1), function(n) {
#   beta_bray[which(rownames(beta_bray)==group1[n]),
#             which(colnames(beta_bray) %in% group2)]
# }))
# mean(beta_group12)
