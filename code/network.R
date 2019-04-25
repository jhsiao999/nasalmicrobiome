# library(devtools)
# install_github('zdk123/SpiecEasi')
library(SpiecEasi)
library(metagenomeSeq)
library(phyloseq)
library(dplyr)
#setwd("/Users/paulsoj1/Desktop/tmp/nasal")
setwd("/Users/joycehsiao/Dropbox/GitHub/nasalmicrobiome/data")
x = readRDS("finalMRobj.rds")
x = filterData(x,present=100)
x = x[-grep("Chlor",fData(x)$class),]
x = x[-grep("TM",fData(x)$class),]
cnts= MRcounts(x)
gom1 = cnts[,which(pData(x)$GOM==1)] %>% t
gom2 = cnts[,which(pData(x)$GOM==2)] %>% t

nc = 3
se.est1 <- spiec.easi(gom1,lambda.min.ratio=1e-2, nlambda=20,
                      icov.select.params=list(rep.num=20, ncores=nc))
se.est2 <- spiec.easi(gom2,lambda.min.ratio=1e-2, nlambda=20,
                      icov.select.params=list(rep.num=20, ncores=nc))

library(igraph)

pdf("~/Desktop/nasal_graphs.pdf",width=8,height=8)
b = MRexperiment2biom(x[,which(pData(x)$GOM==1)])
p = phyloseq::import_biom(b)
gom1g <- adj2igraph(se.est1$refit, vertex.attr=list(name=taxa_names(p)))
phyloseq::plot_network(gom1g,p, type='taxa', color="Rank3", label=NULL)


b2 = MRexperiment2biom(x[,which(pData(x)$GOM==2)])
p2 = phyloseq::import_biom(b2)
gom2g <- adj2igraph(se.est2$refit, vertex.attr=list(name=taxa_names(p)))
phyloseq::plot_network(gom2g,p2, type='taxa', color="Rank3", label=NULL)


dev.off()


dd.gom1 <- degree_distribution(gom1g, cumulative=FALSE)
dd.gom2 <- degree_distribution(gom2g, cumulative=FALSE)
## ave degree
sum(seq_along(dd.gom1)*dd.gom1)-1
## [1] 2.367347
sum(seq_along(dd.gom2)*dd.gom2)-1
## [1] 1.37037
## plot degree distributions
plot(seq_along(dd.gom1)-1, dd.gom1, type='b', xlim=c(0,8),ylim=c(0,0.55),
     ylab="Frequency", xlab="Degree", col='red')
points(seq_along(dd.gom2)-1, dd.gom2 , type='b')
legend("topright", c("Cluster2", "Cluster1"), col=c("black", "red"), pch=1, lty=1)
