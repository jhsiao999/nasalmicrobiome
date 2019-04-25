# for each K,
# we selected the best seed for membership assignment
# see gom-seeds.Rmd for the analysis

# this docoument outputs tables containing sample IDs for Nauder


clust_2 <- readRDS("output/gom-k2-best-seed-67.rds")
clust_3 <- readRDS("output/gom-k3-best-seed-67.rds")
clust_4 <- readRDS("output/gom-k4-best-seed-59.rds")
clust_5 <- readRDS("output/gom-k5-best-seed-307.rds")

# clust_2_tab <- data.frame(sampleID=rownames(clust_2[[1]]$omega),
#                           clust=apply(clust_2[[1]]$omega, 1,
#                                        function(x) which.max(x)))
# write.csv(clust_2_tab,
#             file = "output/clust2_sampleID.csv",
#           row.names = F)

clust_3_tab <- data.frame(sampleID=rownames(clust_3[[1]]$omega),
                          clust=apply(clust_3[[1]]$omega, 1,
                                      function(x) which.max(x)))
write.csv(clust_3_tab,
          file = "output/clust3_sampleID.csv",
          row.names = F)



clust_4_tab <- data.frame(sampleID=rownames(clust_4[[1]]$omega),
                          clust=apply(clust_4[[1]]$omega, 1,
                                      function(x) which.max(x)))
write.csv(clust_4_tab,
          file = "output/clust4_sampleID.csv",
          row.names = F)


clust_5_tab <- data.frame(sampleID=rownames(clust_5[[1]]$omega),
                          clust=apply(clust_5[[1]]$omega, 1,
                                      function(x) which.max(x)))
write.csv(clust_5_tab,
          file = "output/clust5_sampleID.csv",
          row.names = F)


