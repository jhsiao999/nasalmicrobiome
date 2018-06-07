library(knitr)
library(kableExtra)
library(dplyr)
library(metagenomeSeq)
library(CountClust)
library(parallel)


#Read in filtered data.
MRobj = readRDS("data/nasal_filtered_normed_batchcorrected.rds")
counts <- MRcounts(MRobj,norm=FALSE,log=FALSE)

# get 100 random seeds
source("code/print-prime.R")
seeds <- prime(543)

for (k in 2:5) {
  res <- vector("list", length(seeds))
  for (i in 1:100) {
    set.seed(seeds[i])
    res[[i]] <- FitGoM(t(counts), K=k, tol=1e-6)[[1]]
  }
  saveRDS(res, file = paste0("output/gom-seeds-k-",k,".rds"))
}




