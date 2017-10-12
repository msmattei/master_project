
# load packages -----------------------------------------------------------
library(isa2)
library(gplots)

# Specify OS system: ------------------------------------------------------
### for windows
pathOS <- "C://"
### for linux
pathOS <- "/media/mirjam/OS/"

# ISA (and other) functions  ----------------------------------------------
source(paste0(pathOS, 'Mimi/Stage_CBG/2.EXPRESSION_MODULE/expression_module/ISA_functions.R'))
source(paste0(pathOS, 'Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project/functions.R'))

# load data ---------------------------------------------------------------
## use of prenormalized serum maybe not necessary?? 



# run isa -----------------------------------------------------------------

## serum
# create a row.seeds vector (reproducible results!!)
set.seed(10) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
row.seeds <- generate.seeds(length=nrow(serum), count=100, sparsity=sparsity)
## run isa
serum.isa <- isa.run(data = serum, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)
serum.isa.no.norm <- isa.run(data = as.matrix(serum_log), thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)

# module information
serum.modules <- isaModules(serum.isa)
serum.modules.no.norm <- isaModules(serum.isa.no.norm)
# identity between modules created using pre-normalized data or non-normalized data
identity.norm.non.serum <- identity(data1 = serum, data2 = serum_log, data.isa1 = serum.isa, 
                              data.isa2 = serum.isa.no.norm, modules1 = serum.modules, 
                              modules2 = serum.modules.no.norm, sel = 0)
# visualization
heatmap.2(identity.norm.non.serum, main = "Identity between module created using \n prenormalized or raw serum",
          notecol="black", density.info="none", trace="none", 
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, xlab = "non-normalized", ylab = "pre-normalized")


## urine
# create a row.seeds vector (reproducible results!!)
set.seed(122) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
row.seeds <- generate.seeds(length=nrow(urine), count=100, sparsity=sparsity)
## run isa
urine.isa <- isa.run(data = urine, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)
urine.isa.no.norm <- isa.run(data = as.matrix(urine_log), thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)

## module information
urine.modules <- isaModules(urine.isa)
urine.modules.no.norm <- isaModules(urine.isa.no.norm)
## identity between modules created using pre-normalized data or non-normalized data
identity.norm.non.urine <- identity(data1 = urine, data2 = urine_log, data.isa1 = urine.isa, 
                              data.isa2 = urine.isa.no.norm, modules1 = urine.modules, 
                              modules2 = urine.modules.no.norm, sel = 0)
# visualization
heatmap.2(identity.norm.non.urine, 
          main = "Identity between module created using \n prenormalized or raw urine data",
          notecol="black", density.info="none", trace="none", 
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, xlab = "non-normalized", ylab = "pre-normalized")





