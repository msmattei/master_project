
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
load(file = "../data/Rdata/metabo_data2017-10-12.Rdata") ## normalized-transformed data, normalized data (*_log)

# run isa -----------------------------------------------------------------
## Increase the number of subjects per module, reduce the total number of modules (to have less, bigger modules)

## serum
# create a row.seeds vector (reproducible results!!)
set.seed(10) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
row.seeds <- generate.seeds(length=nrow(serum), count=100, sparsity=sparsity)
## run isa
serum.isa <- isa.run(data = serum, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)
serum.isa.no.norm <- isa.run(data = as.matrix(serum_log), thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)

## urine
# create a row.seeds vector (reproducible results!!)
set.seed(122) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
row.seeds <- generate.seeds(length=nrow(urine), count=100, sparsity=sparsity)
## run isa
urine.isa <- isa.run(data = urine, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)
urine.isa.no.norm <- isa.run(data = as.matrix(urine_log), thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds)



# module information ------------------------------------------------------
# serum
serum.modules <- isaModules(serum.isa)
serum.modules.no.norm <- isaModules(serum.isa.no.norm)
# urine
urine.modules <- isaModules(urine.isa)
urine.modules.no.norm <- isaModules(urine.isa.no.norm)

## visualization
# serum
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_hist_serum",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
hist(serum.modules$colGroups, col = "red", main = "Peaks by modules", xlab = "serum", cex.lab=1.6, 
     cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(serum.modules$rowGroups, col = "orange", main = "Individuals by modules", xlab = "serum", 
     cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
dev.off()

# urine
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_hist_urine",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
hist(urine.modules$colGroups, col = "blue", main = "Peaks by modules", xlab = "urine", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(urine.modules$rowGroups, col = "yellow", main = "Individuals by modules", xlab = "urine", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
dev.off()

#  ------------------------------------------------------------------------

## identity between modules created using pre-normalized data or non-normalized data
# serum
identity.norm.non.serum <- identity(data1 = serum, data2 = serum_log, data.isa1 = serum.isa, 
                                    data.isa2 = serum.isa.no.norm, modules1 = serum.modules, 
                                    modules2 = serum.modules.no.norm, sel = 0)
# visualization
heatmap.2(identity.norm.non.serum, main = "Identity between module created using \n prenormalized or raw serum",
          notecol="black", density.info="none", trace="none", 
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, xlab = "non-normalized", ylab = "pre-normalized")

#urine
identity.norm.non.urine <- identity(data1 = urine, data2 = urine_log, data.isa1 = urine.isa, 
                                    data.isa2 = urine.isa.no.norm, modules1 = urine.modules, 
                                    modules2 = urine.modules.no.norm, sel = 0)
# visualization
heatmap.2(identity.norm.non.urine, 
          main = "Identity between module created using \n prenormalized or raw urine data",
          notecol="black", density.info="none", trace="none", 
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, xlab = "non-normalized", ylab = "pre-normalized")


# combine urine and serum dataset -----------------------------------------

# urine and serum combination
uri.ser <- cbind(urine, serum)
colnames(uri.ser) <- c(paste0("u.", colnames(urine)), paste0("s.", colnames(serum)))
# uri.ser.isa <- isa(as.matrix(uri.ser))
# to re-run the isa function with adequate parameters!! (12.10.1017)

# combination of urine and serum datasets
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_hist_uri_ser",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
hist(uri.ser.modules$colGroups, col = "blue", main = "Peaks by modules", xlab = "Serum and urine combined", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(uri.ser.modules$rowGroups, col = "yellow", main = "Individuals by modules", xlab = "Serum and urine combined", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
dev.off()



# Identity between modules ------------------------------------------------
## identity function (in ISA_functions.R file) used to calculate the percentage of identity between modules 
## either of the same biological fluid or of modules constructed using different samples (for example urine and serum)

# Overlap of individuals or ppms present in urine modules
identity.urine <- identity(data1 = urine, data2 = urine, data.isa1 = urine.isa, data.isa2 = urine.isa, 
                           modules1 = urine.modules, modules2 = urine.modules, sel = 0)
identity.urine.col <- identity(data1 = urine, data2 = urine, data.isa1 = urine.isa, data.isa2 = urine.isa, 
                               modules1 = urine.modules, modules2 = urine.modules, sel = 0, Col = TRUE)
# Overlap of individuals or ppms present in serum modules
identity.serum <- identity(data1 = serum, data2 = serum, data.isa1 = serum.isa, data.isa2 = serum.isa, 
                           modules1 = serum.modules, modules2 = serum.modules, sel = 0)
identity.serum.col <- identity(data1 = serum, data2 = serum, data.isa1 = serum.isa, data.isa2 = serum.isa, 
                               modules1 = serum.modules, modules2 = serum.modules, sel = 0, Col = TRUE)

# look for similarity of modules between urine and serum 
# (selection of row and columns containings at least one value > 0.3)
identity.urine.serum <- identity(data1 = urine, data2 = serum, data.isa1 = urine.isa, data.isa2 = serum.isa, 
                                 modules1 = urine.modules, modules2 = serum.modules, sel = 0.3)
# Overlap of individuals in modules created using the combined dataset (urine-serum)
identity.uri.ser <- identity(data1 = uri.ser, data2 = uri.ser, data.isa1 = uri.ser.isa, data.isa2 = uri.ser.isa, 
                             modules1 = uri.ser.modules, modules2 = uri.ser.modules, sel = 0)
identity.uri.ser.sel <- identity(data1 = uri.ser, data2 = uri.ser, data.isa1 = uri.ser.isa, data.isa2 = uri.ser.isa, 
                                 modules1 = uri.ser.modules, modules2 = uri.ser.modules, sel = 0.8)


## Visualization
#urine
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/individuals_similarity_urine",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.urine,
          main = "Urinary module having similar individuals",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/ppm_similarity_urine",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.urine.col,
          main = "Urinary module having similar ppm values",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

# example of similar urine modules: 4 subjects of module 149 are present in module 147 (5 subjects). 
# 21 of the 23 ppms of module 147 are present in module 149 (24 ppms)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_urine_ex",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
isa2image(urine, data.isa = urine.isa, type = "isa", n = 147)
isa2image(urine, data.isa = urine.isa, type = "isa", n = 149)
dev.off()

# serum
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/individuals_similarity_serum",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.serum,
          main = "Serum module having similar individuals",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/ppm_similarity_serum",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.serum.col,
          main = "Serum module having similar ppm values",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks, labRow = FALSE, labCol = FALSE)
dev.off()

## urine versus serum modules
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_urine_vs_serum",
                  Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.urine.serum, cellnote = identity.urine.serum, notecex= 0.7,
          main = "Matching individuals in serum and urine modules", margins = c(7,7), xlab = "Serum", ylab = "Urine",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5, cexRow = 1.2, cexCol = 1.2,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks)
dev.off()




# PPA ---------------------------------------------------------------------
ser_ur <- list(as.matrix(t(serum)), as.matrix(t(urine)))
metabo.ppa.def <- ppa(ser_ur) # resulting in too many modules
metabo.ppa <- ppa(ser_ur, thr.row1 = seq(1.5, 2.5, by = 0.1), thr.row2 = seq(1.5, 2.5, by = 0.1), thr.col = seq(2,3, by = 0.1))

save(metabo.ppa.def, file = paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ppa20170812.RData"))
load(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ppa20170812.RData"))
metabo.module.info <- isaModules(data = metabo.ppa.def, type = "ppa")

hist(metabo.module.info$colGroups)
hist(metabo.module.info2$colGroups)
hist(metabo.module.info$row1Groups)

## ppa plot

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa.def,
          n = 29, name1 = "serum", name2 = "urine")

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa2,
          n = 44, name1 = "serum", name2 = "urine")




