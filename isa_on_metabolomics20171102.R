## metabolomics data analysis using 
## - serum.nmr.zscore.all.colaus1.20160830.csv
## - urine.nmr.bin.zscored.colaus1.20170227.csv

# Specify OS system: ------------------------------------------------------

### for windows
# pathOS <- "C://"
### for linux
pathOS <- "/media/mirjam/OS/"

# Import package ----------------------------------------------------------
library(isa2)
library(gplots)
library(ggplot2)

# import ISA (and other) functions  ---------------------------------------
source(paste0(pathOS, 'Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project/ISA_functions.R'))
source(paste0(pathOS, 'Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project/functions.R'))

# Import Data -------------------------------------------------------------
## serum data
serum <- as.matrix(read.csv("../data/metabolome/serum.nmr.zscore.all.colaus1.20160830.csv", h = F))
rownames(serum) <- c("id", paste("id", serum[2:nrow(serum),1], sep = "_"))
serum <- serum[,-1]
colnames(serum) <- paste("ppm", serum[1,], sep ="_")
serum <- serum[-1,]
serum <- t(na.omit(t(serum)))

## urine data
# urine <- as.matrix(read.csv("../data/metabolome/urine.nmr.bin.zscored.colaus1.20170227.csv", h = F))
urine <- as.matrix(read.csv("../data/metabolome/urine.nmr.bin.decov.colaus1.20171115.csv", h = F))
# urine <- as.matrix(read.csv("../data/metabolome/urine.nmr.focus700.zscored.colaus1.20161205.csv", h = F))
rownames(urine) <- c("id", paste("id", urine[2:nrow(urine),1], sep = "_"))
urine <- urine[,-1]
colnames(urine) <- paste("ppm", urine[1,], sep ="_")
urine <- urine[-1,]
urine <- t(na.omit(t(urine))) # removing samples
# urine <- na.omit(urine) # removing features ot a good idea)

## select matching id (790 individuals having serum and urine metabolomics)
serum <- serum[rownames(serum) %in% rownames(urine),]
serum <- serum[order(rownames(serum)),]
urine <- urine[rownames(urine) %in% rownames(serum),]
urine <- urine[order(rownames(urine)),]


# Data Normalization ------------------------------------------------------
## data are already normalised!!
# save(urine, urine_log, serum, serum_log, file = paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/Rdata/metabo_data", 
#                                                        Sys.Date(), ".Rdata"))

# # Look to the data!! ------------------------------------------------------
# # Know the story around the data!
# # Ask concrete questions!
# ## Data summary and visualization
# # heatmap
# image(as.matrix(serum), xlab = "Individuals", ylab = "Features", main = "serum.nmr.zscore.all.colaus1.20160830.csv")
# image(as.matrix(urine), xlab = "Individuals", ylab = "Features", main = "urine.nmr.bin.zscored.colaus1.20170227.csv")
# # histograms
# hist(as.numeric(as.matrix(serum)), main = "Intensity distribution of \n 
#      serum.nmr.zscore.all.colaus1.20160830.csv", xlab = "Intensity", 80)
# hist(as.numeric(as.matrix(urine)), main = "Intensity distribution of \n 
#      urine.nmr.bin.zscored.colaus1.20170227.csv", xlab = "Intensity", 100)
# 

## phenotype
phenotype <- read.table("../data/phenotype/phenotype_mirjam_20171101.txt", h = T)
all_variables <- colnames(phenotype)[2:ncol(phenotype)]
variables     <- c("gluc", "trig", "hdlch", "chol", "ldlch", "sbp", "cr", "cru",
                   "gfr", "insulin")
co_variables  <- all_variables[! all_variables %in% variables]


# isa analysis ------------------------------------------------------------
# create a row.seeds vector (reproducible results!!)
set.seed(10) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
row.seeds <- generate.seeds(length=nrow(urine), 
                            count=100, sparsity=sparsity)
# serum.isa    <- isa.run(data = serum, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds, cor.limit = 0.6)
# urine.isa    <- isa.run(data = urine, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds, cor.limit = 0.6)
# urine.isa2   <- isa.run(data = urine, thr.col = seq(2.0, 3.5, by = 0.1), row.seeds = row.seeds, cor.limit = 0.6)
urine.isa3   <- isa.run(data = urine, thr.col = seq(2.4, 3.6, by = 0.2), 
                        row.seeds = row.seeds, cor.limit = 0.6)
urine.isa0.9 <- isa.run(data = urine, thr.col = seq(2.4, 3.6, by = 0.2), 
                        row.seeds = row.seeds, cor.limit = 0.9)

## after the combination of different threshold analysis
## isa using the urine.nmr.bin.zscored.colaus1.20170227.csv file
urine.isa20171122 <- isa.run(data = urine, thr.col = seq(3,5.4, by = 0.4), 
                             thr.row = seq(0.2, 3.8, by = 0.2), 
                             row.seeds = row.seeds, cor.limit = 0.6)


## modules information
urine.modules20171122 <- isaModules(urine.isa20171122)
par(mfrow=c(1,2))
hist(urine.modules20171122$colGroups, col = "blue", main = "Number of features", 
     xlab = "urine", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
hist(urine.modules20171122$rowGroups, col = "yellow", main = "Number of ndividuals", 
     xlab = "urine", cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
par(mfrow=c(1,1))

## metabomatching input
input <- metabomatching_input(urine, urine.isa20171122, 
                              ppm = as.numeric(sapply(strsplit(colnames(urine), "_"), "[[", 2)))
dir.create("../result/metabomatching/ps.urine.bin.zscored.20171122")
write.table(input, "../result/metabomatching/ps.urine.bin.zscored.20171122/urine20171122.pseudospectrum.tsv", quote = F, sep = "\t", row.names = F)

## gwas input
input.gwas <- last.iteration(urine, urine.isa20171122, 
                              ppm = as.numeric(sapply(strsplit(colnames(urine), "_"), "[[", 2)))
dir.create("../result/metabomatching/ps.urine.bin.zscored.20171122")
write.table(input, "../result/metabomatching/ps.urine.bin.zscored.20171122/urine20171122.pseudospectrum.tsv", quote = F, sep = "\t", row.names = F)





# run isa -----------------------------------------------------------------
files <- c("serum", "urine")
path1  <- "../result/metabomatching/"
path2 <- "../result/linear_regression/"

for(data in files){
  # create a row.seeds vector (reproducible results!!)
  set.seed(10) # to have the same random seeds! not using the random.seeds function!
  sparsity <- rep(2* (c(1,5,25,125)), length=100)
  row.seeds <- generate.seeds(length=nrow(eval(parse(text = data))), 
                              count=100, sparsity=sparsity)
  ## run isa
  data.isa   <- isa.run(data = as.matrix(eval(parse(text = data))), 
                        thr.col = seq(1.5, 2, by = 0.1), 
                        row.seeds = row.seeds, cor.limit = 0.6)
  ## metabomatching input
  data.input <- metabomatching_input(eval(parse(text = data)), data.isa, 
                                     as.numeric(sapply(strsplit(colnames(eval(parse(text = data))), "_"),
                                                       "[[", 2)))
  
  write.table(data.input, file = paste0(path1, 
                                        data, ".tsv"), 
              quote = F, sep = "\t", row.names = F)
  ## correlation results
  score.id <- as.data.frame(last.iteration(as.matrix(eval(parse(text = data))), data.isa))
  score.id$ID <- rownames(eval(parse(text = data)))
  ## phenotype data (select only women or men)
  phen <- phenotype[na.omit(match(phenotype$ID, score.id$ID)),]
  
  linear_regression_result <- matrix(NA, ncol(score.id)-1, 10)
  for(j in 1:length(variables)){
    for (i in 1:(ncol(score.id)-1)){
      p <- summary(lm(phen[,variables[j]] ~ score.id[,i] + sex + age + bmi, data = phen))$coefficients[2,4]
      linear_regression_result[i,j] <- round(-log(p),2)
    }
  }
  
  colnames(linear_regression_result) <- variables
  rownames(linear_regression_result) <- paste0("Module", 1:(ncol(score.id)-1))
  
  png(paste0(path2, data, Sys.Date(), ".png"), 
      width = 23, height = 20, units = 'cm', res = 300)
  heatmap.2(linear_regression_result, cellnote = linear_regression_result, 
            notecex= 0.7, margins = c(7,7), main = paste0("Linear regression result ", data), 
            cexRow = 0.8, cexCol = 1.2, xlab = "", ylab = "",
            notecol="black", density.info="none", trace="none", dendrogram = "none", 
            col = c("forestgreen", "red"), breaks = c(0, 3, max(linear_regression_result)), 
            Colv="NA", Rowv = "NA", keysize = 0.5, key =F)
  dev.off()
  
}


# isa results analysis ----------------------------------------------------

urine.modules <- isaModules(urine.isa)
serum.modules <- isaModules(serum.isa)

input <- metabomatching_input(urine, urine.isa3, 
                 ppm = as.numeric(sapply(strsplit(colnames(urine), "_"), "[[", 2)))
dir.create("../result/metabomatching/ps.urine20171114.default")
write.table(input, "../result/metabomatching/ps.urine20171114.default/urine.pseudospectrum.tsv", quote = F, sep = "\t", row.names = F)
input <- na.omit(metabomatching_input(serum, serum.isa, 
                 ppm = as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2))))
write.table(input, "../result/metabomatching/serum.tsv", quote = F, sep = "\t", row.names = F)


urine.isa0.9
cor.limit <- 0.7
abs(cor(urine.isa0.9$rows))[1:10, 1:6]
abs(cor(urine.isa0.9$columns))[1:10, 1:6]
cm <- pmin(abs(cor(urine.isa0.9$rows)), abs(cor(urine.isa0.9$columns)))
cm[ lower.tri(cm, diag=TRUE) ] <- 0


# correlation between score!!! 
# Not between individuals and fetures present in the modules!!
# when correlation limit is set, it refears to how similar the scores are bwtween two modules!
# how well it correlate

id1 <- isaRowNames(urine, urine.isa0.9,50)
id2 <- isaRowNames(urine, urine.isa0.9,51)
c(sum(id1 %in% id2), 
  sum(id2 %in% id1))[which.max(c(sum(id1 %in% id2), 
                                 sum(id2 %in% id1)))]/c(length(id2), length(id1))[which.max(c(sum(id1 %in% id2), 
                                        sum(id2 %in% id1)))]
cor(urine.isa0.9$rows[,50], urine.isa0.9$rows[,51])
plot(urine.isa0.9$rows[,50], urine.isa0.9$rows[,51])

ppm1 <- isaColNames(urine, urine.isa0.9,50)
ppm2 <- isaColNames(urine, urine.isa0.9,51)
c(sum(ppm1 %in% ppm2), 
  sum(ppm2 %in% ppm1))[which.max(c(sum(ppm1 %in% ppm2), 
                                 sum(ppm2 %in% ppm1)))]/c(length(ppm2), length(ppm1))[which.max(c(sum(ppm1 %in% ppm2), 
                                                                                              sum(ppm2 %in% ppm1)))]
cor(urine.isa0.9$columns[,50], urine.isa0.9$columns[,51])
plot(urine.isa0.9$columns[,50], urine.isa0.9$columns[,51])


## find group of modules that are similar (correlation between 0.6 and 0.9)
## 
data <- "urine"
isaresults <- urine.isa0.9
i=0
while(ncol(isaresults$rows) >= 1){
  i <- i + 1
  corr_modules       <- pmin(abs(cor(isaresults$rows)), abs(cor(isaresults$columns)))
  corr_modules[lower.tri(corr_modules, diag=TRUE) ] <- 0
  
  similar_modules     <- corr_modules[1,] >= 0.6
  similar_modules[1]  <- TRUE
  isaresult           <- list()
  isaresult$rows      <- isaresults$rows[,similar_modules, drop=FALSE]
  isaresult$columns   <- isaresults$columns[,similar_modules,drop=FALSE]
  isaresult$seeddata  <- isaresults$seeddata[similar_modules,]
  isaresults$rows     <- isaresults$rows[,!similar_modules,drop=FALSE]
  isaresults$columns  <- isaresults$columns[,!similar_modules,drop=FALSE]
  isaresults$seeddata <- isaresults$seeddata[!similar_modules,]
  data.input          <- metabomatching_input(eval(parse(text = data)), isaresult, 
                                             as.numeric(sapply(strsplit(colnames(eval(parse(text = data))), "_"),
                                                               "[[", 2)))
  if(ncol(data.input) == 4){
    colnames(data.input) <- c("shift", "beta", "se", "p")
  } 
  path <- sprintf("../result/metabomatching/modules_groups/ps.group%03d", i)
  dir.create(path)
  setwd(path)
  write.table(data.input, paste0(sprintf("group%03d", i), ".pseudospectrum.tsv"), 
              quote = F, sep = "\t", row.names = F)
  setwd("/media/mirjam/OS/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project")
}


# Run isa with combination of thresholds ----------------------------------

### !!!!!!! change directory!! #####
cor.limit <- 0.6
# setwd(paste0("../result/isa/urine/urine.nmr.bin.zscored/corr", cor.limit))
setwd(paste0("../result/isa/urine/urine.nmr.bin.decov/corr", cor.limit))
# setwd(paste0("../result/isa/urine/urine.nmr.focus700.zscored/corr", cor.limit))
# setwd(paste0("../result/isa/serum/corr", cor.limit))

set.seed(10) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
row.seeds <- generate.seeds(length=nrow(urine), 
                            count=100, sparsity=sparsity)
# thr.row.list <- seq(0.4, 6, by = 0.2)
# thr.col.list <- seq(0.4, 6, by = 0.2)

thr.row.list <- c(6.2, 6.4)
thr.col.list <- c(6.2, 6.4)

normed.data <- isa.normalize(urine)
for (thr.col in thr.col.list){
  for (thr.row in thr.row.list){
    isaresults <- isa.iterate(normed.data,
                              row.seeds=row.seeds,
                              thr.row=thr.row,
                              thr.col=thr.col,
                              direction=c("updown", "updown"))
    ## Make it unique for every threshold combination
    isaresults.uni <- isa.unique(normed.data, isaresults, cor.limit = cor.limit)
    if(ncol(isaresults.uni$rows) == 0){
      isaresults <- isaresults
    } else {
      isaresults <- isaresults.uni
    }
    
    ## Filter according to robustness
    isaresults.rob <- isa.filter.robust(urine,
                        normed.data=normed.data,
                        isares=isaresults,
                        row.seeds=row.seeds)
    if(ncol(na.omit(isaresults.rob$rows)) == 0){
      isaresults <- isaresults
    } else {
      isaresults <- isaresults.rob
    }
    isa.modules <- isaModules(isaresults)
    if(nrow(isa.modules)==100){
      isa.modules <- isa.modules[1,]
    }
    write.table(isa.modules, paste0("module_info_col", thr.col, 
                                               "row", thr.row, ".txt"), 
                quote = F, sep = "\t", row.names = F)
    save(isaresults, file = paste0("isa_col", thr.col, 
                            "row", thr.row, ".Rdata"))
  }
}

# par(mfrow = c(length(thr.col.list), length(thr.row.list)))
# par(mar = c(0, 0, 0, 0), oma = c(5, 5, 2, 2))
# my.colors <- colorRampPalette(c("darkred", "yellow"))(n = 30)
# for (file in files) {
#   res <- read.table(paste0("../result/isa/", file), h = T, sep = "\t")
#   mn <- round(mean(res$colGroups))
#   boxplot(res$colGroups, axes = FALSE, type = "n", ylim = c(0,60), col = my.colors[mn])
# }
# 
# mtext(thr.row.list, side = 1, outer = TRUE, line = 0.5, 
#       cex = 0.8, font = 2, at = seq(0.055, 0.95, length.out = 9))
# mtext(thr.col.list, side = 2, outer = TRUE, line = -1, 
#       cex = 0.8, font = 2, at = seq(0.015, 0.95, length.out = 11))
# 
# mtext("Row thresholds", side = 1, outer = TRUE, cex = 1, line = 2.2, col = "grey20")
# mtext("Columns thresholds", side = 2, outer = TRUE, cex = 1, line = 1,
#       col = "grey20")

files <- list.files(".", pattern = ".txt")
features <- matrix(NA, length(thr.row.list), length(thr.col.list))
colnames(features) <- thr.col.list
rownames(features) <- thr.row.list
individuals <- matrix(NA, length(thr.row.list), length(thr.col.list))
colnames(individuals) <- thr.col.list
rownames(individuals) <- thr.row.list
n.modules <- matrix(NA, length(thr.row.list), length(thr.col.list))
colnames(n.modules) <- thr.col.list
rownames(n.modules) <- thr.row.list
modules.sd.col <- matrix(NA, length(thr.row.list), length(thr.col.list))
colnames(modules.sd.col) <- thr.col.list
rownames(modules.sd.col) <- thr.row.list
modules.sd.row <- matrix(NA, length(thr.row.list), length(thr.col.list))
colnames(modules.sd.row) <- thr.col.list
rownames(modules.sd.row) <- thr.row.list
for (file in files) {
  res <- read.table(file, h = T, sep = "\t")
  mn.col <- round(mean(res$colGroups), 2)
  mn.row <- round(mean(res$rowGroups), 2)
  n <- nrow(res)
  sd.col <- round(sd(res$colGroups), 2)
  sd.row <- round(sd(res$rowGroups), 2)
  thr.col <- match(strsplit(strsplit(file, "col")[[1]][2], "row")[[1]][1], thr.col.list)
  thr.row <- match(strsplit(strsplit(file, "row")[[1]][2], ".txt")[[1]][1], thr.row.list)
  features[thr.row, thr.col] <- mn.col
  individuals[thr.row, thr.col] <- mn.row
  n.modules[thr.row, thr.col] <- n
  modules.sd.col[thr.row, thr.col] <- sd.col
  modules.sd.row[thr.row, thr.col] <- sd.row
}


my.colors <- c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59",
               "#ef6548", "#d7301f", "#990000")
## number of features
png(paste0("n_features", cor.limit, Sys.Date(), ".png"), width = 23, height = 20, units = 'cm', res = 300)
# class: < 10
# 10 - 20
# 20 - 50
# 50 - 100
# 100 - 200
# 200 - 300
# > 300                   

heatmap.2(features, cellnote = round(features), notecex= 0.7, notecol="black",
          Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
          keysize = 0.5, key = F, margins = c(7,7), main = "Number of features", 
          ylab = "Row thresholds", xlab = "Col threscholds", cexRow = 1.2, 
          cexCol = 1.2, col = my.colors) ## , breaks = c(0, 10, 20, 50, 100, 200, 300, 400)
dev.off()
## number of individuals
png(paste0(paste0("n_individuals", cor.limit, Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(individuals, cellnote = round(individuals), notecex= 0.7, notecol="black",
          Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
          keysize = 0.5, key = F, margins = c(7,7), main = "Number of individuals", 
          ylab = "Row thresholds", xlab = "Col threscholds", cexRow = 1.2, 
          cexCol = 1.2, col = my.colors) ## , breaks = c(0, 5, 10, 20, 50, 100, 200, 300)
dev.off()
## number of modules
png(paste0(paste0("n_modules", cor.limit, Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(n.modules, cellnote = n.modules, notecex= 0.7, notecol="black",
          Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
          keysize = 0.5, key = F, margins = c(7,7), main = "Number of modules", 
          ylab = "Row thresholds", xlab = "Col threscholds", cexRow = 1.2, 
          cexCol = 1.2, col = my.colors)
dev.off()
## sd columns
png(paste0(paste0("sd_col", cor.limit, Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(modules.sd.col, cellnote = round(modules.sd.col), notecex= 0.7, notecol="black",
          Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
          keysize = 0.5, key = F, margins = c(7,7), main = "Standard deviation (# features)", 
          ylab = "Row thresholds", xlab = "Col threscholds", cexRow = 1.2, 
          cexCol = 1.2, col = my.colors) ## , breaks = c(0, 1, 5, 10, 20, 30, 50, 100)
dev.off()

## sd rows
png(paste0(paste0("sd_row", cor.limit, Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(modules.sd.row, cellnote = round(modules.sd.row), notecex= 0.7, notecol="black",
          Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
          keysize = 0.5, key = F, margins = c(7,7), main = "Standard deviation (# individuals)", 
          ylab = "Row thresholds", xlab = "Col threscholds", cexRow = 1.2, 
          cexCol = 1.2, col = my.colors) ## , breaks = c(0, 1, 3, 5, 10, 20, 30, 50)
dev.off()

setwd("../../../../../master_project/")
