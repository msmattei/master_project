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
urine.nmr.bin.zscored <- as.matrix(read.csv("../data/metabolome/urine.nmr.bin.zscored.colaus1.20170227.csv", h = F))
urine.nmr.bin.decov   <- as.matrix(read.csv("../data/metabolome/urine.nmr.bin.decov.colaus1.20171115.csv", h = F))
urine.nmr.focus700.zscored    <- as.matrix(read.csv("../data/metabolome/urine.nmr.focus700.zscored.colaus1.20161205.csv", h = F))
# rownames(urine) <- c("id", paste("id", urine[2:nrow(urine),1], sep = "_"))
# urine <- urine[,-1]
# colnames(urine) <- paste("ppm", urine[1,], sep ="_")
# urine <- urine[-1,]
# urine <- t(na.omit(t(urine))) # removing samples
# ## urine <- na.omit(urine) # removing features not a good idea)

## select matching id (790 individuals having serum and urine metabolomics)
# serum <- serum[rownames(serum) %in% rownames(urine),]
# serum <- serum[order(rownames(serum)),]
# urine <- urine[rownames(urine) %in% rownames(serum),]
# urine <- urine[order(rownames(urine)),]

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
# serum.isa    <- isa.run(data = serum, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds, "0", substr(cor.limit,3,3) = 0.6)
# urine.isa    <- isa.run(data = urine, thr.col = seq(1.5, 2, by = 0.1), row.seeds = row.seeds, "0", substr(cor.limit,3,3) = 0.6)
# urine.isa2   <- isa.run(data = urine, thr.col = seq(2.0, 3.5, by = 0.1), row.seeds = row.seeds, "0", substr(cor.limit,3,3) = 0.6)
urine.isa3   <- isa.run(data = urine, thr.col = seq(2.4, 3.6, by = 0.2), 
                        row.seeds = row.seeds, "0", substr(cor.limit,3,3) = 0.6)
urine.isa0.9 <- isa.run(data = urine, thr.col = seq(2.4, 3.6, by = 0.2), 
                        row.seeds = row.seeds, "0", substr(cor.limit,3,3) = 0.9)

## after the combination of different threshold analysis
## isa using the urine.nmr.bin.zscored.colaus1.20170227.csv file
urine.isa20171122 <- isa.run(data = urine, thr.col = seq(3,5.4, by = 0.4), 
                             thr.row = seq(0.2, 3.8, by = 0.2), 
                             row.seeds = row.seeds, "0", substr(cor.limit,3,3) = 0.6)


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




# family modules ----------------------------------------------------------

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

set.seed(10) # to have the same random seeds! not using the random.seeds function!
sparsity <- rep(2* (c(1,5,25,125)), length=100)
thr.row.list <- format(seq(0.4, 6.4, by = 0.2), 1)
thr.col.list <- format(seq(0.4, 6.4, by = 0.2), 1)
urine.list <- c("urine.nmr.bin.zscored", "urine.nmr.bin.decov", "urine.nmr.focus700.zscored")

for(i in 1:3){
  urine <- eval(parse(text = urine.list[i]))
  rownames(urine) <- c("id", paste("id", urine[2:nrow(urine),1], sep = "_"))
  urine <- urine[,-1]
  colnames(urine) <- paste("ppm", urine[1,], sep ="_")
  urine <- urine[-1,]
  for(cor.limit in c(0.3, 0.6, 0.9)){
    setwd(paste0("../result/isa/urine/", urine.list[i], "/corr", cor.limit))
    row.seeds <- generate.seeds(length=nrow(urine), count=100, sparsity=sparsity)
    normed.data <- isa.normalize(urine)
    for (thr.col in thr.col.list){
      for (thr.row in thr.row.list){
        isaresults <- isa.iterate(normed.data,
                                  row.seeds=row.seeds,
                                  thr.row=as.numeric(thr.row),
                                  thr.col=as.numeric(thr.col),
                                  direction=c("updown", "updown"))
        ## Make it unique for every threshold combination
        isaresults <- isa.unique(normed.data, isaresults, cor.limit = cor.limit)
        
        ## Filter according to robustness
        isaresults <- isa.filter.robust(urine,
                                        normed.data=normed.data,
                                        isares=isaresults,
                                        row.seeds=row.seeds)
        if(ncol(na.omit(isaresults$rows)) == 0){
          isa.modules <- matrix(NA, 0, 5)
          colnames(isa.modules) <- c("colGroups", "rowGroups", "rob", "thr.col", "thr.row")
        } else {
          isa.modules <- isaModules(isaresults)
        }
        
        write.table(isa.modules, paste0("info_modules/module_info_col", thr.col, 
                                        "row", thr.row, ".txt"), 
                    quote = F, sep = "\t", row.names = F)
        save(isaresults, file = paste0("Rdata/isa_col", thr.col, 
                                       "row", thr.row, ".Rdata"))
      }
    }
    setwd("../../../../../master_project/")
  }
  
}


# boxplot multiple --------------------------------------------------------

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

# thr.row.list <- seq(0.4, 6.4, by = 0.2)
# thr.col.list <- seq(0.4, 6.4, by = 0.2)


# plots -------------------------------------------------------------------

my.colors <- c("white", "#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59",
               "#ef6548", "#d7301f", "#990000")
my.colors2 <- c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59",
               "#ef6548", "#d7301f", "#990000")

for(cor.limit in c(0.3, 0.6, 0.9)){
  paths <- c(paste0("../result/isa/urine/urine.nmr.bin.zscored/corr", cor.limit, "/info_modules/"),
             paste0("../result/isa/urine/urine.nmr.bin.decov/corr", cor.limit, "/info_modules/"),
             paste0("../result/isa/urine/urine.nmr.focus700.zscored/corr", cor.limit, "/info_modules/"))
    for(path in paths){
    files <- list.files(path, pattern = ".txt")
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
      res <- read.table(paste0(path, file), h = T, sep = "\t")
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
    features[is.na(features)] <- NA
    individuals[is.na(individuals)] <- NA
    
    ## number of features
    png(paste0(strsplit(path, "info_modules/")[[1]], "n_features", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(features, cellnote = round(features), notecex = 0.8, notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Number of features", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors, breaks = c(0, 0.9, 10, 20, 50, 100, 200, 300, 400))
    dev.off()
    ## number of individuals
    png(paste0(strsplit(path, "info_modules/")[[1]], "n_individuals", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(individuals, cellnote = round(individuals), notecex = 0.8, notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Number of individuals", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors, breaks = c(0, 0.9, 10, 20, 50, 100, 200, 300, 400))
    dev.off()
    ## number of modules
    png(paste0(strsplit(path, "info_modules/")[[1]], "n_modules", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(n.modules, cellnote = n.modules, notecex = 0.8, notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Number of modules", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors, breaks = c(0, 0.9, 5, 10, 15, 20, 30, 40, 50))
    dev.off()
    
    ## sd columns
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_col", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.col, cellnote = round(modules.sd.col), notecex = 0.8, notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Standard deviation (# features)", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors2)
    dev.off()
    
    ## sd rows
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_row", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.row, cellnote = round(modules.sd.row), notecex = 0.8, notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Standard deviation (# individuals)", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors2)
    dev.off()
    
    ## sd columns corrected by features
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_col_correct", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.col/features, cellnote = round(1/(modules.sd.col/features)), notecex = 0.8, 
              notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Standard deviation (# features)", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors2)
    dev.off()
    
    ## sd rows corrected by number of individuals
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_row_correct", "0", substr(cor.limit,3,3), Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.row/individuals, cellnote = round(1/(modules.sd.row/individuals)), 
              notecex = 0.8, notecol="black",
              Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), main = "Standard deviation (# individuals)", 
              ylab = "Row thresholds", xlab = "Column thresholds", cexRow = 1.2, 
              cexCol = 1.2, col = my.colors2)
    dev.off()
    
    
    #### small figures
    ## number of features
    png(paste0(strsplit(path, "info_modules/")[[1]], "n_features", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(features, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), 
              col = my.colors, breaks = c(0, 0.9, 10, 20, 50, 100, 200, 300, 400), 
              labRow = NA, labCol = NA,
              srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 2,
              sepcolor = "white", colsep = c(7, 13, 19, 25), 
              rowsep = c(7, 13, 19, 25))
    dev.off()
    ## number of individuals
    png(paste0(strsplit(path, "info_modules/")[[1]], "n_individuals", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(individuals, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), 
              col = my.colors, breaks = c(0, 0.9, 10, 20, 50, 100, 200, 300, 400), 
              labCol = NA, labRow = NA, 
              srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 2,
              sepcolor = "white", colsep = c(7, 13, 19, 25), 
              rowsep = c(7, 13, 19, 25))
    dev.off()
    ## number of modules
    if(length(grep("urine.nmr.focus700.zscored", path))==1){
      png(paste0(strsplit(path, "info_modules/")[[1]], "n_modules", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
          width = 23, height = 20, units = 'cm', res = 300)
      heatmap.2(n.modules, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
                keysize = 0.5, key = F, margins = c(5,5), 
                col = my.colors, breaks = c(0, 0.9, 5, 10, 15, 20, 30, 40, 50),
                labRow = c(0.4, NA,NA,NA,NA,NA, 
                           1.6, NA,NA,NA,NA,NA, 
                           2.8, NA,NA,NA,NA,NA, 
                           4.0, NA,NA,NA,NA,NA, 
                           5.2, NA,NA,NA,NA,NA, 
                           6.4), 
                labCol = NA,
                srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 3.5, cexCol = 2,
                sepcolor = "white", colsep = c(7, 13, 19, 25), 
                rowsep = c(7, 13, 19, 25)) 
      dev.off()
    } else {
      png(paste0(strsplit(path, "info_modules/")[[1]], "n_modules", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
          width = 23, height = 20, units = 'cm', res = 300)
      heatmap.2(n.modules, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
                keysize = 0.5, key = F, margins = c(5,5), 
                col = my.colors, breaks = c(0, 0.9, 5, 10, 15, 20, 30, 40, 50),
                labRow = NA, labCol = NA,
                srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 4, cexCol = 2,
                sepcolor = "white", colsep = c(7, 13, 19, 25), 
                rowsep = c(7, 13, 19, 25)) 
      dev.off()
    }
    ## sd columns
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_col", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.col, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), col = my.colors2, 
              labRow = NA, labCol = NA,
              srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 2,
              sepcolor = "white", colsep = c(7, 13, 19, 25), 
              rowsep = c(7, 13, 19, 25))
    dev.off()
    
    ## sd rows
    if(length(grep("urine.nmr.bin.zscored", path))==1){
      png(paste0(strsplit(path, "info_modules/")[[1]], "sd_row", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
          width = 23, height = 20, units = 'cm', res = 300)
      heatmap.2(modules.sd.row, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
                keysize = 0.5, key = F, margins = c(5,5), 
                col = my.colors2, labRow = NA, labCol = c(0.4, NA,NA,NA,NA,NA, 
                                                          1.6, NA,NA,NA,NA,NA, 
                                                          2.8, NA,NA,NA,NA,NA, 
                                                          4.0, NA,NA,NA,NA,NA, 
                                                          5.2, NA,NA,NA,NA,NA, 
                                                          6.4),
                srtCol = 90, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 3.5,
                sepcolor = "white", colsep = c(7, 13, 19, 25), 
                rowsep = c(7, 13, 19, 25))
      dev.off()
    } else {
      png(paste0(strsplit(path, "info_modules/")[[1]], "sd_row", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
          width = 23, height = 20, units = 'cm', res = 300)
      heatmap.2(modules.sd.row, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
                keysize = 0.5, key = F, margins = c(5,5), 
                col = my.colors2, labRow = NA, labCol = NA,
                srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 2,
                sepcolor = "white", colsep = c(7, 13, 19, 25), 
                rowsep = c(7, 13, 19, 25))
      dev.off()
    }
    
    ## sd columns corrected by features
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_col_correct", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.col/features, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), 
              col = my.colors2, labRow = NA, labCol = NA,
              srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 2,
              sepcolor = "white", colsep = c(7, 13, 19, 25), 
              rowsep = c(7, 13, 19, 25)) 
    dev.off()
    
    ## sd rows corrected by number of individuals
    png(paste0(strsplit(path, "info_modules/")[[1]], "sd_row_correct", "0", substr(cor.limit,3,3), "small", Sys.Date(), ".png"), 
        width = 23, height = 20, units = 'cm', res = 300)
    heatmap.2(modules.sd.row/individuals, Rowv = NA, Colv = NA, trace = "none", dendrogram = "none", 
              keysize = 0.5, key = F, margins = c(5,5), 
              col = my.colors2, labRow = NA, labCol = NA,
              srtCol = 0, sepwidth = c(0.2,0.2), cexRow = 2, cexCol = 2,
              sepcolor = "white", colsep = c(7, 13, 19, 25), 
              rowsep = c(7, 13, 19, 25))
    dev.off()
    
    
  }
}



# threshold combination (multiple) ----------------------------------------


# analysis threshold combination ------------------------------------------
# extract isa scores as metabomatching input
urine.list <- c("urine.nmr.bin.zscored", "urine.nmr.bin.decov", "urine.nmr.focus700.zscored")

for(i in 1:3){
  urine <- eval(parse(text = urine.list[i]))
  rownames(urine) <- c("id", paste("id", urine[2:nrow(urine),1], sep = "_"))
  urine <- urine[,-1]
  colnames(urine) <- paste("ppm", urine[1,], sep ="_")
  urine <- urine[-1,]
  for(cor.limit in c(0.3)){
    path <- paste0("../result/isa/urine/", urine.list[i], "/corr", cor.limit, "/Rdata/")
    rfiles <- list.files(path, pattern = ".Rdata")
    for(file in rfiles){
      load(paste0(path, file))
      if(ncol(isaresults$rows) != 0){
        input <- metabomatching_input(urine, isaresults, 
                                      ppm = as.numeric(sapply(strsplit(colnames(urine), "_"), "[[", 2)))
        dir.create(paste0("../result/metabomatching/", urine.list[i], "/corr",
                          cor.limit, "/ps.", strsplit(file, ".Rdata")[[1]]))
        write.table(input, paste0("../result/metabomatching/", urine.list[i], 
                                  "/corr", cor.limit, "/ps.", strsplit(file, ".Rdata")[[1]],
                                  "/", strsplit(file, ".Rdata")[[1]], ".pseudospectrum.tsv"), 
                    quote = F, sep = "\t", row.names = F)
        
      }
    }
  }
}
            
