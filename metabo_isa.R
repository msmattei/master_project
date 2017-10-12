
# Specify OS system: ------------------------------------------------------

### for windows
pathOS <- "C://"
### for linux
pathOS <- "/media/mirjam/OS/"

# Import package ----------------------------------------------------------
library(isa2)
library(gplots)
library(ggplot2)


# ISA (and other) functions  ----------------------------------------------
source(paste0(pathOS, 'Mimi/Stage_CBG/2.EXPRESSION_MODULE/expression_module/ISA_functions.R'))
source(paste0(pathOS, 'Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project/functions.R'))

# Grafical parameters -----------------------------------------------------
## green, yellow and red scale used for heatmap plot (used in the case of a graphical representation of the percentage of identity for example)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
col_breaks = c(seq(0,0.33,length=100),             # for red
               seq(0.34,0.66,length=100),          # for yellow
               seq(0.67,1,length=100))             # for green


# Import Data -------------------------------------------------------------

## serum data
serum <- read.csv(paste0(pathOS, "Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/metabolomics/serum.nmr.focus.all.colaus1.20160830.csv"), h = F)
ppm <- paste("ppm", serum[1,2:ncol(serum)], sep ="_")
tp  <- paste("id", serum[2:nrow(serum),1], sep = "_")
serum <- serum[2:nrow(serum), 2:ncol(serum)]
colnames(serum) <- ppm
rownames(serum) <- tp

## urine data
urine <- read.csv(paste0(pathOS, "/Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/metabolomics/urine.nmr.focus.all.colaus1.20161205.csv"), h = F)
ppm <- paste("ppm", urine[1,2:ncol(urine)], sep ="_")
tp  <- paste("id", urine[2:nrow(urine),1], sep = "_")
urine <- urine[2:nrow(urine), 2:ncol(urine)]
colnames(urine) <- ppm
rownames(urine) <- tp

rm("ppm", "tp")
## select matching id (790 individuals having serum and urine metabolomics)
serum <- serum[rownames(serum) %in% rownames(urine),]
serum <- serum[order(rownames(serum)),]
urine <- urine[rownames(urine) %in% rownames(serum),]
urine <- urine[order(rownames(urine)),]


# Data Normalization ------------------------------------------------------

## Data Normalization (z-score normalization)
# serum
serum[serum<1]=1 # to avoid negative numbers
serum_log <- log10(serum) # log-transformed data
## z-score normalization for features
serum_feature_scaled <- scale(serum_log)
## z-score normalization for individulas
serum_ind_scaled <- t(scale(t(serum_log)))
## z-score normalization first for individuals and second for the features
serum_ind_feat_scaled <- scale(serum_ind_scaled)
## z-score normaizazion first accoring to features and second for individuals
serum_feat_ind_scaled <- t(scale(t(serum_feature_scaled)))

### urine
urine[urine<1]=1 # to avoid negative numbers
urine_log <- log10(urine)# log-transformed data
## z-score normalization for features
urine_feature_scaled <- scale(urine_log)
## z-score normalization for individulas
urine_ind_scaled <- t(scale(t(urine_log)))
## z-score normalization first for individuals and second for the features
urine_ind_feat_scaled <- scale(urine_ind_scaled)
## z-score normaizazion first accoring to features and second for individuals
urine_feat_ind_scaled <- t(scale(t(urine_feature_scaled)))

### used data -> first normalized by individuals, second by features
serum <- serum_ind_feat_scaled
urine <- urine_ind_feat_scaled

#### Phenotype data:
pheno_raw <- read.csv(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.raw.colaus1.20161116.csv", h = F, sep = ",", stringsAsFactors = T))
pheno_transf <- read.csv(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.transformed.colaus1.20161116.csv", h = F, sep = ","))
pheno_names <- read.csv(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/trait_names.raw.colaus1.20161116.csv", h = F, sep = ","))
pheno <- pheno_raw
colnames(pheno) <- c("ID", as.character(pheno_names$V1))
rm(pheno_names, pheno_raw)

save(urine, serum, pheno, file = paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/metabo_pheno_data.Rdata"))

# Look to the data!! ------------------------------------------------------
# Know the story around the data!
# Ask concrete questions!
# Always look at the data!
# Transformation if needed!
load(file = "../data/metabo_pheno_data.Rdata")
str(urine)

## Data summary and viszualization
dim(urine)
dim(serum)

table(urine<0)

plot(sort(as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2)), 
          decreasing = T), rev(serum[1,]), type = "h")
plot(sort(as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2)), 
          decreasing = T), rev(serum[2,]), type = "h")

summary(urine$ppm_9.28)

summary(urine)

############ visualization #################
## raw data
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/raw_serum",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum), main = "Raw data (serum)")
dev.off()
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/raw_urine",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine), main = "Raw data (urine)")
dev.off()

## log transformed serum data image
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/log_serum",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_log), main = "Log transformed data (serum)")
dev.off()

## z-score normalization for features (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_scaled_serum",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_feature_scaled), main = "z-transformation for features (serum)")
dev.off()

## z-score normalization for individulas (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_scaled_serum",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_ind_scaled), main = "z-transformation for individuals (serum)")
dev.off()

## z-score normalization first for individuals and second for the features (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_feat_scaled_serum",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_ind_feat_scaled), main = "z-transformation for individuals and features (serum)")
dev.off()
## z-score normalization first for features and second for the individuals (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_ind_scaled_serum",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(serum_feat_ind_scaled), main = "z-transformation for features and individuals (serum)")
dev.off()



### urine
## log transformed urine data image
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/log_urine",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_log), main = "Log transformed data (urine)")
dev.off()

## z-score normalization for features (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_scaled_urine",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_feature_scaled), main = "z-transformation for features (urine)")
dev.off()

## z-score normalization for individulas (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_scaled_urine",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_ind_scaled), main = "z-transformation for individuals (urine)")
dev.off()

## z-score normalization first for individuals and second for the features (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/ind_feat_scaled_urine",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_ind_feat_scaled), main = "z-transformation for individuals and features (urine)")
dev.off()
## z-score normalization first for features and second for the individuals (with image)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/feat_ind_scaled_urine",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(as.matrix(urine_ind_feat_scaled), main = "z-transformation for features and individuals (urine)")
dev.off()


## Correlation between features
# serum
fluid = "serum"
# urine
fluid = "urine"
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_raw",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = fluid))), main = paste0("Features correlation of raw ", fluid, " data"))
dev.off()
# z-score normalized data by features
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_feat_norm",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_feature_scaled")))), main = paste0("Features correlation of features scaled ", fluid, " data"))
dev.off()
# z-score normalized data by individuals
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_ind_norm",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_ind_scaled")))), main = paste0("Features correlation of individuals scaled ", fluid, " data"))
dev.off()
# z-score normalized data first by individuals second by features
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_ind_feat_norm",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_ind_feat_scaled")))), main = paste0("Features correlation of individuals and features scaled ", fluid, " data"))
dev.off()
# z-score normalized data first by features second by indiviuals
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/data_info/", fluid, "_corr_feat_ind_norm",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
image(cor(eval(parse(text = paste0(fluid, "_feat_ind_scaled")))), main = paste0("Features correlation of features and individuals scaled ", fluid, " data"))
dev.off()






# isa run ---------------------------------------------------------

# 27 07 2017 parameters
## serum
serum.isa <- isa(as.matrix(serum), thr.row = seq(2, 3.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))
## urine
urine.isa <- isa(as.matrix(urine), thr.row = seq(2, 3.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))

## save analysis:
# save(serum.isa, urine.isa, file = paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_modules20170727.RData")
# load(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_modules20170502.RData")
# load(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_modules20170727.RData")


# urine and serum combination
uri.ser <- cbind(urine, serum)
colnames(uri.ser) <- c(paste0("u.", colnames(urine)), paste0("s.", colnames(serum)))
uri.ser.isa <- isa(as.matrix(uri.ser))



# isa unique --------------------------------------------------------------
## Increase the number of subjects per module, reduce the total number of modules (to have less, bigger modules)
serum.norm <- isa.normalize(as.matrix(serum_ind_feat_scaled))
urine.norm <- isa.normalize(as.matrix(urine_ind_feat_scaled))

## August 2017, result saved in ser_ur_isa20170827.RData (created using the server)
#serum
serum.isa <- isa(as.matrix(serum_ind_feat_scaled), thr.row = seq(1, 2.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))
## urine
urine.isa <- isa(as.matrix(urine_ind_feat_scaled), thr.row = seq(1, 2.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))

load(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/ser_ur_isa20170827.RData"))

serum.uni <- isa.unique(serum.norm, serum.isa, method = c("cor"), 
                        ignore.div = TRUE, cor.limit = 0.7, neg.cor = TRUE, drop.zero = TRUE)

urine.uni <- isa.unique(urine.norm, urine.isa, method = c("cor"), ignore.div = TRUE, 
                        cor.limit = 0.7, neg.cor = TRUE, drop.zero = TRUE)



# Module analysis ---------------------------------------------------------

## module info
serum.modules   <- isaModules(serum.isa, type = "isa")
urine.modules   <- isaModules(urine.isa, type = "isa")
# uri.ser.modules <- isaModules(uri.ser.isa)
serum.modules.uni <- isaModules(serum.uni, type = "isa")
urine.modules.uni <- isaModules(urine.uni, type = "isa")

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


head(serum.modules)
## serum -> columns
## urine -> rows of identity.urine.serum

## change row and columns names the same as serum.modules (or urine.modules)
colnames(identity.urine.serum) <- gsub(pattern = "mod", replacement = "Module", colnames(identity.urine.serum))
rownames(identity.urine.serum) <- gsub(pattern = "mod", replacement = "Module", rownames(identity.urine.serum))

## correct the identity value of two modules by the number of individuals that are present in the modules
id_es <- matrix(NA, nrow(identity.urine.serum), ncol(identity.urine.serum))
for(j in 1:ncol(identity.urine.serum)){
  n = serum.modules[colnames(identity.urine.serum)[j],][2]
  for(i in 1:nrow(identity.urine.serum)){
    m = urine.modules[rownames(identity.urine.serum)[i],][2]
    id_es[i,j] <- identity.urine.serum[i,j]/(m+n)[,1]
  }
}

colnames(id_es) <- colnames(identity.urine.serum)
rownames(id_es) <- rownames(identity.urine.serum)

heatmap.2(round(id_es/(max(id_es)), 2), cellnote = round(id_es/(max(id_es)), 2), notecex= 0.7,
          main = "Matching individuals in serum and urine modules", margins = c(7,7), xlab = "Serum", ylab = "Urine",
          notecol="black", density.info="none", trace="none", col=my_palette,
          dendrogram = "none", Colv="NA", Rowv = "NA", keysize = 0.5, cexRow = 1.2, cexCol = 1.2,
          key.title = "Identity", key.xlab = NA, breaks = col_breaks)
par(mfrow = c(1,2))
isa2image(serum, data.isa = serum.isa, type = "isa", n = 197)
isa2image(urine, data.isa = urine.isa, type = "isa", n = 30)


# example of similar modules
## plot modules (example of two modules with similar individuals)
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_es",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,2))
isa2image(data = serum, type = "isa", data.isa = serum.isa, n = 11, cex = 1)
isa2image(data = urine, type = "isa", data.isa = urine.isa, n = 1, cex = 1)
dev.off()

png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_es2_",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
par(mfrow = c(1,2))
isa2image(data = serum, type = "isa", data.isa = serum.isa, n = 115, cex = 1)
isa2image(data = urine, type = "isa", data.isa = urine.isa, n = 257, cex = 1)
dev.off()

# Combined urine and serum dataset
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_uri_ser_joint",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.uri.ser,
          main = "Identity between modules \n created using the combined urine and serum dataset",
          notecol="black", density.info="none", trace="none", col=my_palette, dendrogram = "none",
          Colv="NA", Rowv = "NA", keysize = 0.5, key.title = "Identity", key.xlab = NA, breaks = col_breaks,
          labRow = FALSE, labCol = FALSE)

dev.off()
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/module_similarity_uri_ser_joint_sel",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
heatmap.2(identity.uri.ser.sel, cellnote = identity.uri.ser.sel, notecex= 0.5,
          main = "Identity between modules \n created using the combined urine and serum dataset",
          notecol="black", density.info="none", trace="none", col=my_palette, dendrogram = "none",
          Colv="NA", Rowv = "NA", keysize = 0.5, key.title = "Identity", key.xlab = NA, breaks = col_breaks)

dev.off()

par(mfrow = c(1,2))
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 124)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 24)

isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 43)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 181)

isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 83)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 64)

isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 224)
isa2image(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 193)

isaColNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 193)

-------------

# idea: identity matrix between individuals and ppms -> one minus the other -> max of absolute value
# more differences between concordance module individuals or ppms!

-------------

### correlation matrix between metabolic feature in modules created using urine and serum 
# databases having similar individuals
## module 11 in serum and module 1 in urine (example july 2017)
isaRow = serum.isa$rows[, 115] != 0
isaCol = serum.isa$columns[, 115] != 0
module.ser <- t(as.matrix(serum[isaRow, isaCol, drop=FALSE]))

isaRow = urine.isa$rows[, 257] != 0
isaCol = urine.isa$columns[, 257] != 0
module.ur <- t(as.matrix(urine[isaRow, isaCol, drop=FALSE]))

mod.ur <- module.ur[,colnames(module.ur) %in% colnames(module.ser)]
mod.ser <- module.ser[,colnames(module.ser) %in% colnames(module.ur)]

i=4
plot(as.numeric(sapply(strsplit(rownames(mod.ur), "_"), "[[", 2)), mod.ur[,i], pch = 19, col = "yellow")
points(as.numeric(sapply(strsplit(rownames(mod.ser), "_"), "[[", 2)), mod.ser[,i], pch = 19, col = "red")


ppm.u <- as.numeric(sapply(strsplit(rownames(mod.ur), "_"), "[[", 2))
ppm.s <- as.numeric(sapply(strsplit(rownames(mod.ser), "_"), "[[", 2))

i = 1
plot(mod.ser[,i], mod.ur[,i][sapply(ppm.s, function(x) which.min(abs(ppm.u-x)))])
sapply(ppm.u, function(x) which.min(abs(ppm.s-x)))

# plot(sort(mod.ur[,1]), sort(mod.ser[,1])) no sense!

intensity_correlation <- matrix(NA, nrow(mod.ser), nrow(mod.ur))
for(i in 1:nrow(mod.ser)){
  for(j in 1:nrow(mod.ur)){
    r <- cor(mod.ser[i,], mod.ur[j,], method = "pearson")
    intensity_correlation[i,j] <- round(r,2)
  }
}
rownames(intensity_correlation) <- sapply(strsplit(rownames(mod.ser), "_"), "[[", 2)
colnames(intensity_correlation) <- sapply(strsplit(rownames(mod.ur), "_"), "[[", 2)

png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/corr_matrix",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
 
heatmap.2(intensity_correlation, cellnote = intensity_correlation, notecex= 0.8, margins = c(7,7),
          main = "Intensity correlation", cexRow = 1.2, cexCol = 1.2, xlab = "ppm Urine", ylab = "ppm Serum",
          notecol="black", density.info="none", trace="none", dendrogram = "none",
          Colv="NA", Rowv = "NA", keysize = 0.7, key.title = "Identity", key.xlab = NA)
  dev.off()
png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/corr_plot",
           Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
plot(mod.ser["ppm_2.3066",], mod.ur[ "ppm_1.0186",], pch = 19)
abline(lm(mod.ur[12,]~mod.ser[7,]), col="red", lwd = 2)
dev.off()
hist(urine[,"ppm_1.0186"],20) ## to have an idea about the intensity distribution of this spectral position
hist(serum[,"ppm_2.3066"], 20)

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


# Cross correlation -------------------------------------------------------
ccf(x = as.matrix(data.frame(ppm = es.ppm, height = es.score)), 
    y = as.numeric(as.matrix(data.frame(ppm = hmdb$ppm[hmdb$HMDB.ID == "HMDB00157"], 
                                        height = hmdb$adjusted_height[hmdb$HMDB.ID == "HMDB00157"]))), type = "correlation")


hmdb <- read.table(paste0(pathOS, "/Mimi/Stage_CBG/1.METABOMATCHING/Files/Metabolite/hmdb.20160809.180.slop", h = T))

ppm.corr <- matrix(NA, nrow = nrow(urine.modules), ncol = length(unique(hmdb$HMDB.ID)))
for(i in 1:nrow(urine.modules)){
  ppm.u <- as.numeric(sapply(strsplit(isaColNames(data = urine, type = "isa", data.isa = urine.isa, 
                                                  n = i), "_"), "[[", 2))
  for(j in 1:length(unique(hmdb$HMDB.ID))){
    ppm.hmdb <- hmdb$ppm[hmdb$HMDB.ID  == unique(hmdb$HMDB.ID)[j]]
    ppm.tot=NULL
    for(ppm.i in hmdb$ppm[hmdb$HMDB.ID==unique(hmdb$HMDB.ID)[j]]){
      ppm  <- ppm.u[findInterval((ppm.u), c((ppm.i-0.03), (ppm.i+0.03)), rightmost.closed = T)==1]
      ppm.tot <- unique(c(ppm.tot, ppm))
    }
    perc <- length(ppm.tot)/(length(ppm.u)+length(ppm.hmdb))
    ppm.corr[i,j] <- perc
  }
}
colnames(ppm.corr) <- unique(hmdb$HMDB.ID)
rownames(ppm.corr) <- rownames(urine.modules)

ppm.max <- data.frame(module = paste0("Module", apply(ppm.corr, MARGIN = 2, FUN = which.max)), 
                      perc = round(apply(ppm.corr, MARGIN = 2, FUN = max), 2))

### "pseudospectrum" image using isa score
es <- isaColNames(urine, data.isa = urine.isa, type ="isa", n=252)
es.ppm <- as.numeric(sapply(strsplit(isaColNames(urine, data.isa = urine.isa, type = "isa", n = 252), "_"), "[[", 2))
es.score <- urine.isa$columns[colnames(urine) %in% es, 252]
colnames(urine)[138]
match(es, colnames(urine))


### "pseudospectrum" image using isa score
es <- isaColNames(serum, data.isa = serum.isa, type ="isa", n=250)
serum.isa$columns[colnames(serum) %in% es, 250]
colnames(serum)[138]
match(es, colnames(serum))

plot(serum.isa$columns[colnames(serum) %in% es, 250], type = "h")




plot(c(es.ppm, hmdb$ppm[hmdb$HMDB.ID == "HMDB00157"]), c(es.score, -hmdb$adjusted_height[hmdb$HMDB.ID == "HMDB00157"]), type = "h")
plot(hmdb$ppm[hmdb$HMDB.ID == "HMDB00157"], hmdb$adjusted_height[hmdb$HMDB.ID == "HMDB00157"], type = "h")



