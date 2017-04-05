

# Import package ----------------------------------------------------------
library(isa2)
library(gplots)

# Import Data -------------------------------------------------------------

## serum data
serum <- read.csv("C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/metabolomics/serum.nmr.focus.all.colaus1.20160830.csv", h = F)
ppm <- paste("ppm", serum[1,2:ncol(serum)], sep ="_")
tp  <- paste("id", serum[2:nrow(serum),1], sep = "_")
serum <- serum[2:nrow(serum), 2:ncol(serum)]
colnames(serum) <- ppm
rownames(serum) <- tp

## urine data
urine <- read.csv("C://Mimi/Stage_CBG/2.EXPRESSION_MODULE/Data/metabolomics/urine.nmr.focus.all.colaus1.20161205.csv", h = F)
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

## Data Normalization
# serum
serum[serum<1]=1 # to avoid negative numbers
serum <- log10(serum)# log-transformed data
serum <- serum-(rep(rowMeans(serum), each = ncol(serum))/rep(apply(serum, 1, sd, na.rm = TRUE), each = ncol(serum))) #Normalize subjects
serum <- serum-(rep(colMeans(serum), each = ncol(serum))/rep(apply(serum, 2, sd, na.rm = TRUE), each = ncol(serum))) #Normalize variables
#urine
urine[urine<1]=1 # to avoid negative numbers
urine <- log10(urine)# log-transformed data
urine <- urine-(rep(rowMeans(urine), each = ncol(urine))/rep(apply(urine, 1, sd, na.rm = TRUE), each = ncol(urine))) #Normalize subjects
urine <- urine-(rep(colMeans(urine), each = ncol(urine))/rep(apply(urine, 2, sd, na.rm = TRUE), each = ncol(urine))) #Normalize variables

# isa run ---------------------------------------------------------

## serum
serum.isa <- isa(as.matrix(serum), thr.row = seq(2, 3.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))

## urine
urine.isa <- isa(as.matrix(urine), thr.row = seq(2, 3.5, by = 0.1), thr.col = seq(2, 3.5, by = 0.1))

## urine and serum combination
uri.ser <- cbind(urine, serum)
uri.ser.isa <- isa(as.matrix(uri.ser))

## module info
serum.modules <- isaModules(serum.isa, type = "isa")
urine.modules <- isaModules(urine.isa, type = "isa")
hist(serum.modules$colGroups, col = "red", main = "Peaks by modules", xlab = "serum")
hist(serum.modules$rowGroups, col = "orange", main = "Individuals by modules", xlab = "serum")
hist(urine.modules$colGroups, col = "blue", main = "Peaks by modules", xlab = "urine")
hist(urine.modules$rowGroups, col = "yellow", main = "Individuals by modules", xlab = "urine")


## look for similarity of modules between urine and serum
identity <- matrix(NA, nrow = nrow(urine.modules), ncol = nrow(serum.modules))
for(i in 1:nrow(urine.modules)){
  id.u <- isaRowNames(data = urine, type = "isa", data.isa = urine.isa, n = i)
  for(j in 1:nrow(serum.modules)){
    id.s <- isaRowNames(data = serum, type = "isa", data.isa = serum.isa, n = j)
    perc <- round(ifelse(length(id.u) > length(id.s), sum(id.u %in% id.s)/length(id.u), sum(id.s %in% id.u)/length(id.s)), 2)
    identity[i,j] <- perc
  }
}
rownames(identity) <- paste0("mod.u", seq(1,nrow(urine.modules)))
colnames(identity) <- paste0("mod.s", seq(1,nrow(serum.modules)))

identity.sel <- identity[apply(identity, MARGIN = 1, function(x) any(x >= 0.8)), ]
identity.sel <- identity[apply(identity, MARGIN = 2, function(x) any(x >= 0.8)), ]

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
col_breaks = c(seq(0,0.39,length=100),               # for red
               seq(0.4,0.79,length=100),           # for yellow
               seq(0.8,1,length=100))             # for green

heatmap.2(identity, 
          cellnote = identity,
          notecex= 0.7,
          main = "Matching individuals in serum and urine modules",
          notecol="black",
          density.info="none",
          trace="none",       
          col=my_palette,
          dendrogram = "none",
          Colv="NA",
          Rowv = "NA",
          keysize = 0.5,
          key.title = "Identity",
          key.xlab = NA,
          breaks = col_breaks)

## plot modules (example of two modules with same individuals)
par(mfrow = c(1,2))
isa2image(data = serum, type = "isa", data.isa = serum.isa, n = 8)
isa2image(data = urine, type = "isa", data.isa = urine.isa, n = 9)
isa2image(data = serum, type = "isa", data.isa = serum.isa, n = 4)
isa2image(data = urine, type = "isa", data.isa = urine.isa, n = 6)

# PPA ---------------------------------------------------------------------
ser_ur <- list(as.matrix(t(serum)), as.matrix(t(urine)))
metabo.ppa.def <- ppa(ser_ur)
metabo.ppa <- ppa(ser_ur, thr.row1 = seq(2, 3.5, by = 0.1), thr.row2 = seq(2, 3.5, by = 0.1), thr.col = seq(2,3.5, by = 0.1))

metabo.module.info <- isaModules(data = metabo.ppa, type = "ppa")

hist(metabo.module.info$colGroups)
hist(metabo.module.info2$colGroups)
hist(metabo.module.info$row1Groups)

## ppa plot

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa,
          n = 29, name1 = "serum", name2 = "urine")

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa2,
          n = 44, name1 = "serum", name2 = "urine")


# Cross correlation -------------------------------------------------------

hmdb <- read.table("C://Mimi/Stage_CBG/1.METABOMATCHING/Files/Metabolite/Table_of_peaks180.txt",h= T)
cross.corr <- matrix(NA, nrow = nrow(urine.modules), ncol = length(unique(hmdb$HMDB_ID)))
for(i in 1:nrow(urine.modules)){
  ppm.u <- as.numeric(sapply(strsplit(isaColNames(data = urine, type = "isa", data.isa = urine.isa, n = i),
                                    "_"), "[[", 2))
  for(id in unique(hmdb$HMDB_ID)){
    ppm <- hmdb$ppm[hmdb$HMDB_ID == id]
    perc <- round(ifelse(length(id.u) > length(id.s), sum(id.u %in% id.s)/length(id.u), sum(id.s %in% id.u)/length(id.s)), 2)
    identity[i,j] <- perc
  }
}



Table_of_Peaks <- read.table("C:/Mimi/Stage_CBG/1.METABOMATCHING/Files/Metabolite/hmdb.20160809.180.slop", h=T)
snp <- colnames(read.table("C:/Mimi/Stage_CBG/1.METABOMATCHING/Data/Genotype/hit.snp.genotypes.colaus.urine.20160523.txt", h = T, sep = ","))[-1]
rico_snp <- c('rs4327428','rs37369','rs17169536','rs4488133','rs10774021','rs7314056','rs676882','rs6510300','rs17273533')
synonym <- as.data.frame(fread("C:/Mimi/Stage_CBG/1.METABOMATCHING/Files/Metabolite/synonym20160722.txt", h=T))
synonym <- synonym[!duplicated(synonym$HMDB.ID),]

# Calculate score with mirjam's pseudospectrum ------------------------------

setwd("C:/Mimi/Stage_CBG/Files/Results/score/with.mirjam.ps/")
for(g in snp){
  ps <- read.table(paste("C:/Mimi/Stage_CBG/1.METABOMATCHING/Files/Results/Correlation_results/pseudospectrum/mirjam.", 
                         ".pseudospectrum.tsv", sep = g), h = T)
  score_final <- NULL
  for(id in unique(Table_of_Peaks$HMDB.ID)){
    k=0
    ppm_ps=NULL
    for(i in Table_of_Peaks$ppm[Table_of_Peaks$HMDB.ID==as.character(id)]){
      ppm  <- ppm.u[findInterval((ppm.u), c((i-0.03), (i+0.03)), rightmost.closed = T)==1]
      ppm_ps <- unique(c(ppm_ps, ppm))
    }
    z  <- (ps$beta/ps$se)[ps$shift %in% ppm_ps]
    if(length(z) == 0){
      chisq <- 0
    } else {
      chisq <- sum(z^2)
    }
    k=k+length(z)
    score <- -log10(pgamma(q = chisq, shape = k/2, scale = 2, lower.tail = F))
    score_final <- cbind(score_final, score)  
  }
  score_final <- as.data.frame(t(score_final))
  rownames(score_final) <- unique(Table_of_Peaks$HMDB.ID)
  final_score_log <- merge(data.frame(CAS.ID = unique(Table_of_Peaks$CAS.Number), HMDB.ID = unique(Table_of_Peaks$HMDB.ID), idDat = rep(1, nrow(score_final)), V1 = round(score_final, 4)), synonym)
  final_score_log2 <- final_score_log[order(final_score_log$V1, decreasing = T),]
  write.table(final_score_log2, file = paste("mirjam", g, "score.tsv", sep = "."), sep = "\t", row.names = FALSE, col.names = F, quote = F)
}

# scratch -----------------------------------------------------------------



