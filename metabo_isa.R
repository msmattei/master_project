

# Import package ----------------------------------------------------------
library(isa2)

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
metabo.ppa <- ppa(ser_ur)
metabo.ppa2 <- ppa(ser_ur, thr.row1 = seq(1,6, by = 0.5), thr.row2 = seq(1,6, by = 0.5), thr.col = seq(0,3, by = 0.5))

metabo.module.info <- isaModules(data = metabo.ppa, type = "ppa")
metabo.module.info2 <- isaModules(data = metabo.ppa2, type = "ppa")
hist(metabo.module.info$colGroups)
hist(metabo.module.info2$colGroups)
hist(metabo.module.info$row1Groups)

## ppa plot

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa,
          n = 29, name1 = "serum", name2 = "urine")

isa2image(data = serum, data2 = urine, type = "ppa", data.isa = metabo.ppa2,
          n = 44, name1 = "serum", name2 = "urine")


# scratch -----------------------------------------------------------------



