# serum <- serum-(rep(rowMeans(serum), each = ncol(serum))/rep(apply(serum, 1, sd, na.rm = TRUE), each = ncol(serum))) #Normalize subjects
# serum <- serum-(rep(colMeans(serum), each = ncol(serum))/rep(apply(serum, 2, sd, na.rm = TRUE), each = ncol(serum))) #Normalize variables
# apply(serum, 2, sd)

# urine <- urine-(rep(rowMeans(urine), each = ncol(urine))/rep(apply(urine, 1, sd, na.rm = TRUE), each = ncol(urine))) #Normalize subjects
# urine <- urine-(rep(colMeans(urine), each = ncol(urine))/rep(apply(urine, 2, sd, na.rm = TRUE), each = ncol(urine))) #Normalize variables


# scratch -----------------------------------------------------------------

isa2image(urine, data.isa = urine.isa, type = "isa", n = 5)
colnames(ppm.corr)[ppm.corr[5,] == max(ppm.corr[5,])]
par(mfrow = c(2,1))
plot(hmdb$ppm[hmdb$HMDB.ID=="HMDB04983"], hmdb$adjusted_height[hmdb$HMDB.ID=="HMDB04983"], type = "h", 
     xlim = c(0,9))
plot(ppm.u, urine["id_996",unlist(lapply(ppm.u, grep, colnames(urine)))], type = "h", xlim = c(0,9))

max(ppm.corr)

urine[5, 1:10]
metabolome["102",1:10]
isa2image(urine, data.isa = urine.isa, type = "isa", n = 1)


# age classes extraction function:
classes <- function(data, intervals, steps, age, sex, sex_variable) {
  classes_ma <- matrix(NA, 0, 3)
  for(interval in intervals){
    if (grepl("<", interval)) {
      classes_ma <- rbind(classes_ma, c(sum(data[,colnames(data) %in% age] < na.omit(as.numeric(unlist(strsplit(interval, "<"))))[1], na.rm = T), 
                                        sum(data[,colnames(data) %in% age] < na.omit(as.numeric(unlist(strsplit(interval, 
                                                                                                                "<"))))[1] & data[,colnames(data) %in% sex] == sex_variable[1], na.rm = T), 
                                        sum(data[,colnames(data) %in% age] < na.omit(as.numeric(unlist(strsplit(interval, 
                                                                                                                "<"))))[1] & data[,colnames(data) %in% sex] == sex_variable[2], na.rm = T)))
    } else if (grepl(">", interval)) {
      classes_ma <- rbind(classes_ma, c(sum(data[,colnames(data) %in% age] > na.omit(as.numeric(unlist(strsplit(interval,
                                                                                                                ">"))))[1], na.rm = T), 
                                        sum(data[,colnames(data) %in% age] > na.omit(as.numeric(unlist(strsplit(interval,
                                                                                                                ">"))))[1] & data[,colnames(data) %in% sex] == sex_variable[1], na.rm = T), 
                                        sum(data[,colnames(data) %in% age] > na.omit(as.numeric(unlist(strsplit(interval,
                                                                                                                ">"))))[1] & data[,colnames(data) %in% sex] == sex_variable[2], na.rm = T)))
    } else {
      classes_ma <- rbind(classes_ma, c(sum(data[,colnames(data) %in% age] %in% seq(as.numeric(unlist(strsplit(interval, 
                                                                                                               "-")))[1], as.numeric(unlist(strsplit(interval, "-")))[2], by = steps)), 
                                        sum(data[,colnames(data) %in% age]  %in% seq(as.numeric(unlist(strsplit(interval, 
                                                                                                                "-")))[1], as.numeric(unlist(strsplit(interval, "-")))[2], 
                                                                                     by = steps) & data[,colnames(data) %in% sex] == sex_variable[1]), 
                                        sum(data[,colnames(data) %in% age]  %in% seq(as.numeric(unlist(strsplit(interval,
                                                                                                                "-")))[1], as.numeric(unlist(strsplit(interval, "-")))[2], 
                                                                                     by = steps) & data[,colnames(data) %in% sex] == sex_variable[2])))
    }
  }
  return(as.data.frame(classes_ma))
}

intervals <- c("<35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", ">75")
classes_age <- classes(data = pheno, intervals = intervals, steps = 0.1, sex = "SEX", age = "AGE", sex_variable = c(0,1))
classes_age <- cbind(classes_age, intervals)
colnames(classes_age) <- c("all", "sex 0", "sex 1", "age")


pheno$age_class <- cut(pheno$AGE, seq(30, 80, by = 5))





pheno$age_class <- cut(pheno$AGE, seq(30, 80, by = 5))

library(ggplot2)
library(plyr) 
install.packages("plotrix")
library(plotrix)

pyramid.plot(as.numeric(classes_age[,2]), as.numeric(classes_age[,3]), labels = intervals)

ggplot(na.omit(pheno), aes(x = age_class, fill = as.factor(SEX))) +
  geom_bar(position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme_bw()

n1 <- ggplot(na.omit(pheno), aes(x = age_class, fill = as.factor(SEX))) + 
  geom_bar() + 
  scale_y_continuous(breaks = seq(-1500, 1500, 500), labels = paste0(as.character(c(seq(15, 0, -5), seq(5, 15, 5))), "m")) +
  coord_flip() + 
  scale_fill_brewer(palette = "Set1") + 
  theme_bw()

n1

ggplot(data=as.data.frame(classes_age)) +
  geom_bar(aes(age,sex0,group=sex,fill=sex), stat = "identity",subset(pheno,pheno$SEX=="0")) +
  geom_bar(aes(age,-sex1,group=sex,fill=sex), stat = "identity",subset(pheno,pheno$SEX=="1")) +
  scale_y_continuous(breaks=seq(-100,40,10),labels=abs(seq(-100,40,10))) +
  coord_flip()


ggplot(classes_age) +
  geom_col(aes(age, classes_age$`sex 0`, fill="#0000CC"), width = 0.5) +
  geom_col(aes(age, -classes_age$`sex 1`), fill = "#0000CC", width = 0.5) +
  scale_fill_manual(values=c("#FF0000", "#0000CC")) +
  coord_flip() + 
  theme_bw()



Z <- t(as.matrix(1:100))

pal.1 <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
pal.2 <- colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=3)
my_color <- colorRampPalette(c("blue", "red"))(n = 10000)

image(Z, col=my_color)
image(Z, col=my_color, las=1, xlab="",ylab="",xaxt="n",yaxt="n")

x=matrix(rnorm(100),nrow=10)*100
xmin=0; xmax=100;
collist<-c("#053061",'#2166AC','#4393C3','#92C5DE','#D1E5F0','#F7F7F7','#FDDBC7','#F4A582',
           '#D6604D','#B2182B','#67001F')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorUsed <- my_color[round(1+(min(module.n)-min(uri.ser))*10000/(max(uri.ser)-min(uri.ser))
                                ) : round( (max(module.n)-min(uri.ser))*10000/(max(uri.ser)-min(uri.ser)) )]

min(module.n)
min(uri.ser)




table(isaColNames(urine, data.isa = urine.isa, type = "isa", n = 149
) %in% isaColNames(urine, data.isa = urine.isa, type = "isa", n = 147))

table(isaRowNames(urine, data.isa = urine.isa, type = "isa", n = 147
) %in% isaRowNames(urine, data.isa = urine.isa, type = "isa", n = 149))


## joint urine serum datasets, plot two modules at the same time
isaRow = uri.ser.isa$rows[, 124] != 0
isaCol = uri.ser.isa$columns[, 124] != 0
module.n <- t(as.matrix(uri.ser[isaRow, isaCol, drop=FALSE]))
module.n <- module.n[unique(c(isaColNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 24), 
                              rownames(module.n))), ]
image(module.n, axes = F, main = "Module 24 and 124", col=ColorUsed)
mtext(text=colnames(module.n), side=2, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = 0.7)
mtext(text=rownames(module.n), side=1, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = 0.5)

isaColNames(uri.ser, data.isa = uri.ser.isa, type = "isa", 
            n = 24) %in% isaColNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 124)

table(isaRowNames(uri.ser, data.isa = uri.ser.isa, type = "isa", 
                  n = 447) %in% isaRowNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 432))
length(isaRowNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 447))
length(isaRowNames(uri.ser, data.isa = uri.ser.isa, type = "isa", n = 432))



## example of ppms
ppm.ser <- as.numeric(sapply(strsplit(isaColNames(serum, data.isa = serum.isa, type = "isa", n = 11), "_"), "[[", 2))
ppm.ur  <- as.numeric(sapply(strsplit(isaColNames(urine, data.isa = urine.isa, type = "isa", n = 1), "_"), "[[", 2))


#correlation matrix representation (heatmap)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
col_breaks = c(seq(0,0.39,length=100),
               seq(0.4,0.79,length=100),
               seq(0.8,1,length=100))


heatmap.2(sqrt(corr.matrix^2), 
          cellnote = round(corr.matrix, 2),
          notecex= 0.7,
          main = "Correlation Matrix",
          notecol="black",
          density.info="none",
          trace="none",       
          dendrogram = "none",
          Colv="NA",
          Rowv = "NA",
          keysize = 0.5,
          key.title = "Identity",
          key.xlab = NA, 
          xlab = "urinary ppm", 
          ylab = "serum ppm", cexRow = 0.7, cexCol = 0.7,
          col = my_palette, breaks = col_breaks)





# metabomatching test using R ---------------------------------------------

module.es <- data.frame(serum.isa$columns[,1])
rownames(module.es) <- colnames(serum)
colnames(module.es) <- "isaScore"
module.es$ppm <- as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2))

# isa score as pseudospectrum -------------------------------------------------

Table_of_Peaks <- read.table(paste0(pathOS, "Mimi/Stage_CBG/1.METABOMATCHING/Files/Metabolite/hmdb.20160809.180.slop"), h = T)

score_final <- matrix(0, 1, 0)
for(j in unique(Table_of_Peaks$HMDB.ID)){
  k=0
  ppm_ps=NULL
  for(i in Table_of_Peaks$ppm[Table_of_Peaks$HMDB.ID==j]){
    ppm  <- module.es$ppm[findInterval(module.es$ppm, 
                                       c((i-0.03), (i+0.03)), rightmost.closed = T)==1]
    ppm_ps <- unique(c(ppm_ps, ppm))
  }
  z  <- module.es[module.es$ppm %in% ppm_ps,][,1]
  if(length(z) == 0){ ## use nrow when multiple column!!
    chisq <- data.frame(matrix(0, 1, 1))
  } else {
    chisq <- data.frame(sum(z^2)) ## when multiple modules use: data.frame(apply(z^2, 2, sum))
  }
  colnames(chisq) <- j
  k=k+length(z) ## use nrow for multiple modules
  score <- -log10(pgamma(q = chisq[[j]], shape = k/2, scale = 2, lower.tail = F))
  score_final <- cbind(score_final, score)  
}

score_final <- t(as.data.frame(score_final))
rownames(score_final) <- unique(Table_of_Peaks$HMDB.ID)
colnames(score_final) <- colnames(linear.regression.z)[-1]


## top 8 score for each SNP

final_score_log <- as.data.frame(score_final)

top8 <- data.frame(top1=0, top2=0, top3=0, top4=0, top5=0, top6=0, top=0, top8=0)

for(i in 1:ncol(final_score_log)){
  n <- unlist(unique(sapply(head(sort(final_score_log[[i]], decreasing = T), 8), grep, final_score_log[[i]])))
  top8[(2*i-1),] <- rownames(final_score_log)[n]
  top8[(2*i),]   <- round(head(sort(final_score_log[[i]], decreasing = T), 8), 4)
}

rownames(top8) <- make.names(rep(colnames(linear.regression.z)[-1], each=2), unique = T)
colnames(top8) <- c("top1", "top2", "top3", "top4", "top5", "top6", "top7", "top8")

write.table(top8, paste0(pathOS, "Mimi/Stage_CBG/1.METABOMATCHING/Files/Results/top_180.20160818.txt", quote = F)
            
            
# visualization -----------------------------------------------------------

ps.plot <- data.frame(ppm = c(Table_of_Peaks[Table_of_Peaks$HMDB.ID == "HMDB00197",]$ppm, 
                              module.es$ppm), h = c(Table_of_Peaks[Table_of_Peaks$HMDB.ID == "HMDB00197",]$adjusted_height,
                                                    -(module.es$isaScore)))

plot(ps.plot$ppm, ps.plot$h, type = "h")                      



