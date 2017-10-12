## use isa score as pseudospectrum input to estimate which 
## metabolite is more likely to be representative of a module
## 






## multiply the gene vector with the whole matrix to have the 
## correspondent score for each idividuals! 
## Like doing an extra iteration without thresholding
## 

module.n # exemplary module (#5) issued from the isa run using serum data.
score_features <- as.matrix(data.isa$columns[, 5])

## multiplication of the normalized matrix with gene vector!
score_id <- s.isa.norm$Ec%*%score_features

plot(score_id, data.isa$rows[,5])

# additiopnal iteration ---------------------------------------------------

score_non_norm <- last.iteration(as.matrix(serum_log), data.isa2, Col = T)
score_norm <- last.iteration(serum, data.isa, Col = T)

metabomatching_input <- data.frame(shift = as.numeric(sapply(strsplit(colnames(serum), "_"), "[[", 2)))
for (i in 1:ncol(score)){
  metabomatching_input <- cbind(metabomatching_input, 
                                score[,i],
                                rep(1, ncol(serum)),
                                pt((score[,i]), df = 1, lower.tail = F))
  colnames(metabomatching_input)[colnames(metabomatching_input) %in% c("score[, i]", "rep(1, ncol(serum))",
                                                                       "pt((score[, i]), df = 1, lower.tail = F)")] <- c(paste0("beta/m", i), 
                                                                                                                         paste0("se/m", i),
                                                                                                                         paste0("p/m", i))
}
# write.table(metabomatching_input, file = "../result/metabomatching/serum//module_non_normalized.tsv", 
# quote = F, sep = "\t", row.names = F)

## checking normality for each the p-value distribution
for(i in 1:ncol(score_non_norm)){
  png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/score/non_norm/m",
                    i, "_qqplot", Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
  qqnorm(score_non_norm[,i], main = paste0("Non-normalized m", i, " Normal Q-Q Plot"))
  qqline(score_non_norm[,i], col = 2, lwd = 2)
  text(-2, max(score_non_norm[,i], na.rm = T)-(max(score_non_norm[,i], na.rm = T)/10), 
       paste0("Shapiro-Wilk normality test \n p = ", shapiro.test(score_non_norm[,i][1:5000])$p))
  dev.off()
}
