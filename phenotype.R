


# Phenotype ---------------------------------------------------------------
pheno           <- read.csv(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.raw.colaus1.20161116.csv"), h = F, sep = ",", stringsAsFactors = T)
pheno_transf    <- read.csv(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/traits.transformed.colaus1.20161116.csv"), h = F, sep = ",")
pheno_names     <- read.csv(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/data/trait_names.raw.colaus1.20161116.csv"), h = F, sep = ",")
colnames(pheno) <- c("ID", as.character(pheno_names$V1))
colnames(pheno_transf) <- c("ID", as.character(pheno_names$V1))
pheno$ID        <- paste0("id", pheno$ID)
pheno_transf$ID <- paste0("id_", pheno_transf$ID)
rm(pheno_names)

## phenotype information
# classes extraction function:
classes <- function(data, from, to, by, variable_name, variable_name2) {
  data_df <- data.frame(data)
  data_df$variable_class <- cut(as.numeric(unlist(data[colnames(data) %in% variable_name])), seq(from, to, by))
  return(as.data.frame(table(data_df[,colnames(data_df) %in% c("variable_class", variable_name2)])))
}

# Gender distribution with age
age_classes <- classes(data = pheno, from = 30, to = 80, by = 5, variable_name = "AGE", variable_name2 = "SEX")

# Pyramidal plot
ggplot() +
  geom_col(aes(age_classes$variable_class[age_classes$SEX == 0], 
               age_classes$Freq[age_classes$SEX == 0], fill = "Sex 0", width = 0.5)) +
  geom_col(aes(age_classes$variable_class[age_classes$SEX == 1], 
               -age_classes$Freq[age_classes$SEX == 1], fill = "Sex 1", width = 0.5)) +
  scale_fill_manual(values=c("#FF0000", "#0000CC")) +
  labs(title="Age and Gender distribution", y ="Number of subjects", x = "Age classes", fill = "Gender") +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=16, face="bold"),axis.title=element_text(size=16, face="bold"))

pheno$age_class <- cut(pheno$AGE, seq(30, 80, by = 5))
boxplot(pheno$GLUC~pheno$age_class)
boxplot(pheno$SBP~pheno$age_class)
boxplot(pheno$BMI~pheno$age_class)
boxplot(pheno$TRIG~pheno$age_class)
boxplot(pheno$HDLCH~pheno$age_class)
boxplot(pheno$CHOL~pheno$age_class)
boxplot(pheno$LDLCH~pheno$age_class)
boxplot(pheno$ADTRN~pheno$age_class)

boxplot(pheno$GLUC~pheno$SEX)
boxplot(pheno$SBP~pheno$SEX)
boxplot(pheno$BMI~pheno$SEX)
boxplot(pheno$TRIG~pheno$SEX)
boxplot(pheno$HDLCH~pheno$SEX)
boxplot(pheno$CHOL~pheno$SEX)
boxplot(pheno$LDLCH~pheno$SEX)
boxplot(pheno$ADTRN~pheno$SEX)

plot(pheno$AGE, pheno$GLUC)

model <- lm(CHOL ~ SEX + AGE + BMI + PHYACT, data = pheno)
model2 <- glm(CHOL ~ SEX + AGE + BMI + PHYACT, data = pheno)
summary(lm(CHOL ~ as.factor(SEX) + AGE + BMI + as.factor(PHYACT) + as.factor(SMK), data = pheno))
summary(lm(SBP ~ SEX + AGE + BMI + PHYACT + SMK, data = pheno))
summary(lm(TRIG ~ SEX + AGE + BMI + PHYACT, data = pheno))
summary(lm(HDLCH ~ SEX + AGE + BMI + PHYACT, data = pheno))
summary(lm(LDLCH ~ SEX + AGE + BMI + PHYACT, data = pheno))
summary(lm(ADTRN ~ SEX + AGE + BMI + PHYACT, data = pheno))

summary(lm(CHOL ~ BMI, data = pheno))

## Looking if some module are related to any kind of phenotype
## selecti individuals present in the metabolomics data
phen <- pheno
phen <- phen[phen$ID %in% rownames(urine),]
phen <- phen[order(phen$ID),]
phen <- unique(phen)
attach(phen)

phen <- pheno_transf[pheno_transf$ID %in% rownames(urine),]
phen <- unique(phen)

### Variable distribution
## Look at the densities distribution, at the normal qqplot and perform the shapiro-wilk test to check normality.
covariables <- c("GLUC", "TRIG", "HDLCH", "CHOL", "LDLCH", "SBP")

for(trait in covariables){
  png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
                    trait, Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
  plot(density(na.omit(phen[,trait])), main = paste0(trait, " Density Plot"))
  dev.off()
  png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
                    trait, "qqplot", Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
  qqnorm(phen[,trait], main = paste0(trait, " Normal Q-Q Plot"))
  qqline(phen[,trait], col = 2)
  text(-2, max(phen[,trait], na.rm = T)-(max(pheno_transf[,trait], na.rm = T)/10), 
       paste0("Shapiro-Wilk normality test \n p = ", shapiro.test(phen[,trait])$p))
  dev.off()
}

## transformed data:
for(trait in covariables[c(1,2,6)]){
  png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
                    trait, "tranf_", Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
  plot(density(na.omit(pheno_transf[,trait])), main = paste0("Transformed ", trait, " Density Plot"))
  dev.off()
  png(paste0(paste0(pathOS, "/Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/result/figure/normality_test/",
                    trait, "tranf_", "qqplot", Sys.Date(), ".png")), width = 23, height = 20, units = 'cm', res = 300)
  qqnorm(pheno_transf[,trait], main = paste0("Transformed ", trait, " Normal Q-Q Plot"))
  qqline(pheno_transf[,trait], col = 2)
  text(-2, max(pheno_transf[,trait], na.rm = T)-(max(pheno_transf[,trait], na.rm = T)/10), 
       paste0("Shapiro-Wilk normality test \n p = ", shapiro.test(pheno_transf[,trait][1:5000])$p))
  dev.off()
}

# serum:
id.ser <- isaRowNames(serum, data.isa = serum.isa, type = "isa", n = 1)
id.ur  <- isaRowNames(urine, data.isa = urine.isa, type = "isa", n = 1)

age.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "AGE")
sex.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "SEX")
chol.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "CHOL")
GLUC.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "GLUC")
SBP.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "SBP")
BMI.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "BMI")
TRIG.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "TRIG")
HDLCH.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "HDLCH")
LDLCH.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "LDLCH")
SMK.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "SMK")
PHYACT.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "PHYACT")
ADTRN.test <- t.test.variable(phenotype.data = phen, metabo.data = serum, metabo.isa = serum.isa, variable1 = "ADTRN")

p_val_result <- data.frame(p_val.age = age.test$p_val, p_val.sex = sex.test$p_val, p_val.chol = chol.test$p_val, 
                           p_val.GLUC = GLUC.test$p_val, p_val.SBP = SBP.test$p_val, p_val.BMI = BMI.test$p_val, 
                           p_val.TRIG = TRIG.test$p_val, p_val.HDLCH = HDLCH.test$ p_val, p_valLDLCH = LDLCH.test$p_val, 
                           p_val.SMK = SMK.test$p_val, p_val.PHYACT = PHYACT.test$p_val, p_val.ADTRN = ADTRN.test$p_val)

rownames(p_val_result) <- paste0("Module", 1:nrow(serum.modules))
p_val_result <- p_val_result[apply(p_val_result, MARGIN = 1, function(x) any(x <= 0.05/nrow(serum.modules))), ]



plot(phen$CHOL)
points(phen$CHOL[phen$ID %in% isaRowNames(serum, data.isa = serum.isa, type = 'isa', n = 93)], col = "red", pch = 19)

t.test(phen$AGE[phen$ID %in% id.ser], phen$AGE[!phen$ID %in% id.ser])
t.test(phen$AGE[phen$ID %in% id.ur], phen$AGE[!phen$ID %in% id.ur])$p.value

plot(phen$AGE, phen$GLUC, pch = 19)
points(phen$AGE[phen$ID %in% id.ser], phen$GLC[phen$ID %in% id.ser], col = "red", pch = 19)




# Phenotype - Module correaltion ------------------------------------------

# phenotype = dependent variable ("GLUC", "TRIG", "HDLCH", "CHOL", "LDLCH", "SBP")
# modules = independent variables. use the score to analyse the correlation with the phenotype!
# co-variables = "AGE", "SEX", "BMI", "SMK", "PHYACT", "ADTRN" 


# Multiple linear regression model with phenotype ~ module's score + covariables
covariables <- c("GLUC", "TRIG", "HDLCH", "CHOL", "LDLCH", "SBP", "ADTRN")

## selection of modules having higher robustness
urine.modules.uni[order(urine.modules.uni$rob, decreasing = T),]
urine.mod.sel <- urine.modules.uni[urine.modules.uni$rob > 40,]

## calculate score for every subjects given a certain number of features
## serum normalized
score_id <- last.iteration(as.matrix(serum), data.isa)

linear_regression_result <- matrix(NA, ncol(score_id), 7)
for(j in 1:length(covariables)){
  for (i in 1:ncol(score_id)){
    p <- summary(lm(phen[,covariables[j]] ~ score_id[,i] + SEX + AGE + BMI + PHYACT + SMK, data = phen))$coefficients[,4] [2]
    linear_regression_result[i,j] <- round(-log(p),2)
  }
}

colnames(linear_regression_result) <- covariables
rownames(linear_regression_result) <- paste0("Module", 1:ncol(score_id))

heatmap.2(linear_regression_result, cellnote = linear_regression_result, 
          notecex= 0.7, margins = c(7,7), main = "Linear regression result \n Normalized serum data", 
          cexRow = 0.8, cexCol = 1.2, xlab = "", ylab = "",
          notecol="black", density.info="none", trace="none", dendrogram = "none", 
          col = c("forestgreen", "red"), breaks = c(0, 3, max(linear_regression_result)), 
          Colv="NA", Rowv = "NA", keysize = 0.5, key =F)


## serum non-normalized
score_id <- last.iteration(as.matrix(serum_log), data.isa2)

linear_regression_result <- matrix(NA, ncol(score_id), 7)
for(j in 1:length(covariables)){
  for (i in 1:ncol(score_id)){
    p <- summary(lm(phen[,covariables[j]] ~ score_id[,i] + SEX + AGE + BMI + PHYACT + SMK, data = phen))$coefficients[,4] [2]
    linear_regression_result[i,j] <- round(-log(p),2)
  }
}

colnames(linear_regression_result) <- covariables
rownames(linear_regression_result) <- paste0("Module", 1:ncol(score_id))

heatmap.2(linear_regression_result, cellnote = linear_regression_result, 
          notecex= 0.7, margins = c(7,7), main = "Linear regression result \n non-normalized serum data", 
          cexRow = 0.8, cexCol = 1.2, xlab = "", ylab = "",
          notecol="black", density.info="none", trace="none", dendrogram = "none", 
          col = c("forestgreen", "red"), breaks = c(0, 3, max(linear_regression_result)), 
          Colv="NA", Rowv = "NA", keysize = 0.5, key =F)



# HDL and tot Chol significant
id <- isaRowNames(urine, data.isa = urine.uni, type = "isa", n = i)

plot(phen$AGE, phen$TRIG)
points(phen$AGE[phen$ID %in% id], phen$TRIG[phen$ID %in% id], col = "red", pch = 19)


## serum
serum.mod.sel <- serum.modules.uni[serum.modules.uni$rob > 30,]

linear_regression_result <- matrix(NA, ncol(serum.uni$rows), 6)
for(j in 1:length(covariables)){
  for (i in 1:ncol(serum.uni$rows)){
    p <- summary(lm(phen[,covariables[j]] ~ serum.uni$rows[,i] + SEX + AGE + BMI + PHYACT + SMK + ADTRN, data = phen))$coefficients[,4] [2]
    linear_regression_result[i,j] <- round(p,2)
  }
}


colnames(linear_regression_result) <- covariables
rownames(linear_regression_result) <- paste0("Module", 1:ncol(serum.uni$rows))

linear_regression_result <- linear_regression_result[order(serum.modules.uni$rowGroups, decreasing = T),]


heatmap.2(linear_regression_result, cellnote = linear_regression_result, notecex= 0.7, margins = c(7,7),
          main = "Linear regression result", cexRow = 0.8, cexCol = 1.2, xlab = "", ylab = "",
          notecol="black", density.info="none", trace="none", dendrogram = "none", 
          col = colorRampPalette(c("red", "green"))(n = 199), breaks = c(seq(0,0.05,length=100),             # for red
                                                                         seq(0.3,1,length=100)), Colv="NA", Rowv = "NA", 
          keysize = 0.7, key.title = "p-value", key.xlab = NA)

serum.modules.uni[order(serum.modules.uni$rowGroups, decreasing = T),]
