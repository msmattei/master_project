## run test using only female
## 

# load packages -----------------------------------------------------------
library(isa2)
library(gplots)

# Specify OS system: ------------------------------------------------------
### for windows
# pathOS <- "C://"
### for linux
pathOS <- "/media/mirjam/OS/"

# ISA (and other) functions  ----------------------------------------------
source(paste0(pathOS, 'Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project/ISA_functions.R'))
source(paste0(pathOS, 'Mimi/UNI/Master/MLS_BIOINFORMATICS/Master_Project_MLS/master_project/functions.R'))

# load data ---------------------------------------------------------------
## use of prenormalized serum maybe not necessary?? 
load(file = "../data/Rdata/metabo_data2017-10-12.Rdata") ## normalized-transformed data 
load(file = "../data/Rdata/pheno_data2017-10-20.Rdata") ## transformed phenotype data
sex     <- read.csv("../data/phenotype/all/sex.01.20170306.csv", h = F, sep = ",", stringsAsFactors = T)
sex$V1  <- paste0("id_", sex$V1)
sex     <- as.data.frame(cbind(sex$V2[na.omit(match(sex$V1, rownames(serum)))], serum))
serum.w <- sex[sex$V1 == 0, 2:ncol(sex)]
serum.m <- sex[sex$V1 == 1, 2:ncol(sex)]
urine.w <- sex[sex$V1 == 0, 2:ncol(sex)]
urine.m <- sex[sex$V1 == 1, 2:ncol(sex)]

files <- c("serum.w", "serum.m", "urine.w", "urine.m")
path1  <- "../result/metabomatching/gender/"
path2 <- "../result/linear_regression/gender/"

files <- c("serum", "serum_log", "urine", "urine_log")
path1  <- "../result/metabomatching/"
path2 <- "../result/linear_regression/"
# run isa -----------------------------------------------------------------
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


