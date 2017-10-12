
# Data information and load -----------------------------------------------

setwd("C:/Mimi/Stage_CBG/Data/NMR_data/hmdb_nmr_peak_lists")

txtfiles = list.files(path = "C:/Mimi/Stage_CBG/Data/NMR_data/hmdb_nmr_peak_lists", pattern = "HMDB*")
f1 = "List_of_peaks"
f2 = "List_of_Multiplets"
f3 = "List_of_Assignements"

PeakTableLable        = "of Peaks"
MultipletsTableLable  = "of Mu"
AssignmentsTableLable = "of Assig"


# Table of Peaks ----------------------------------------------------------

Table_of_Peaks  = NULL
HMDB.ID         = NULL
ppm             = NULL
adjusted_height = NULL

for (i in 1:length(txtfiles))	{
  HMDBtxtfile       <- readLines(txtfiles[i], warn = F)
    n_p <- grep(PeakTableLable, HMDBtxtfile)
    n_m <- grep(MultipletsTableLable, HMDBtxtfile)
    n_a <- grep(AssignmentsTableLable, HMDBtxtfile)
      tryCatch({
        if(HMDBtxtfile[1]=="No.\t(ppm)\tHeight"){
          f1 <- read.table(text = HMDBtxtfile[1 : (min(n_m, n_a)-1)], h=T, check.names = F)      
          } else {
            if(n_p+1 > max(n_m, n_a)){
              f1 <- read.table(text = HMDBtxtfile[(n_p+1) : (max(c(n_m, n_a, length(HMDBtxtfile))[c(n_m, n_a, length(HMDBtxtfile)) > n_p])-1)], h=T, check.names = F)
            } else {
              if(min(n_m, n_a) < (n_p+1)){
                f1 <- read.table(text = HMDBtxtfile[(n_p+1) : (max(n_m, n_a)-1)], h=T, check.names = F)
              } else {
                f1 <- read.table(text = HMDBtxtfile[(n_p+1) : (min(n_m, n_a)-1)], h=T, check.names = F)
              }
            }
        }
        if (colnames(f1[1])!="No.") {
            rbind(colnames(f1), f1)
            colnames(f1) <- c("No.", "(ppm)", "(Hz)", "Height")
          } else {
            f1=f1
        }
      HMDB.ID_i         <- rep(substr(txtfiles[i],1, 9), nrow(f1))
        ppm_i             <- f1$`(ppm)`
        if	(length(f1$Height)==0)	{
          adjusted_height_i=rep(NA, length(HMDB.ID_i))
          }	else	{
            adjusted_height_i=(f1$Height)/max(f1$Height)
        }
        HMDB.ID           <- c(HMDB.ID, HMDB.ID_i)
        ppm               <- c(ppm, ppm_i)
        adjusted_height   <- c(adjusted_height, adjusted_height_i)
    }, error=function(e){})
Table_of_Peaks <- data.frame(HMDB.ID, ppm, adjusted_height = round(adjusted_height, digits = 4))
}

Table_of_Peaks <- Table_of_Peaks[complete.cases(Table_of_Peaks),]
Table_of_Peaks <- unique(Table_of_Peaks)
Table_of_Peaks <- Table_of_Peaks[!Table_of_Peaks$adjusted_height<=0,]


# Table of Multiplets --------------------------------------------------------------

Table_of_Multiplets = NULL
HMDB.ID             = NULL
Multiplet           = NULL
ppm_start           = NULL
ppm_end             = NULL
Hs                  = NULL

for (i in 1:length(txtfiles))	{
  
  HMDBtxtfile       <- scan(txtfiles[i], what = character(), sep = "\t", strip.white=TRUE)
  tryCatch({
    #print(txtfiles[i])
      n_p <- grep(PeakTableLable, HMDBtxtfile)
      n_a <- grep(AssignmentsTableLable, HMDBtxtfile)
      n_m <- grep(MultipletsTableLable, HMDBtxtfile)
      n   <- length(HMDBtxtfile)
          f <- HMDBtxtfile[(n_m+1):(ifelse(n_m<n_a, ifelse(min(n_p,n_a)>n_m, min(n_p, n_a)-1, n_a-1), n))]
          f <- f[!f==""]
          if (f[10] %in% c("(ppm)", "Multiplet1 (ppm)") & length(f)/10 == floor(length(f)/10)) {
            f2 <- matrix(f, ncol = 10, byrow = T)
          } else {
          if (f[6] %in% c("(ppm)", "Multiplet1 (ppm)") & length(f)/6 == floor(length(f)/6)) {
            f2 <- matrix(f, ncol = 6, byrow = T)
            } else {
              if (f[7] %in% c("(ppm)", "Multiplet1 (ppm)") & length(f)/7 == floor(length(f)/7)) {
                f2 <- matrix(f, ncol = 7, byrow = T)
               } else {
                if(f[8] %in% c("(ppm)", "Multiplet1 (ppm)") & length(f)/8 == floor(length(f)/8)){
                  f2 <- matrix(f, ncol = 8, byrow = T)
                } else
                  f2 <- NULL
               }
            }
          }
          colnames(f2) <- f2[1,]
          ifelse(test = nrow(f2)==2, yes = f2 <- t(data.frame(f2[-1,], check.names = F)), no = f2 <- data.frame(f2[-1,], check.names = F))
          ifelse(colnames(f2)[3]=="(ppm)", f2 <- NULL, f2 <- f2)
          f2 <- as.data.frame(f2)
    HMDB.ID_i           <- rep(substr(txtfiles[i],1, 9), nrow(f2))
    Multiplet_i         <- as.character(f2$Multiplet1)
    Hs_i                <- as.numeric(levels(f2$Hs))[f2$Hs]
    ppm                 <- do.call("rbind", strsplit(as.character(f2[,colnames(f2)[colnames(f2) %in% c("(ppm)", "Multiplet1 (ppm)")]]), "..", fixed=TRUE))
    ppm_start           <- c(ppm_start, as.numeric(ppm[,1]))
    ppm_end             <- c(ppm_end, as.numeric(ppm[,2]))
    HMDB.ID             <- c(HMDB.ID, HMDB.ID_i)
    Multiplet           <- c(Multiplet, Multiplet_i)
    Hs                  <- c(Hs, Hs_i)
  }, error=function(e){}) # cat("ERROR :",conditionMessage(e), "\n")
  Table_of_Multiplets <- data.frame(HMDB.ID, Multiplet, Hs, ppm_start, ppm_end) 
}

Table_of_Multiplets <- Table_of_Multiplets[complete.cases(Table_of_Multiplets),]
Table_of_Multiplets <- unique(Table_of_Multiplets)


# Table of Assignements ---------------------------------------------------

Table_of_Assignments = NULL
HMDB.ID              = NULL
Multiplet            = NULL
Exp_ppm              = NULL


for (i in 1:length(txtfiles))	{
  
  HMDBtxtfile       <- scan(txtfiles[i], what = character(), sep = "\t", strip.white=TRUE)
  tryCatch({
    print(txtfiles[i])
        f <- HMDBtxtfile[(grep("of Assi", HMDBtxtfile)+1):length(HMDBtxtfile)]
        f <- f[!f==""]
        f3 <- matrix(f, ncol = 4, byrow = T)
        colnames(f3) <- f3[1,]
        ifelse(test = nrow(f3)==2, yes = f3 <- t(data.frame(f3[-1,], check.names = F)), no = f3 <- data.frame(f3[-1,], check.names = F))
    HMDB.ID_i            <- rep(as.numeric(substr(txtfiles[i],5, 9)), nrow(f3))
    Multiplet_i          <- f3$Multiplet
    Exp_ppm_i            <- (f3$`Exp. Shift (ppm)`)
    HMDB.ID              <- c(HMDB.ID, HMDB.ID_i)
    Multiplet            <- c(Multiplet, Multiplet_i)
    Exp_ppm              <- c(Exp_ppm, Exp_ppm_i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  Table_of_Assignments <- data.frame(HMDB.ID, Multiplet, Exp_ppm) 
}


# Add CAS number ----------------------------------------------------------

cas_ID               <- read.csv("C:/Mimi/Stage_CBG/Files/Metabolite/2098.csv")
Table_of_Peaks       <- merge(cas_ID, Table_of_Peaks, all.y = T)
Table_of_Peaks       <- Table_of_Peaks[order(Table_of_Peaks$HMDB.ID, Table_of_Peaks$ppm),]
Table_of_Multiplets  <- merge(cas_ID, Table_of_Multiplets, all.y = T)
Table_of_Multiplets  <- Table_of_Multiplets[order(Table_of_Multiplets$HMDB.ID, Table_of_Multiplets$Multiplet),]
Table_of_Assignments <- merge(cas_ID, Table_of_Assignments)

# Write table -------------------------------------------------------------

write.table(Table_of_Peaks, "C:/Mimi/Stage_CBG/Files/Metabolite/hmdb.20160809.import.mirjam.slop", quote = T, row.names = F)
write.table(Table_of_Multiplets, "C:/Mimi/Stage_CBG/Files/Metabolite/hmdb.20160809.import.mirjam.slom", quote = T, row.names = F)
# write.table(Table_of_Assignments, "C:/Mimi/Stage_CBG/Files/Metabolite/hmdb.20160602.import.mirjam.sloa.txt", quote = F, row.names = F)


# Selecting 180 uriary metabolites ----------------------------------------

slop.r <- read.table("C:/Mimi/Stage_CBG/Package/New/colaus.test/umdb.20160526.slop")
Table_of_Peaks$V4 <- as.numeric(substr(Table_of_Peaks$HMDB.ID, 5, 9))
HMDB_ID180 <- unique(Table_of_Peaks[,c(1,5)][Table_of_Peaks[,5] %in% unique(slop.r$V4),])
# write.table(HMDB_ID180, "C:/Mimi/Stage_CBG/Files/Metabolite/hmdb180.txt", row.names = F)
# HMDB_ID180 <- read.table("C:/Mimi/Stage_CBG/Files/Metabolite/hmdb180.txt", h=T)

Table_of_Peaks180 <- Table_of_Peaks[Table_of_Peaks$HMDB.ID %in% HMDB_ID180$HMDB.ID,]
Table_of_Multiplets180 <- Table_of_Multiplets[Table_of_Multiplets$HMDB.ID %in% HMDB_ID180$HMDB.ID,]
write.table(Table_of_Peaks180, "C:/Mimi/Stage_CBG/Files/Metabolite/hmdb.20160809.180.slop", quote = T, row.names = F)
write.table(Table_of_Multiplets180, "C:/Mimi/Stage_CBG/Files/Metabolite/hmdb.20160809.180.slom", quote = T, row.names = F)


# end ---------------------------------------------------------------------


