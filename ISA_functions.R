# IMPORTANT: original code issued from https://github.com/gaborcsardi/ISA

## useful isa functions

# Modules informations -----------------------------------------------------
# extract information about each module. Contains
# number of columns, number of rows, robustness and threshold used (for each module)
# possible to use this function for isa and for ppa

isaModules <- function(data.isa, type = "isa") {
  if(type == "isa"){
    colGroups <- matrix(NA, ncol(data.isa$rows))
    rowGroups <- matrix(NA, ncol(data.isa$rows))
    for (i in 1:ncol(data.isa$rows)) {
      colGroups[i] <- sum(data.isa$columns[,i] != 0)
      rowGroups[i] <- sum(data.isa$rows[,i] != 0)
    }
    modules <- data.frame(colGroups, rowGroups, rob = round(data.isa$seeddata$rob, 2),  
                          thr.col = data.isa$seeddata$thr.col, thr.row = data.isa$seeddata$thr.row)
    colnames(modules) <- c("colGroups", "rowGroups", "rob", "thr.col", "thr.row")
    rownames(modules) <- paste0("Module", 1:nrow(modules))
  }
  if(type == "ppa"){
    colGroups <- matrix(NA, ncol(data.isa$rows1))
    row1Groups <- matrix(NA, ncol(data.isa$rows1))
    row2Groups <- matrix(NA, ncol(data.isa$rows2))
    for (i in 1:ncol(data.isa$rows1)) {
      colGroups[i] <- sum(data.isa$columns[,i] != 0)
      row1Groups[i] <- sum(data.isa$rows1[,i] != 0)
      row2Groups[i] <- sum(data.isa$rows2[,i] != 0)
    }
    modules <- data.frame(colGroups, row1Groups, row2Groups, rob = data.isa$seeddata$rob,  
                          thr.col = data.isa$seeddata$thr.col,
                          thr.row1 = data.isa$seeddata$thr.row1, 
                          thr.row2 = data.isa$seeddata$thr.row2)
    colnames(modules) <- c("colGroups", "row1Groups","row2Groups", "rob", "thr.col", "thr.row1", "thr.row2")
    rownames(modules) <- paste0("Module", 1:nrow(modules))
  }
  return(modules)
}

# Extract informations ----------------------------------------------------
# extract row names for a module of interest
isaRowNames <- function(data, data.isa, n, data2 = NULL, type = "isa"){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- as.matrix(data[isaRow, isaCol, drop=FALSE])
    return(rownames(module.n))
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- data[isaRow1, isaCol, drop=FALSE]
    module.n2 <- data2[isaRow2, isaCol, drop=FALSE]
    return(list(data1 = rownames(module.n1), data2 = rownames(module.n2)))
  }
}

# extract column names from a module of interest
isaColNames <- function(data, data.isa, n, data2 = NULL, type = "isa"){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- as.matrix(data[isaRow, isaCol, drop=FALSE])
    return(colnames(module.n))
  }
  if(type == "ppa"){
    isaCol = data.isa$columns[, n] != 0
    return(colnames(module.n1), data2 = rownames(module.n2))
  }
}

# extracting the score from a module of interest (for row and columns)
isaScore <- function(data, data.isa, n){
    isaRow = data.isa$rows[, n] != 0
    isaNames <- rownames(data)[isaRow]
    scoreRow <- data.frame(RowScore = data.isa$row[isaRow, n], row.names = isaNames)

    isaCol = data.isa$columns[, n] != 0
    isaNames <- colnames(data)[isaCol]
    scoreCol <- data.frame(ColScore = data.isa$column[isaCol, n], row.names = isaNames)

  return(list(scoreRow, scoreCol))
}

# ISA run -----------------------------------------------------------------
# run the different isa steps with the parameteres that better fits to my needs
# i.e. aiming to have:
# 1) a way to reproduce the results (need to create row.seeds in a reproducible way)
# 2) Choose the correlation value to consider two module as equal 
#    (cor.limit, from the isa.unique function)
# 

isa.run <- function(data, 
                    thr.row=seq(0.5,1.5,by=0.5),
                    thr.col=seq(0.5,1.5,by=0.5),
                    row.seeds, 
                    direction=c("updown", "updown"),
                    cor.limit = 0.9) {
  ## Normalize the matrix (using isa.normalize() function)
  normed.data <- isa.normalize(data)
  
  ## Determine thresholds (expand.grid() create a data frame with all possible combination of thresholds)
  thr.list <- expand.grid(thr.row=thr.row, thr.col=thr.col)
  thr.list <- unlist(apply(thr.list, 1, list), recursive=FALSE) # make a list with the possible thresholds combination
  
  ## Do the ISA, for all thresholds (lapply())
  isaresults <- lapply(thr.list, function(x)
    isa.iterate(normed.data,
                row.seeds=row.seeds,
                thr.row=x["thr.row"],
                thr.col=x["thr.col"],
                direction=direction))
  
  ## Make it unique for every threshold combination
  isaresults <- lapply(isaresults, function(x)
    isa.unique(normed.data, x, cor.limit = cor.limit))
  
  ## Filter according to robustness
  isaresults <- lapply(isaresults, function(x)
    isa.filter.robust(data=data,
                      normed.data=normed.data,
                      isares=x,
                      row.seeds=row.seeds))
  
  ## Merge them
  result <- list()
  result$rows <- do.call(cbind, lapply(isaresults, "[[", "rows"))
  result$columns <- do.call(cbind, lapply(isaresults, "[[", "columns"))
  result$seeddata <- do.call(rbind, lapply(isaresults, "[[", "seeddata"))
  result$rundata <- isaresults[[1]]$rundata
  result$rundata$N <- sum(sapply(isaresults, function(x) x$rundata$N))
  
  ## Another filtering
  result <- isa.unique(normed.data, result, cor.limit = cor.limit)
  
  ## We are done
  result
}

# Module visualization ----------------------------------------------------
# heatmap of one module of interest

isa2image <- function(data, 
                      data.isa, 
                      n, 
                      type = "isa", 
                      data2 = NULL, 
                      name1 = NULL, 
                      name2 = NULL, 
                      cex = 0.6, 
                      color1 = "red", 
                      color2 = "white", 
                      all = FALSE){
  if(type == "isa"){
    colors <- colorRampPalette(c(color1, color2))(n = 10000)
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- as.matrix(data[isaRow, isaCol, drop=FALSE])
    ColorUsed <- colors[round(1+(min(module.n, na.rm = T)-min(data, na.rm = T))*10000/
                                (max(data, na.rm = T)-min(data, na.rm = T))) : round(
                                  (max(module.n, na.rm = T)-min(data, na.rm = T))*10000/(max(data, na.rm = T)-min(data, na.rm = T)) )]
    if(all == TRUE){
      allRow <- c(rownames(module.n), rownames(data)[!rownames(data) %in% rownames(module.n)])
      allCol <- c(colnames(module.n), colnames(data)[!colnames(data) %in% colnames(module.n)])
      image(data[allRow, allCol], axes = F, col=colors, main = paste("Module",  n))
      v <- nrow(module.n)/nrow(data)
      h <- ncol(module.n)/ncol(data)
      abline(h = h, v = v)
    } else {
      image(t(module.n), axes = F, main = paste("Module",  n), col=ColorUsed)
      mtext(text=rownames(module.n), side=2, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = cex)
      mtext(text=colnames(module.n), side=1, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = cex)
    }
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- data[isaRow1, isaCol, drop=FALSE]
    module.n2 <- data2[isaRow2, isaCol, drop=FALSE]
    par(mfrow = c(1,2))
    image(module.n1, axes = F, main = paste("Module",  n, name1))
    mtext(text=colnames(module.n1), side=2, line=0.3, at=seq(0,1,l=ncol(module.n1)), las=2, cex = cex)
    mtext(text=rownames(module.n1), side=1, line=0.3, at=seq(0,1,l=nrow(module.n1)), las=2, cex = cex)
    image(module.n2, axes = F, main = paste("Module",  n, name2))
    mtext(text=colnames(module.n2), side=2, line=0.3, at=seq(0,1,l=ncol(module.n2)), las=2, cex = cex)
    mtext(text=rownames(module.n2), side=1, line=0.3, at=seq(0,1,l=nrow(module.n2)), las=2, cex = cex)
    par(mfrow = c(1,1))
    if(all == TRUE){
      isaRow <- c(rownames(data)[data.isa$rows[, n] != 0], rownames(data)[data.isa$rows[, n] == 0])
      isaCol <- c(colnames(data)[data.isa$columns[, n] != 0], colnames(data)[data.isa$columns[, n] == 0])
      module.n   <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
      image(module.n, axes = F, main = paste("Module",  n))
      mtext(text=colnames(module.n), side=2, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = cex)
      mtext(text=rownames(module.n), side=1, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = cex)
    }
    
  }
}



# additionary iteration step ----------------------------------------------
## to calculate the score for each row (or column) in the modules
## multiply the gene vector with the whole matrix to have the 
## correspondent score for each idividuals! 
## Like doing an extra iteration without thresholding
## 
## IMPORTANT: the output will be a z-score (not the direct output of the isa function)
## z score can be used to calculate the corresponding p-value (using the pt() function)


last.iteration <- function(data, data.isa, Col = FALSE, type = "isa", data2 = FALSE){
  data.norm <- isa.normalize(data)
  if(type == "isa"){
    if(Col == TRUE) {
      score <- na.multiply(data.norm$Er, data.isa$rows)
      score <- apply(score, MARGIN = 2, function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T))
    } else {
      score <- na.multiply(data.norm$Ec, data.isa$columns)
      score <- apply(score, MARGIN = 2, function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T))
    }
  }
  if(type == "ppa"){
    if(data2 == TRUE){
      rows <- "rows2"
    } else {
      rows <- "rows1"
    }
    if(Col == TRUE) {
      score <- na.multiply(data.norm$Er, data.isa$rows)
      score <- apply(score, MARGIN = 2, function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T))
    } else {
      score <- na.multiply(data.norm$Ec, data.isa$columns)
      score <- apply(score, MARGIN = 2, function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T))
    }
  }
  return(score)
}


