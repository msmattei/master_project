## useful isa functions
# Modules informations -----------------------------------------------------
isaModules <- function(data.isa, type = c("isa", "ppa")) {
  if(type == "isa"){
    colGroups <- matrix(NA, ncol(data.isa$rows))
    rowGroups <- matrix(NA, ncol(data.isa$rows))
    for (i in 1:ncol(data.isa$rows)) {
      colGroups[i] <- sum(data.isa$columns[,i] != 0)
      rowGroups[i] <- sum(data.isa$rows[,i] != 0)
    }
    modules <- data.frame(colGroups, rowGroups, freq = data.isa$seeddata$freq,  thr.row = data.isa$seeddata$thr.row, thr.col = data.isa$seeddata$thr.col)
    colnames(modules) <- c("colGroups", "rowGroups", "freq", "thr.row", "thr.col")
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
    modules <- data.frame(colGroups, row1Groups, row2Groups)
    colnames(modules) <- c("colGroups", "row1Groups","row2Groups")
    rownames(modules) <- paste0("Module", 1:nrow(modules))
  }
  print(modules)
}
# Extract module informations ----------------------------------------------------
isaRowNames <- function(data, data2 = NULL, data.isa, type = "isa", n){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
    return(colnames(module.n))
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- t(data)[isaRow1, isaCol, drop=FALSE]
    module.n2 <- t(data2)[isaRow2, isaCol, drop=FALSE]
    return(list(data1 = colnames(module.n1), data2 = colnames(module.n2)))
  }
}
isaColNames <- function(data, data2 = NULL, data.isa, type = "isa", n){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
    return(rownames(module.n))
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- t(data)[isaRow1, isaCol, drop=FALSE]
    module.n2 <- t(data2)[isaRow2, isaCol, drop=FALSE]
    return(list(data1 = rownames(module.n1), data2 = rownames(module.n2)))
  }
}

# Module visualization ----------------------------------------------------

isa2image <- function(data, data2 = NULL, data.isa, type = "isa", n, name1 = NULL, name2 = NULL, cex = 0.6, all = FALSE){
  if(type == "isa"){
    isaRow = data.isa$rows[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
    image(module.n, axes = F, main = paste("Module",  n))
    mtext(text=colnames(module.n), side=2, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = cex)
    mtext(text=rownames(module.n), side=1, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = cex)
    # if(all == TRUE){
    #   allCol <- c(colnames(module.n), rownames(data)[!rownames(data) %in% colnames(module.n)])
    #   allRow <- c(rownames(module.n), colnames(data)[!colnames(data) %in% rownames(module.n)])
    #   image(as.matrix(t(data[allCol, allRow])), axes = F)
    # }
  }
  if(type == "ppa"){
    isaRow1 = data.isa$rows1[, n] != 0
    isaRow2 = data.isa$rows2[, n] != 0
    isaCol = data.isa$columns[, n] != 0
    module.n1 <- t(data)[isaRow1, isaCol, drop=FALSE]
    module.n2 <- t(data2)[isaRow2, isaCol, drop=FALSE]
    par(mfrow = c(1,2))
    image(module.n1, axes = F, main = paste("Module",  n, name1))
    mtext(text=colnames(module.n1), side=2, line=0.3, at=seq(0,1,l=ncol(module.n1)), las=2, cex = cex)
    mtext(text=rownames(module.n1), side=1, line=0.3, at=seq(0,1,l=nrow(module.n1)), las=2, cex = cex)
    image(module.n2, axes = F, main = paste("Module",  n, name2))
    mtext(text=colnames(module.n2), side=2, line=0.3, at=seq(0,1,l=ncol(module.n2)), las=2, cex = cex)
    mtext(text=rownames(module.n2), side=1, line=0.3, at=seq(0,1,l=nrow(module.n2)), las=2, cex = cex)
    par(mfrow = c(1,1))
  }
  if(all == TRUE){
    isaRow <- c(rownames(data)[data.isa$rows[, n] != 0], rownames(data)[data.isa$rows[, n] == 0])
    isaCol <- c(colnames(data)[data.isa$columns[, n] != 0], colnames(data)[data.isa$columns[, n] == 0])
    module.n   <- t(as.matrix(data[isaRow, isaCol, drop=FALSE]))
    image(module.n, axes = F, main = paste("Module",  n))
    mtext(text=colnames(module.n), side=2, line=0.3, at=seq(0,1,l=ncol(module.n)), las=2, cex = cex)
    mtext(text=rownames(module.n), side=1, line=0.3, at=seq(0,1,l=nrow(module.n)), las=2, cex = cex)
  }
}




