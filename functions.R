
# Identity function -------------------------------------------------------
# this function check the identity between modules created 
# using different datasets, or parameters
# 

## Identity function
identity <- function(data1, data2, data.isa1, data.isa2, modules1, modules2, sel = 0, Col = FALSE){
  id.matr <- matrix(NA, nrow = nrow(modules1), ncol = nrow(modules2))
  if(Col == TRUE){
    for(i in 1:nrow(modules1)){
      id1 <- isaColNames(data = data1, type = "isa", data.isa = data.isa1, n = i)
      for(j in 1:nrow(modules2)){
        id2 <- isaColNames(data = data2, type = "isa", data.isa = data.isa2, n = j)
        perc <- round(ifelse(length(id1) > length(id2), sum(id1 %in% id2)/length(id1), sum(id2 %in% id1)/length(id2)), 2)
        id.matr[i,j] <- perc
      }
    }
  } else {
    for(i in 1:nrow(modules1)){
      id1 <- isaRowNames(data = data1, type = "isa", data.isa = data.isa1, n = i)
      for(j in 1:nrow(modules2)){
        id2 <- isaRowNames(data = data2, type = "isa", data.isa = data.isa2, n = j)
        perc <- round(ifelse(length(id1) > length(id2), sum(id1 %in% id2)/length(id1), sum(id2 %in% id1)/length(id2)), 2)
        id.matr[i,j] <- perc
      }
    }
  }
  rownames(id.matr) <- paste0("mod", seq(1,nrow(modules1)))
  colnames(id.matr) <- paste0("mod", seq(1,nrow(modules2)))
  
  identity.sel <- id.matr[apply(id.matr, MARGIN = 1, function(x) any(x %in% seq(sel, 0.99, by = 0.001))), ]
  identity.sel <- identity.sel[,apply(identity.sel, MARGIN = 2, function(x) any(x %in% seq(sel, 0.99, by = 0.001)))]
  
  return(identity.sel)
  
}


# other functions ---------------------------------------------------------

## t.test variables 

t.test.variable <- function(phenotype.data, metabo.data, metabo.isa, variable1){
  p_val = NULL
  for(i in 1:nrow(metabo.isa$seeddata)){
    id <- isaRowNames(data = metabo.data, data.isa = metabo.isa, type = "isa", n = i)
    if(length(id)==1) {
      p <- NA
    } else {
      if(length(na.omit(phenotype.data[colnames(phenotype.data) %in% variable1][phenotype.data$ID %in% id,])) == 1) {
        p <- NA
      } else {
        p <- format(t.test(na.omit(phenotype.data[colnames(phenotype.data) %in% variable1][phenotype.data$ID %in% id,]), 
                           phenotype.data[colnames(phenotype.data) %in% variable1][!phenotype.data$ID %in% id,])$p.value, 
                    digits = 2)
      }
    }
    p_val <- c(p_val, as.numeric(p))
  }
  return(data.frame(modules = paste0("Module", 1:nrow(metabo.isa$seeddata)), p_val))  
}