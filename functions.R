## matrix multiplication in the presence of NA: (Gabor's function)

na.multiply <- function(A, B) {
  M <- !is.na(A)
  modA <- A
  modA[!M] <- 0
  v <- modA %*% B
  w <- sqrt(M %*% B^2)
  w2 <- sqrt(apply(B^2, 2, sum))
  ret <- v/w * rep(w2, each=nrow(v))
  ret[ w==0 ] <- 0
  ret
}

## metabomatching input for phenomenal
## calls last.iteration function (to calculate the z-score for each of the metabolome features)

metabomatching_input <- function(data, data.isa, ppm, Col = TRUE){
  score <- last.iteration(as.matrix(data), data.isa, Col = Col)
  input <- data.frame(shift = as.numeric(ppm))
  for (i in 1:ncol(score)){
    input <- cbind(input, 
                   score[,i], 
                   rep(1, ncol(data)),
                   2*pnorm(-abs(score[,i])))
    colnames(input)[colnames(input) %in% c("score[, i]", "rep(1, ncol(data))",
                                           "2 * pnorm(-abs(score[, i]))")] <- c(sprintf("beta/m%03d", i), 
                                                                                sprintf("se/m%03d", i),
                                                                                sprintf("p/m%03d", i))
  }
  return(input)  
}


# Identity function -------------------------------------------------------
# this function check the identity between modules created 
# using different datasets (or using different parameters)
# 

## Identity function
identity <- function(data1, data2, data.isa1, data.isa2, sel = 0, Col = FALSE){
  id.matr <- matrix(NA, nrow = ncol(data.isa1$rows), ncol = ncol(data.isa2$rows))
  if(Col == TRUE){
    for(i in 1:ncol(data.isa1$rows)){
      id1 <- isaColNames(data = data1, type = "isa", data.isa = data.isa1, n = i)
      for(j in 1:ncol(data.isa2$rows)){
        id2 <- isaColNames(data = data2, type = "isa", data.isa = data.isa2, n = j)
        perc <- round(ifelse(length(id1) > length(id2), sum(id1 %in% id2)/length(id1), sum(id2 %in% id1)/length(id2)), 2)
        id.matr[i,j] <- perc
      }
    }
  } else {
    for(i in 1:ncol(data.isa1$rows)){
      id1 <- isaRowNames(data = data1, data.isa = data.isa1, n = i, type = "isa")
      for(j in 1:ncol(data.isa2$rows)){
        id2 <- isaRowNames(data = data2, type = "isa", data.isa = data.isa2, n = j)
        perc <- round(ifelse(length(id1) > length(id2), sum(id1 %in% id2)/length(id1), sum(id2 %in% id1)/length(id2)), 2)
        id.matr[i,j] <- perc
      }
    }
  }
  rownames(id.matr) <- paste0("mod", seq(1,ncol(data.isa1$rows)))
  colnames(id.matr) <- paste0("mod", seq(1,ncol(data.isa2$rows)))
  
  identity.sel <- id.matr[apply(id.matr, MARGIN = 1, function(x) any(x >= sel)), ]
  identity.sel <- identity.sel[,apply(identity.sel, MARGIN = 2, function(x) any(x >=sel))]
  
  return(identity.sel)
  
}


# other functions ---------------------------------------------------------

## classes extraction function: (for example to extract the age classes)
classes <- function(data, from, to, by, variable_name, variable_name2) {
  data_df <- data.frame(data)
  data_df$variable_class <- cut(as.numeric(unlist(data[colnames(data) %in% variable_name])), seq(from, to, by))
  return(as.data.frame(table(data_df[,colnames(data_df) %in% c("variable_class", variable_name2)])))
}


## t.test variables (not used anymore. We decided to use the linear regression with the z-score for each individuals)

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
