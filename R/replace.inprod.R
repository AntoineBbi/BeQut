replace.inprod <-
  function (body.model, name_model, Data, param) {
    mt <- deparse(body.model, width.cutoff = 200L)
    if(name_model == "jags_lqm"){
      # regression from X
      ncX <- Data$ncX
      rplc <- paste(paste("beta[", 1:ncX, "] * X[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[1:ncX], X[i, 1:ncX])", rplc, mt, fixed = TRUE)
    }
    if(name_model %in% c("jags_lqmm","jags_qrjm")){
      # regression from X
      ncX <- Data$ncX
      rplc <- paste(paste("beta[", 1:ncX, "] * X[j, ", 1:ncX, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(beta[1:ncX], X[j, 1:ncX])", rplc, mt, fixed = TRUE)
      # regression on random effects U
      ncU <- Data$ncU
      rplc <- paste(paste("b[i, ", 1:ncU, "] * U[j, ", 1:ncU, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(b[i, 1:ncU], U[j, 1:ncU])", rplc, mt, fixed = TRUE)
    }
    # regression in survival submodel
    if(name_model == "jags_qrjm"){
      ncZ <- Data$ncZ
      rplc <- paste(paste("alpha[", 1:ncZ, "] * Z[i, ", 1:ncZ, "]", sep = ""), collapse = " + ")
      mt <- gsub("inprod(alpha[1:ncZ], Z[i, 1:ncZ])", rplc, mt, fixed = TRUE)
      if (param %in% c("value")) {
        rplc <- paste(paste("beta[", 1:ncX, "] * Xtime[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(beta[1:ncX], Xtime[i, 1:ncX])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("b[i, ", 1:ncU, "] * Utime[i, ", 1:ncU, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(b[i, 1:ncU], Utime[i, 1:ncU])", rplc, mt, fixed = TRUE)
        #
        rplc <- paste(paste("beta[", 1:ncX, "] * Xs[K * (i - 1) + k, ", 1:ncX, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(beta[1:ncX], Xs[K * (i - 1) + k, 1:ncX])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("b[i,", 1:ncU, "] * Us[K * (i - 1) + k, ", 1:ncU, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(b[i, 1:ncU], Us[K * (i - 1) + k, 1:ncU])", rplc, mt, fixed = TRUE)
      }
      if (param %in% c("sharedRE")) { # why add '&& ncU>1' in condition of this if ?
        ncU <- Data$ncU
        rplc <- paste(paste("alpha.assoc[", 1:ncU, "] * b[i,  ", 1:ncU, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(alpha.assoc[1:ncU], b[i, 1:ncU])", rplc, mt, fixed = TRUE)
      }
    }
    # if(name_model == "jags_mlqmm"){
    #   # regression from X
    #   ncX <- Data$ncX
    #   rplc <- paste(paste("beta[", 1:ncX, ", qq] * X[j, ", 1:ncX, "]", sep = ""), collapse = " + ")
    #   mt <- gsub("inprod(beta[1:ncX, qq], X[j, 1:ncX])", rplc, mt, fixed = TRUE)
    # }
    # end
    c("model", mt)
  }
