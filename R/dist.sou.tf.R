dist.sou = function (expMat = NULL) {

  object_n <- ncol(expMat)
  #gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      #V11 <- var(expMat[,i])
      #V22 <- var(expMat[,j])
      #V12 <- cov(expMat[,i], expMat[,j])

      dis.mat[j,i] <- -log(cor(expMat[,i],expMat[,j]))

    }

  }

  dis.mat

}

