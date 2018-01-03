TFdist.sou = function (bsMat = NULL) {

  object_n <- ncol(bsMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- -log(cor(bsMat[,i],bsMat[,j]))

    }

  }

  dis.mat

}

