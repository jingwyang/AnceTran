#'
#' @title Function for estimating pair-wise TF binding distances based on sOU model
#'
#' @name distances
#'
#' @rdname distances
#' @description sOU method to estimate pair-wise TF binding distances
#'
#' @param bsMat a TF-binding score matrix: column corresponds to binding score value;
#' row corresponds to othologous genes
#' @return A distance matrix shows the sOU distance for each paried specis.
#' @export
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

