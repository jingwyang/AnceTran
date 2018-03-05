#' @title Function for estimating pair-wise TF binding distances based on sOU model
#'
#' @name distances
#'
#' @rdname distances
#' @description Stationary OU method is applied to estimate pair-wise TF binding distances
#' @param bsMat a TF-binding score matrix: column corresponds to binding score value;
#' row corresponds to othologous genes
#' @return A distance matrix shows the transcriptome distance for each paried species.
#' @export
TFdist.sou = function (bsMat = NULL) {

  object_n <- ncol(bsMat)

  dis.mat <- matrix(0, nrow = object_n, ncol = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- -log(cor(bsMat[,i],bsMat[,j]))

    }

  }

  dis.mat

}

