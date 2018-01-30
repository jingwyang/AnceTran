#' @title Neighbor-joining tree
#'
#' @name NJ
#'
#' @rdname nj
#'
#' @description The famous neighbor-joining tree estimation function from Saitou and Nei (1987).
#' @param x a distance matrix shows the transcriptome distances for each paried species
#' @export
NJ = function (x) {


  reorder(nj(x), "postorder")

}
