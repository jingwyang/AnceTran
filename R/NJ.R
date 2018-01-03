#' @title Neighbor-joining
#'
#' @name NJ
#'
#' @rdname nj
#'
#' @description The famous neighbor-joining tree estimation function from Saitou and Nei (1987).
#' @param x a distance matrix shows the sOU distances for each paried specis
#' @export
NJ = function (x) {


  reorder(nj(x), "postorder")

}
