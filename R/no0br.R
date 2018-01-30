#' @title No zero branch length
#'
#' @name no0br
#' @rdname no0br
#'
#' @description This function does a small tweak on the tree to remove zero-length branch
#' between root and MRCA of the ingroup, and replace any negative branch lengths
#' with a small length
#'
#' @param phy an rooted tree
#'
#' @return returns the tree without zero-length or negative-length branch
#'
#' @export
no0br = function(phy) {

  if (!inherits(phy, "phylo"))
    stop(paste0(date(),"tree input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths"))

  if (!is.rooted(phy))
    stop(paste0(date(),": tree is not rooted"))

  if (phy$edge.length[1] != 0)
    stop(paste0(date(),": root has no zero-length branch, nothing to do"))

  tmp <- sort(phy$edge.length)
  tmp1 <- tmp[tmp>0]

  if (any(tmp < 0)) {

    warning(paste0(date(),": there are negative-length branches in the tree, replacing them with a small length"))
    phy$edge.length[phy$edge.length < 0] <- tmp1[1] / 1e3 ### add a small length for negative-length branches

  }

  phy$edge.length[1] <- tmp1[1]
  phy$edge.length[length(phy$edge.length)] <-
  phy$edge.length[length(phy$edge.length)] - phy$edge.length[1]

  phy

}
