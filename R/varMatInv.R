#' @title Generate an inversed variance-covariance matrix from transcriptome profiles across species
#'
#' @name varMatInv
#' @rdname varMatInv
#'
#' @description This function generate an inversed variance-covariance matrix from expression
#' or TF-binding profiles of one-to-one orthologous genes across species
#' @param dismat a distance matrix for paired transcriptome profiles
#' @param tran_table transcriptome value matrix extracted from objects of class
#' \code{taxaTF} or class \code{taxaExp}
#' @param phy an rooted transcriptome character tree
#' @return returns an inversed variance-covariance matrix
#' @export
varMatInv = function(dismat , tran_table, phy) {

  if (!inherits(phy, "phylo"))
    stop(paste0(date(),"tree input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths which is a necessity for \"varMatInv\""))

  if (!all(row.names(dismat) %in% phy$tip.label ))
    stop(paste0(date(),": taxa or tf names do not match perfectly with tree tip labels, please check them."))

  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)

  if (n_tip != n_node + 1)
    stop(paste0(date(),"tree is not rooted, please make sure tree is properly rooted. "))

  ### extract distances from the tree
  nodes_dist <- dist.nodes(phy)
  corrmat <- apply(nodes_dist, c(1,2), function(x) exp(-x))

  ### stationary variance
  #exp_table <- exptabTE(objects, taxa = "all", subtaxa = subtaxa, logrithm = T)
  stat_var <- mean(apply(tran_table, 2, var))
  var_corrmat <- stat_var * corrmat

  solve(var_corrmat)


}
