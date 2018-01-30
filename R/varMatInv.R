#' @title Generate an inversed variance-covariance matrix from transcriptome profiles across species
#'
#' @name varMatInv
#' @rdname varMatInv
#'
#' @description This function generate an inversed variance-covariance matrix from expression
#' or TF-binding profiles of one-to-one orthologous genes across species
#'
#' @param objects a vector of objects of class \code{taxonTF} or an object of class \code{taxaTF}
#' @param phy an rooted expression or TF-binding character tree
#' @param taxa one single character or a vector of characters sepcifying taxa to generate
#' an inversed variance-covariance matrix.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaTF} will be matched and included ("all" by default).
#' @param subtaxa one single character specifying sub taxa to be included in generating
#' an inversed variance-covariance matrix
#'
#' @return returns an inversed variance-covariance matrix
#' @export
varMatInv = function(dismat , tran_table, phy, taxa = "all", tf) {

  if (!inherits(phy, "phylo"))
    stop(paste0(date(),"tree input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths which is a necessity for \"varMatInv\""))

  if (length(tf) > 1 || tf == "all")
    stop(paste0(date(),": only one subtaxon are allowed here"))

  #dismat <- expdist(objects, taxa = taxa, subtaxa = subtaxa, method = "sou") ### using -ln(rho) to estimate pairwise expression distance

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
