varMatInv = function(dismat , exp_table, phy, taxa = "all", subtaxa) {

  if (!inherits(phy, "phylo"))
    stop(paste0(date(),"tree input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths which is a necessity for \"varMatInv\""))

  if (length(subtaxa) > 1 || subtaxa == "all")
    stop(paste0(date(),": only one subtaxon are allowed here"))

  #dismat <- expdist(objects, taxa = taxa, subtaxa = subtaxa, method = "sou") ### using -ln(rho) to estimate pairwise expression distance

  if (!all(row.names(dismat) %in% phy$tip.label ))
    stop(paste0(date(),": taxa or subtaxa names do not match perfectly with tree tip labels, please check them."))

  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)

  if (n_tip != n_node + 1)
    stop(paste0(date(),"tree is not rooted, please make sure tree is properly rooted. "))

  ### extract distances from the tree
  nodes_dist <- dist.nodes(phy)
  corrmat <- apply(nodes_dist, c(1,2), function(x) exp(-x))

  ### stationary variance
  #exp_table <- exptabTE(objects, taxa = "all", subtaxa = subtaxa, logrithm = T)
  stat_var <- mean(apply(exp_table, 2, var))
  var_corrmat <- stat_var * corrmat

  solve(var_corrmat)


}
