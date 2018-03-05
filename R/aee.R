#' @title Ancestral Transcriptome State Inference
#'
#' @name aee
#' @rdname aee
#'
#'
#' @description This function is to inference ancestral transcriptome state
#' and related statistical uncertainty based on RNA-seq expression data
#' or ChIP-seq TF-binding data.
#'
#' @param x a vector of known expression or TF binding score values of genes, preferably log2-transformed
#' @param phy a phylogenetic tree in the form of object "phylo"
#' @param mat a matrix generated from "varMatInv" function
#' @param select indicate if descendents of the node or all tips should be used
#' @param CI a logical specifying whether to return the 95% confidence intervals
#' of the estimated ancestral transcriptome value.
#'
#' @return returns a list containing estimated ancestral transcriptome states
#' as well as other requested parameters
#' @export
aee = function(x, phy, mat, select = c("descendents","all") , CI = TRUE) {

  ### checking input formats
  if (!inherits(phy,"phylo"))
    stop(paste0(date(),": \"phy\" input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths which is a necessity for \"aee\""))

  if (!is.null(names(x))) {
    if (all(names(x) %in% phy$tip.label))
      x <- x[phy$tip.label]
    else warining(paste0(date(),
                         "characters do not match perfectly between expression profile vector names and tree tip labels,
                         only using tree tip labels in the following analysis"))
  }

  #browser()

  ### checking if the tree is rooted
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)

  if (n_tip != n_node + 1)
    stop(paste0(date(),"tree is not rooted, please make sure tree is properly rooted. "))

  ancestral <- list()

  tran <- numeric(length = n_tip + n_node) ### expression values vector initiate

  tran[1:n_tip] <- if (is.null(names(x))) x else as.numeric(x[phy$tip.label]) ### given expression values of tips
  tran[(n_tip+1):(n_tip+n_node)] <- NA ### ancestral nodes expression values to be estimated

  tr_edges <- phy$edge

  ###

  #select <- match.arg(select)

  ### using children nodes to estimate ancestral nodes' expression

  if (0) { ### direct child nodes problematic? start

    while (any (is.na(tran[(n_tip+2):(n_tip+n_node)]))) {
      ### looping while ancestral expression values are not computed at all internal nodes except for root

      for (i in (n_tip+2):(n_tip+n_node)) {

        if (is.na(tran[i])) { ### node that has not been computed yet

          child_nodes <- tr_edges[tr_edges[,1] == i, 2]

          if (any(is.na(tran[child_nodes]))) ### bypass empty children nodes
            next

          mu <- mean(tran[child_nodes])

          beta <- unlist(lapply(child_nodes, function(x) - mat[x,i] / mat[i,i]))
          beta0 <- mu * (1 - sum(beta))

          tran[i] <- beta0 + sum(beta * tran[child_nodes])
        }

      }

    }

    ### for the root node

    child_nodes <- tr_edges[tr_edges[,1] == n_tip+1,2]

    mu <- mean(tran[child_nodes])

    beta <- unlist(lapply(child_nodes, function(x) - mat[x,n_tip+1] / mat[n_tip+1,n_tip+1]))
    beta0 <- mu * (1 - sum(beta))

    tran[n_tip+1] <- beta0 + sum(beta * tran[child_nodes])

  } ### problematic? end


  if (select == "descendents") {
    ### estimating ancestral expression with all decendant tips
    for (i in (n_tip+1):(n_tip+n_node)) {

      child_nodes <- getDescendants(phy, i)

      child_tips <- child_nodes[child_nodes < (n_tip+1)] ### only tips

      mu <- mean(tran[child_tips])

      beta <- unlist(lapply(child_tips, function(x) - mat[x,i] / mat[i,i]))
      beta0 <- mu * (1 - sum(beta))

      tran[i] <- beta0 + sum(beta * tran[child_tips])

    }

  } else {

    all_tips <- 1:n_tip

    for (i in (n_tip+1):(n_tip+n_node)) {

      mu <- mean(tran[all_tips])

      beta <- unlist(lapply(all_tips, function(x) - mat[x,i] / mat[i,i]))
      beta0 <- mu * (1 - sum(beta))

      tran[i] <- beta0 + sum(beta * tran[all_tips])

    }

  }

  ancestral$est <- tran[(n_tip+1):(n_tip+n_node)]
  ### if calculating confidence interval

  if (CI) {

    ci95 <- matrix(nrow = n_node, ncol = 2)

    for (i in (n_tip+1):(n_tip+n_node)) { ### for every node

      tmp = sqrt(1 / mat[i,i]) * qnorm(0.025)
      ci95[(i-n_tip),] = c(tran[i] + tmp, tran[i] - tmp)

    }

    ancestral$ci95 <- ci95
  }

  ancestral
  }
