#' @title TF-binding score table extracted from a \code{taxaTF} class
#'
#' @name TFtab
#'
#' @description Generate an TF binding score table from a \code{taxaTF} class
#'
#' @param objects a vector of objects of class \code{taxonTF} or an object of class \code{taxaTF}
#' @param taxa one single character or a vector of characters specifying main taxa to be included in
#' the TF binding score table. taxa usually corresponds to species tissues or cell line types.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaTF} will be matched and included ("all" by default).
#' @param tf one single character or a vector of characters specifying transcription factor(s) to be included in
#' the TF binding score table.
#' If one single character "all" is given,
#' all the transcription factor included in the \code{taxaTF} will be matched and included ("all" by default).
#' @param rowindex a vector of numbers corresponded to indices of selecting rows
#' @param filtering a logical specifying whether to exclude genes with binding score equal to 0 in all taxaon and tfs.
#' @param normalize a logical specifying whether to perform quantile normalization on sample-specific biases.
#' @param logrithm a logical specifying whether to apply TF binding score log2 tranformation (TRUE by default).
#' @return An TF binding score table: column corresponds to median binding score value of all biological samples
#' within one taxa_TF group; row corresponds to othologous genes
#'
#' @export

TFtab = function (objects = NULL, taxa = "all", tf = "all",
                     rowindex = NULL, filtering = TRUE, normalize = TRUE, logrithm = TRUE){

  if (is.null(objects) || class(objects) != "taxaTF") {
    stop(paste0(date(), ": no valid taxaTF object input!"))
  }

  flag1 <- TRUE
  flag2 <- TRUE

  if (any(grepl("all",taxa, ignore.case = T))) {flag1 = FALSE}
  else { taxa <- gsub("\\s+","",taxa)}

  if (any(grepl("all",tf, ignore.case = T))) {flag2 = FALSE}
  else { tf <- gsub("\\s+","",tf)}

  #browser()

  tfbs_table <- NULL
  sample_names <- NULL
  # subsetting
  objects_n <- length(objects)

  if ( flag1 || flag2)

  {

    for (i in 1:objects_n)

    {
      if (flag1 && flag2) {
        if (any(grepl(objects[[i]]$taxon.name,taxa, ignore.case=T))
            && any(grepl(objects[[i]]$tf.name, tf, ignore.case=T)))
        {
          tfbs_table <- cbind(tfbs_table, apply(objects[[i]]$BindingScore.raw, 1, median))
          sample_names <- c(sample_names,
                            paste0(objects[[i]]$taxon.name,"_",objects[[i]]$tf.name))
        }

      } else {
        if (any(grepl(objects[[i]]$taxon.name,taxa,ignore.case=T))
            ||  any(grepl(objects[[i]]$subTaxon.name, tf, ignore.case=T)))
        {
          tfbs_table <- cbind(tfbs_table, apply(objects[[i]]$BindingScore.raw, 1, median))
          sample_names <- c(sample_names,
                            paste0(objects[[i]]$taxon.name,"_",objects[[i]]$tf.name))
        }
      }
    }
  } else {

    for (i in 1:objects_n) {
      tfbs_table <- cbind(tfbs_table, apply(objects[[i]]$BindingScore.raw, 1, median))
      sample_names <- c(sample_names,
                        paste0(objects[[i]]$taxon.name,"_",objects[[i]]$tf.name))
    }

  }

  if (is.null(tfbs_table)) {

    stop(paste0(date(),": taxa and tf names not found."))

  } else {

    row.names(tfbs_table) = objects[[1]]$gene.names
    colnames(tfbs_table) = sample_names

  }

  if (!is.null(rowindex)) {

    tfbs_table <- tfbs_table[rowindex,]

  }

  if(filtering){
    keep <- rowSums(tfbs_table == 0) < ncol(tfbs_table)
    tfbs_table <- tfbs_table[keep, ]
  }

  if(normalize){
    tfbs_table <- normalizeQuantiles(tfbs_table)
  }

  if (logrithm) {

    tfbs_table <- apply(tfbs_table, c(1,2), function(x) log2(x+1))

  }

  tfbs_table
}
