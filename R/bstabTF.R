bstabTF = function (objects = NULL, taxa = "all", tf = "all",
                     rowindex = NULL, filtering = TRUE, normalize = TRUE, logrithm = TRUE)

{

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
