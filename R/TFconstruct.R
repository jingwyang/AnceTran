#' @title constrct a taxaTF object.
#'
#' @description  \code{TFconstrct} is a constructor function for \code{taxaTF} objects.
#' This fuction takes in a binding score file
#' and constructs a \code{taxaTF} object from which user can extract information
#' for display or for further analysis.
#' @name TEconstruct
#'
#' @rdname TEconstruct
#'
#' @param BSFile a text file contains raw TF binding score values.
#' Row names correspond with gene names,
#' and colum names correspond with species and trasncription factor names.
#' @param taxa one single string or a vector of strings specifying main taxa selected for
#' constructing \code{taxaTF} object.
#' If one single string "all" is given,
#' all the taxa in the row names will be matched and selected ("all" by default).
#' @param tf one single string or a vector of strings sepcifying transcription factor(s) selected for
#' constructing \code{taxaTF} object.
#' If one single string "all" is given,
#' all the transcription factor(s) in the row names will be matched and selected ("all" by default).
#' @param verbose a logical specifying whether to print more information on the screen
#' while constructing \code{taxaTF} objects (FALSE by default).
#'
#'@return returns an object of class \code{TaxaTF} (S3 class, a list of \code{taxonTF} objects).
#'
#' @examples
#'
#' taxa.objects = TFconstruct(BSFile = system.file('extdata/CEBPA_meanIntensity_all_Mouse.txt', package = 'AnceTran'),
#'    taxa = "all", tf = "all", verbose = FALSE)
#'
#' @export
TFconstruct = function(BSFile=NULL, taxa="all", tf="all", verbose=FALSE) {

  # check file handle
  if(is.null(BSFile)){
    stop(paste0(date(),": must provide binding score file path"))
  }

  # check file existance
  if(!file.exists(BSFile)){
    stop(paste0(date(),": fail to open file, check your filename or path"))
  }

  # input
  binding.score.df <- read.table(BSFile, header=T)
  row.names(binding.score.df) <- binding.score.df[, 1]
  binding.score.df <- binding.score.df[, -1]

  # gene number and taxon number
  gene_n <- nrow(binding.score.df)

  # get taxon names from read counts file
  taxon_names <- unique(lapply(colnames(binding.score.df), function(x) unlist(strsplit(x, "_"))[1]))
  taxon_n <- length(taxon_names)

  #normalize<-match.arg(normalize)
  #message(paste0(date(), ": using ", normalize, " to normalize raw TF binding scores"))

  cat("\n")
  message(paste0(date(),": start constructiong TF objects"))

  if (!any(grepl("all", taxa, ignore.case = T))) {

    taxon_names <- gsub("\\s+", "", taxa)
    taxon_n <- length(taxon_names)

  }

  # get transcription factor number
  tf_names <- unique(lapply(colnames(binding.score.df), function(x) unlist(strsplit(x, "_"))[2]))
  tf_n <- length(tf_names)

  if (!any(grepl("all", tf, ignore.case = T))) {

    tf_names <- gsub("\\s+", "", tf)
    tf_n <- length(tf_names)

  }
  message(paste0(date(),": total taxon number ", taxon_n))
  message(paste0(date(),": total transcription factor number ", tf_n))

  title <- lapply(colnames(binding.score.df), function(x) unlist(strsplit(x, "_"))[1]) # taxon names
  subtitle <- lapply(colnames(binding.score.df), function(x) unlist(strsplit(x,"_"))[2]) # TF names
  first_two_names <- unique(paste(title,subtitle,sep="_"))

  index <- intersect(unlist(lapply(taxon_names, function(x) grep(x, first_two_names, ignore.case = T))),
                     unlist(lapply(tf_names,function(x) grep(x, first_two_names, ignore.case = T))))

  objects_names <- first_two_names[index]

  objects_number <- length(objects_names)

  #browser()
  cat("\n")
  # get gene names
  #gene.names <- read.counts.df[,1]

  message(paste0(date(),": now constructing ",objects_number, " TE objects..."))

  if (!verbose) progbar <- txtProgressBar(style = 3)

  # initialization

  taxonTF.objects <- vector("list",length = objects_number)
  # the number of TE objects constructed is based on seleted numnber

  # for each taxon
  objects_counter <- 0
  for (i in 1:objects_number) {

    #browser()
    if (verbose) message(paste0(date(),": proceeding taxon ", objects_names[i]))

    # get all the sample names matching objects names
    # bundle all the biological replicates into one TE object

    #browser()

    ttl <- unlist(strsplit(objects_names[i], "_"))[1] #taxon title
    subttl <- unlist(strsplit(objects_names[i], "_"))[2] # subtaxon title
    #ttl <- lapply(names, function(x) unlist(strsplit(x, "_"))[1]) # taxon names
    #subttl <- lapply(names, function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names
    idx <- grep(objects_names[i],colnames(binding.score.df), ignore.case = T)
    names <- strsplit(colnames(binding.score.df)[idx],"_")
    repttl <- unlist(lapply(names, function(x) unlist(strsplit(x,"_"))[3])) # biological replicates title names


    # get gene names and lengths

    #gene_info = gene.info.df[grep(ttl,colnames(gene.info.df), ignore.case = T)]
    #tmp = apply(gene_info,1,function(x) unlist(strsplit(as.character(x),":")))

    #
    gene_names <- rownames(binding.score.df)  # gene names
    #gene_lengths <- as.integer(tmp[2,]) # gene lengths


    # foreach subtaxon
    bio_rep_n <- length(repttl) # biological replicates number
    omega <- NULL # omega estimated overdispersion parameter

    binding.score.raw <- apply(binding.score.df[idx], c(1,2), as.numeric)

    objects_counter = objects_counter + 1

    if (verbose) message(paste0(date(),": wrapping up into objects"))

    #browser()
    oneObject <- list(BindingScore.raw=binding.score.raw, taxon.name = ttl, tf.name = subttl,
                      gene.num = gene_n, gene.names = gene_names,
                      bioRep.num = bio_rep_n, bioRep.id = repttl, omega = omega
                      )

    class(oneObject) <- "taxonTF"

    taxonTF.objects[[objects_counter]] <- oneObject

    #browser()

    if (verbose) message(paste0(date(), ": ", objects_counter, " TF objects constructed"))

    if (verbose) cat("\n")

    if (!verbose) setTxtProgressBar(progbar, objects_counter/objects_number)


  }
  class(taxonTF.objects) <- "taxaTF"

  attr(taxonTF.objects, "taxa") <- unlist(taxon_names)
  attr(taxonTF.objects, "tf") <- unlist(tf_names)

  cat("\n")
  message(date(),": construction complete.")

  taxonTF.objects

}
