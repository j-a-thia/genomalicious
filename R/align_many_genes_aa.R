#' Alignment of many genes (amino acid sequences)
#'
#' Reads in multiple FASTA files of different genes, performs an alignment,
#' and tests the best evolutionary model. Specifically for amino acid sequences.
#'
#' @param fasta.files Character, vector of paths to FASTA files. Each one
#' containing sequences for a specific gene for multiple taxa.
#'
#' @param method Character, one of "ClustalW", "ClustalOmega", "Muscle".
#' Is passed to the function \code{msa::msa(... method=method)}.
#'
#' @param model.model Character, evolutionary models to test. Any one, or all of
#' the possible model arguments in \code{phangorn::modelTest}. Passed to function
#' \code{phangorn::modelTest(... model=model.model)}. Default is \code{NULL},
#' in which case, no models will be tested.
#'
#' @param model.G Logical, whether a gamma distribution of rates should be
#' tested in evolutionary models. Passed to function
#' \code{phangorn::modelTest(... G=model.G)}. Default is TRUE.
#'
#' @param model.I Logical, whether invariant rates should be tested in
#' evolutionary models. Passed to function
#' \code{phangorn::modelTest(... I=model.I)}. Default is TRUE.
#'
#' @param cores Integer, number of parallel jobs to run. Default is 1.
#'
#' @param gene.names Character, an optional vector of gene names that match
#' the FASTA files listed in \code{fasta.files}. Used to provide name handles
#' for each gene in the final output list. Default is NULL.
#'
#' @details The top evolutionary models per gene are identified using
#' AIC criteria. First, the model with the smallest AIC value is identified,
#' then the delta AIC is calculated for all models as focal model AIC minus the
#' AIC of the best model. Those models with delta AIC <= 10 are cosnidered
#' the top models and non-differentiable.
#'
#' @return Returns a list. Each index is one of the genes from the vector of
#' FASTA files specified in \code{fasta.files}. For each indexed gene, there
#' are additional slots:
#' \enumerate{
#'    \item \code{$gene} = AAStringSet, the imported FASTA sequences.
#'    \item \code{$align} = msa, the multiple sequence alignment for the gene.
#'    \item \code{$model.test} = Data.table, the output of
#'       \code{phangorn::modelTest}.
#'    \item \code{$model.top} = Data.table, the top evolutionary models.
#' }
#'
#' @examples
#' library(genomalicious)
#'
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Path to the demo FASTA file
#' fasta <- paste0(genomaliciousExtData, '/data_COI_aa.fasta')
#'
#' # Multi sequence alignmnet of demo COI data, without model test
#' aln <- align_many_genes_aa(fasta, gene.names='COI')
#'
#' str(aln)
#' names(aln)
#' aln$COI
#'
#' # Same again, but with model test
#' aln.test <- align_many_genes_aa(
#'    fasta.files=fasta, gene.names='COI', model.I=TRUE, model.G=TRUE,
#'    model.model='all'
#'    )
#'
#' aln.test$COI$model.test
#' aln.test$COI$model.top
#'
#'
#' @export
align_many_genes_aa <- function(
    fasta.files, method='ClustalW',
    model.model=NULL,
    model.G=TRUE, model.I=TRUE, cores=1, gene.names=NULL){

  require(ape)
  require(doParallel)
  require(data.table)
  require(phangorn)
  require(msa)
  require(tidyverse)

  if(cores==1){
    alignList <- lapply(fasta.files, function(fa){
      # Read in gene FASTA
      gene <- readAAStringSet(fa)

      # Align
      align <- msa(gene, method = method, order='input')
      align.phydat <- as.phyDat(align)

      # Test evolutionary models
      if(!is.null(model.model)){
        # Fit models
        model.test <- phangorn::modelTest(
          align.phydat,
          model=model.model, G=model.G, I=model.I) %>%
          as.data.table()

        # Best model
        best.mod <- model.test %>%
          as.data.table() %>%
          setorder(AIC) %>%
          .[1,'Model'] %>%
          unlist()

        # AIC of best model
        best.aic <- model.test[Model==best.mod]$AIC

        # Top models based on deltaAIC <= 10
        model.top <- model.test %>%
          as.data.table() %>%
          .[, deltaAIC:=AIC-best.aic] %>%
          .[deltaAIC <= 10] %>%
          setorder(., AIC)
      } else{
        model.test <- NULL
        model.top <- NULL
      }

      # Output
      list(gene=gene, align=align, model.test=model.test, model.top=model.top) %>%
        return()
    })
  } else if(cores>1){
    # Set up cluster
    my.clust <- makeCluster(cores)
    registerDoParallel(my.clust)

    # Iterate gene FASTA files
    alignList <- foreach(fa=fasta.files) %dopar% {
      require(ape)
      require(doParallel)
      require(data.table)
      require(phangorn)
      require(msa)
      require(tidyverse)

      # Read in gene FASTA
      gene <- readDNAStringSet(fa)

      # Align
      align <- msa(gene, method = method, order='input')

      # Test evolutionary models
      model.test <- modelTest(
        align %>% as.phyDat(),
        model=model.model, G=model.G, I=model.I) %>%
        as.data.table()

      # Best model
      best.mod <- model.test %>%
        as.data.table() %>%
        setorder(AIC) %>%
        .[1,'Model'] %>%
        unlist()

      # AIC of best model
      best.aic <- model.test[Model==best.mod]$AIC

      # Top models based on deltaAIC <= 10
      model.top <- model.test %>%
        as.data.table() %>%
        .[, deltaAIC:=AIC-best.aic] %>%
        .[deltaAIC <= 10] %>%
        setorder(., AIC)

      # Output
      list(gene=gene, align=align, model.test=model.test, model.top=model.top)
    }

    stopCluster(my.clust)
  }

  # Add names if specified.
  if(is.null(gene.names)==FALSE){
    names(alignList) <- gene.names
  }

  # Return
  return(alignList) %>% return()
}
