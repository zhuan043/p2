#' Title
#'
#' @param casesM aaa
#' @param ctrlsM bbb
#' @param optimal ccc
#' @param DeriveMed ddd
#' @param ClusteringMed eee
#' @param save_res fff
#'
#' @return
#'
#' @import Hmisc
#' @import mvtnorm
#' @import tmvtnorm
#' @import corpcor
#' @import WGCNA
#' @import doParallel
#' @import ggplot2
#' @import HDtest
#' @import parallel
#' @import foreach
#' @import cluster
#' @import factoextra
#' @import taRifx
#' @import multtest
#' @importFrom stats kmeans
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#'
#'
#'
#' @export
#'
#' @examples
#' library(multtest)
#' data(golub)
#' testcases1 <- DeriveClus(t(golub[1:500,1:27]),t(golub[1:500,28:38]),9,DeriveMed="Cases",ClusteringMed="Pam")


DeriveClus <- function(casesM, ctrlsM, optimal,DeriveMed,ClusteringMed,save_res =TRUE){

  print("Code error session 1")

  datC1 <- casesM
  datC2 <- ctrlsM

#  datC1 <- diffcoexcases
#  datC2 <- diffcoexctrls

  all.equal(colnames(datC1),colnames(datC1))

  Adj_datC1 <- adjacency(datC1,
                         selectCols = NULL,
                         type = "signed",
                         power = 1,
                         corFnc = "cor", corOptions = list(use = "p"),
                         weights = NULL,
                         distFnc = "dist", distOptions = "method = 'euclidean'",
                         weightArgNames = c("weights.x", "weights.y"))

  Adj_datC2 <- adjacency(datC2,
                         selectCols = NULL,
                         type = "signed",
                         power = 1,
                         corFnc = "cor", corOptions = list(use = "p"),
                         weights = NULL,
                         distFnc = "dist", distOptions = "method = 'euclidean'",
                         weightArgNames = c("weights.x", "weights.y"))

  TomdiffC1C2 <- TOMsimilarity(abs(Adj_datC1-Adj_datC2),
                               TOMType = "signed",
                               TOMDenom = "min",
                               suppressTOMForZeroAdjacencies = FALSE,
                               suppressNegativeTOM = FALSE,
                               useInternalMatrixAlgebra = FALSE,
                               verbose = 1,
                               indent = 0)
  dim(TomdiffC1C2)

  Tom_C1 <- TOMsimilarity(Adj_datC1,
                          TOMType = "signed",
                          TOMDenom = "min",
                          suppressTOMForZeroAdjacencies = FALSE,
                          suppressNegativeTOM = FALSE,
                          useInternalMatrixAlgebra = FALSE,
                          verbose = 1,
                          indent = 0)

  Tom_dist_C1 <- TOMdist(Adj_datC1,
                         TOMType = "signed",
                         TOMDenom = "min",
                         suppressTOMForZeroAdjacencies = FALSE,
                         suppressNegativeTOM = FALSE,
                         useInternalMatrixAlgebra = FALSE,
                         verbose = 1,
                         indent = 0)


  Tom_C2 <- TOMsimilarity(Adj_datC2,
                          TOMType = "signed",
                          TOMDenom = "min",
                          suppressTOMForZeroAdjacencies = FALSE,
                          suppressNegativeTOM = FALSE,
                          useInternalMatrixAlgebra = FALSE,
                          verbose = 1,
                          indent = 0)

  Tom_dist_C2 <- TOMdist(Adj_datC2,
                         TOMType = "signed",
                         TOMDenom = "min",
                         suppressTOMForZeroAdjacencies = FALSE,
                         suppressNegativeTOM = FALSE,
                         useInternalMatrixAlgebra = FALSE,
                         verbose = 1,
                         indent = 0)


  Tom_dist_diffC1C2 <- TOMdist(abs(Adj_datC1-Adj_datC2),
                               TOMType = "signed",
                               TOMDenom = "min",
                               suppressTOMForZeroAdjacencies = FALSE,
                               suppressNegativeTOM = FALSE,
                               useInternalMatrixAlgebra = FALSE,
                               verbose = 1,
                               indent = 0)

  colnames(Tom_C1) <- colnames(datC1)
  rownames(Tom_C1) <- colnames(datC1)





  # optimal <- 5

  print("Code error session 2")

  if (DeriveMed=="Cases")
  {
#    clusnum <-fviz_nbclust(Tom_dist_C1, cluster::pam, method = "gap_stat", nboot = 50, print.summary = TRUE)+
#    clusnum <- fviz_nbclust(Tom_dist_C1, cluster::pam, method = "gap_stat", k.max = 20, nboot = 50, print.summary = TRUE)+
 #     labs(subtitle = "Gap statistic method")

 #   clusnum <- fviz_nbclust(Tom_dist_C1, cluster::pam, method = "gap_stat", k.max = 50, nboot = 50, print.summary = TRUE)+
 #     labs(subtitle = "Gap statistic method")

#    pam.res <- pam(Tom_dist_C1, optimal, diss = TRUE)
    if (ClusteringMed=="kmeans"|ClusteringMed=="Kmeans")

    {pam.res <- kmeans(t(casesM),optimal,iter.max = 100)}

    if (ClusteringMed =="Pam"|ClusteringMed =="pam")
    { pam.res <- pam(t(casesM), optimal, diss = FALSE)}

    if (ClusteringMed =="hcut"|ClusteringMed =="Hcut")

    { pam.res <- hcut(t(casesM), optimal, diss = FALSE)}

    addClu <- cbind(Tom_dist_C1, cluster = pam.res$cluster)

    all.equal(colnames(datC1),colnames(addClu[,1:(dim(addClu)[2]-1)]))
    all.equal(colnames(datC1),rownames(addClu[,1:(dim(addClu)[2]-1)]))

    datC1C2 <- cbind(t(datC1),t(datC2))

    all.equal(rownames(datC1C2),rownames(addClu[,1:(dim(addClu)[2]-1)]))
    Clu_datC1C2 <- cbind(datC1C2,addClu[,dim(addClu)[2]])
    Clu_datC1C2 <- as.data.frame(Clu_datC1C2)
    names(Clu_datC1C2)[length(names(Clu_datC1C2))]<-"cluster"


    storeModcases <- list()
    storeModctrls <- list()
  }

  if (DeriveMed=="Ctrls")
  {
  #  optimal <- 11
 #   pam.res <- pam(Tom_dist_C2, optimal,diss = TRUE)
    #kmeans()
    if (ClusteringMed=="kmeans"|ClusteringMed=="Kmeans")

    {pam.res <- kmeans(t(ctrlsM),optimal,iter.max = 100)}

    if (ClusteringMed =="Pam"|ClusteringMed =="pam")
    { pam.res <- pam(t(ctrlsM), optimal, diss = FALSE)}

    if (ClusteringMed =="hcut"|ClusteringMed =="Hcut")

    { pam.res <- hcut(t(ctrlsM), optimal, diss = FALSE)}

 #   clusnum <- fviz_nbclust(Tom_dist_C2, cluster::pam, method = "gap_stat", k.max = 20,nboot = 50, print.summary = TRUE)+
  #    labs(subtitle = "Gap statistic method")

    addClu <- cbind(Tom_dist_C2, cluster = pam.res$cluster)
#    View(as.data.frame(addClu))
 #   dim(as.data.frame(addClu))

    all.equal(colnames(datC2),colnames(addClu[,1:(dim(addClu)[2]-1)]))
    all.equal(colnames(datC2),rownames(addClu[,1:(dim(addClu)[2]-1)]))

    datC1C2 <- cbind(t(datC1),t(datC2))
#    dim(datC1C2)

    all.equal(rownames(datC1C2),rownames(addClu[,1:(dim(addClu)[2]-1)]))
    Clu_datC1C2 <- cbind(datC1C2,addClu[,dim(addClu)[2]])
    Clu_datC1C2 <- as.data.frame(Clu_datC1C2)

    names(Clu_datC1C2)[length(names(Clu_datC1C2))]<-"cluster"

    storeModcases <- list()
    storeModctrls <- list()
  }

  if (DeriveMed=="Diff")
  {
    if (ClusteringMed=="kmeans"|ClusteringMed=="Kmeans")

    {pam.res <- kmeans(Tom_dist_diffC1C2,optimal,iter.max = 100)}

    if (ClusteringMed =="Pam"|ClusteringMed =="pam")
    { pam.res <- pam(Tom_dist_diffC1C2, optimal, diss = TRUE)}

    if (ClusteringMed =="hcut"|ClusteringMed =="Hcut")

    { pam.res <- hcut(Tom_dist_diffC1C2, optimal, diss = TRUE)}

  #  clusnum <- fviz_nbclust(Tom_dist_diffC1C2, cluster::pam, method = "gap_stat", k.max = 20,nboot = 50, print.summary = TRUE)+
  #    labs(subtitle = "Gap statistic method")

    addClu <- cbind(Tom_dist_diffC1C2, cluster = pam.res$cluster)

    all.equal(colnames(TomdiffC1C2),colnames(addClu[,1:(dim(addClu)[2]-1)])) ## Need to check the code
    all.equal(colnames(TomdiffC1C2),rownames(addClu[,1:(dim(addClu)[2]-1)])) ## Need to check the code

    datC1C2 <- cbind(t(datC1),t(datC2))

    all.equal(rownames(datC1C2),rownames(addClu[,1:(dim(addClu)[2]-1)]))
    Clu_datC1C2 <- cbind(datC1C2,addClu[,dim(addClu)[2]])
    Clu_datC1C2 <- as.data.frame(Clu_datC1C2)
    names(Clu_datC1C2)[length(names(Clu_datC1C2))]<-"cluster"


    storeModcases <- list()
    storeModctrls <- list()
  }


  for (i in 1:optimal) {
    # assign(paste("module",i,sep = "_"),Clu_datC1C2[Clu_datC1C2$V81==i,])
    storeModcases[[i]] <-assign(paste("cases",i,"module",sep = "_"),t(Clu_datC1C2[Clu_datC1C2$cluster==i,c(1:dim(casesM)[1])]))
    storeModctrls[[i]] <- assign(paste("ctrls",i,"module",sep = "_"),
                                 t(Clu_datC1C2[Clu_datC1C2$cluster==i,c((dim(casesM)[1]+1):(dim(casesM)[1]+dim(ctrlsM)[1]))]))
    #colnames(paste("cases",i,"module",sep = "_")) <- NULL
    # assign(colnames(paste("cases",i,"module",sep = "_")),NULL)

  }

  if (save_res == TRUE)
  {
  assign("Modcases",storeModcases,envir = .GlobalEnv)
  assign("Modctrls",storeModctrls,envir = .GlobalEnv)
#  assign("Gap_stat_clu",clusnum,envir = .GlobalEnv)
}
#  storeclus <- list()
  invisible(return(list("Modcases"=storeModcases,"Modctrls"=storeModctrls
                        #,"Gap_stat_clu"=clusnum
                        )))
}
