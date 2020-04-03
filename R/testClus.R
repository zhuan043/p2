
#' Title
#'
#' @param Derived_Modules aaaa
#' @param permnum  bbbbb
#' @param compare_med ccccc
#' @param adjpower ddddd
#' @param type eeeee
#' @param save_res fffff
#'
#' @return
#'
#' @import Hmisc
#' @import mvtnorm
#' @import tmvtnorm
#' @import corpcor
#' @import WGCNA
#' @import ggplot2
#' @import HDtest
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
#' @examples
#' library(multtest)
#' data(golub)
#' testcases1 <- DeriveClus(t(golub[1:100,1:27]),t(golub[1:100,28:38]),7,DeriveMed="Cases",ClusteringMed="Pam")
#' trial222 <- testClus(testcases1,100,compare_med="adjvstom",type = "signed",adjpower = 0.5)
#'
#'
testClus <- function(Derived_Modules,permnum,compare_med = "adjvscorr", adjpower=1,type = "signed",save_res =TRUE){

  cores <- 2
  progress <- function(n, tag) if(tag%%100==0){print(tag)}
  opts <- list(progress = progress)
  cl <- makeCluster(cores)
  registerDoSNOW(cl)


  print("code error session 1")

  combres <- list()
  T1storeperm <- list()
  T1storeobs <- list()

#  length(as.list(testcases1$Modcases))
  Modcases <- as.list(Derived_Modules$Modcases)
  Modctrls <- as.list(Derived_Modules$Modctrls)

  optimal <- length(Modcases) # Because cluster 9 only had 1 gene

  storemodres <- matrix(rep(0,optimal*16),nrow = 16, ncol = optimal)

  # wilcoxin tom
  store <- rep(0,optimal)

  # wilcoxin adj
  storeAdj <- rep(0,optimal)

  # t-test tom
  storettest <- rep(0,optimal)

  # t-test adj
  storettestadj <- rep(0,optimal)

  # Mad adj
  storeMadADJ <- rep(0,optimal)

  # Mad tom
  storeMadTOM <- rep(0,optimal)

  # Dispersion Method
  storeDisperInd <- rep(0,optimal)

  storeDisperInd_tom <- rep(0,optimal)

  storecubAdj <- rep(0,optimal)

  storecubTom <- rep(0,optimal)
  ##### Newly Added
  storeMadRankingAdj <- rep(0,optimal)
  storeMadRankingTom <- rep(0,optimal)
  storeDisperRankingAdj <- rep(0,optimal)
  storeDisperRankingTom <- rep(0,optimal)
  storeQuadtestRankingAdj <- rep(0,optimal)
  storeQuadtestRankingTom <- rep(0,optimal)

  if (compare_med == "adj vs corr")
  {

  for (k in 1:optimal) {

    if (dim(as.data.frame(Modcases[k]))[2] > 1)
        {
    print("code error session 2")
    #    k=1
    print(k)

    obs_adjMat_cases <- adjacency(as.data.frame(Modcases[k]),
                                  selectCols = NULL,
                                  type = type,
                                  power = adjpower,
                                  corFnc = "cor", corOptions = list(use = "p"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  weightArgNames = c("weights.x", "weights.y"))


    # check if adj matrix is valid
    checkAdjMat(obs_adjMat_cases, min = 0, max = 1)


    # calcualte adj matrix for controls
    obs_adjMat_ctrls <- adjacency(as.data.frame(Modctrls[k]),
                                  selectCols = NULL,
                                  type = type,
                                  power = adjpower,
                                  corFnc = "cor", corOptions = list(use = "p"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  weightArgNames = c("weights.x", "weights.y"))

    # check adj matrix for controls
    checkAdjMat(obs_adjMat_ctrls, min = 0, max = 1)

    obs_TOM_cases <- cor(as.data.frame(Modcases[k]),method = "pearson")
    obs_TOM_cases <- obs_TOM_cases[lower.tri(obs_TOM_cases)]

    obs_TOM_ctrls <- cor(as.data.frame(Modctrls[k]), method = "pearson")
    obs_TOM_ctrls <- obs_TOM_ctrls[lower.tri(obs_TOM_ctrls)]

    # use lower triangle of observed adj cases matrix
    obs_adj_cases <- obs_adjMat_cases[lower.tri(obs_adjMat_cases)]

    # use lower triangle of observed adj controls matrix
    obs_adj_ctrls <- obs_adjMat_ctrls[lower.tri(obs_adjMat_ctrls)]

    # observed wilconxin test results for tom matrix
    result <- wilcox.test(obs_TOM_cases,obs_TOM_ctrls, paired = TRUE)

    # observed wilconxin test results for adj matrix
    adjresult <-  wilcox.test(obs_adj_cases,obs_adj_ctrls,paired = TRUE)

    # observed t test results for tom matrix
    ttestresult <- t.test(obs_TOM_cases,obs_TOM_ctrls, paired = TRUE)

    # observed t test results for adj matrix
    adjttestresult <- t.test(obs_adj_cases,obs_adj_ctrls,paired = TRUE)


    # observed mad test results for adj matrix
    madADJ <- mean(abs(obs_adj_cases-obs_adj_ctrls))

    # observed mad test results for tom matrix
    madTOM <- mean(abs(obs_TOM_cases-obs_TOM_ctrls))



    # calculate difference between cases correlation vector and controls correlation vector
    cor_diff <- obs_adj_cases - obs_adj_ctrls

    # observed dispersion index
    obs_dis_ind <- sqrt(mean(cor_diff^2))

    # calculate difference between cases correlation vector and controls correlation vector
    cor_diff_tom <- obs_TOM_cases - obs_TOM_ctrls

    # observed dispersion index
    obs_dis_ind_tom <- sqrt(mean(cor_diff_tom^2))

    obs_cubtest_Adj <- mean(cor_diff^4)

    obs_cubtest_Tom <- mean(cor_diff_tom^4)

    #### Newly Added #####
    ranking_adj <- rank(abs(cor_diff))
    ranking_tom <- rank(abs(cor_diff_tom))


    obs_mad_ranking_adj <- mean(abs(obs_adj_cases-obs_adj_ctrls)*ranking_adj)
    obs_mad_ranking_tom <- mean(abs(obs_TOM_cases-obs_TOM_ctrls)*ranking_tom)

    obs_dis_ranking_adj <- sqrt(mean((cor_diff^2)*ranking_adj))
    obs_dis_ranking_tom <- sqrt(mean((cor_diff_tom^2)*ranking_tom))

    obs_quadtest_Ranking_Adj <- mean((cor_diff^4)*ranking_adj)
    obs_quadtest_Ranking_Tom <- mean((cor_diff_tom^4)*ranking_tom)

    # register parallel computing with number of threads (for intel 4th gen and after, one core has two threads)
    #   registerDoParallel(numCores)

    #    permnum <- 1000
    #    G <- dim(eval(as.symbol(paste("cases",k,"module",sep = "_"))))[2]
    G <- dim(as.data.frame(Modcases[k]))[2]

    #colnames(cases_1_module) <- NULL
    #rownames(cases_1_module) <- NULL
    #colnames(ctrls_1_module) <- NULL
    #rownames(ctrls_1_module) <- NULL

    #    tempcases <- eval(as.symbol(paste("cases",k,"module",sep = "_")))
    tempcases <- as.data.frame(Modcases[k])

    #    tempctrls <- eval(as.symbol(paste("ctrls",k,"module",sep = "_")))

    tempctrls <- as.data.frame(Modctrls[k])

    T1obstat <- rbind(result$statistic,ttestresult$statistic,adjresult$statistic,
                      adjttestresult$statistic,madADJ,madTOM,obs_dis_ind,obs_dis_ind_tom,
                      obs_cubtest_Adj,obs_cubtest_Tom,obs_mad_ranking_adj,obs_mad_ranking_tom,obs_dis_ranking_adj,
                      obs_dis_ranking_tom,obs_quadtest_Ranking_Adj,obs_quadtest_Ranking_Tom)

    rownames(T1obstat) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffcubicAdj","DiffcubicTom",
                            "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

    T1storeobs[[k]] <- T1obstat

    pararesults <- foreach(j = 1:permnum,.combine=rbind,.options.snow = opts)%dopar%{
      # set.seed(123)

      library(WGCNA)

      # restore cases module as data frame
      casesPermDF <- as.data.frame(tempcases)
 #     colnames(casesPermDF) <- NULL
 #     rownames(casesPermDF) <- NULL
      # create labels for cases module
      casesPermDF$ind <- "cases"

      # restore controls module as data frame
      ctrlsPermDF <- as.data.frame(tempctrls)
#      colnames(ctrlsPermDF) <- NULL
#      rownames(ctrlsPermDF) <- NULL
      # create lables for controls module
      ctrlsPermDF$ind <- "ctrls"

      # rbind cases module and controls module into a single data frame
      totalDF <- rbind(casesPermDF,ctrlsPermDF)

      # restore this data frame
      newdf <- totalDF

      # shuffle the lables to resample cases and controls
      newdf$ind <- sample(newdf$ind,length(newdf$ind),FALSE)

      # new cases data frame
      newcasesdf <- newdf[which(newdf$ind=="cases"),1:G]

      # new controls data frame
      newctrlsdf <- newdf[which(newdf$ind=="ctrls"),1:G]


      # new adj matrix for cases after the resample
      adjMat_cases <- adjacency(newcasesdf,
                                selectCols = NULL,
                                type = type,
                                power = adjpower,
                                corFnc = "cor", corOptions = list(use = "p"),
                                weights = NULL,
                                distFnc = "dist", distOptions = "method = 'euclidean'",
                                weightArgNames = c("weights.x", "weights.y"))

      # check new adj matrix for cases
      checkAdjMat(adjMat_cases, min = 0, max = 1)


      # new adj matrix for controls after the resample
      adjMat_ctrls <- adjacency(newctrlsdf,
                                selectCols = NULL,
                                type = type,
                                power = adjpower,
                                corFnc = "cor", corOptions = list(use = "p"),
                                weights = NULL,
                                distFnc = "dist", distOptions = "method = 'euclidean'",
                                weightArgNames = c("weights.x", "weights.y"))

      # check new adj matrix for controls
      checkAdjMat(adjMat_ctrls, min = 0, max = 1)


      Perm_cases <- cor(newcasesdf,method = "pearson")
      Perm_cases <- Perm_cases[lower.tri(Perm_cases)]

      Perm_ctrls <- cor(newctrlsdf, method = "pearson")
      Perm_ctrls <- Perm_ctrls[lower.tri(Perm_ctrls)]

      # use lower triangle of permuted tom matrix for cases
  #    Perm_cases <- Tom_cases[lower.tri(Tom_cases)]

      # use lower triangle of permuted tom matrix for controls
  #    Perm_ctrls <- Tom_ctrls[lower.tri(Tom_ctrls)]

      # use lower triangle of permuted adj matrix for cases
      Perm_adj_cases <- adjMat_cases[lower.tri(adjMat_cases)]

      # use lower triangle of permuted adj matrix for controls
      Perm_adj_ctrls <- adjMat_ctrls[lower.tri(adjMat_ctrls)]



      # permuted wilconxin test results for tom matrix
      newresults <- wilcox.test(Perm_cases,Perm_ctrls, paired = TRUE)

      # store wilconxin test statistics
      StatPerm <- newresults$statistic

      # permuted t test results for tom matrix
      newtresults <- t.test(Perm_cases,Perm_ctrls, paired = TRUE)

      # store t test statistics
      ttestPerm <- newtresults$statistic

      # permuted wilconxin test results for adj matrix
      newadjresults <- wilcox.test(Perm_adj_cases,Perm_adj_ctrls,paired = TRUE)

      # store wilconxin test statistics
      AdjPerm <- newadjresults$statistic

      # permuted t test results for adj matrix
      newtresultsadj <- t.test(Perm_adj_cases,Perm_adj_ctrls,paired = TRUE)

      # store t test statistics
      ttestadjPerm <- newtresultsadj$statistic

      # permuted mad test for adj matrix
      PermmadADJ <- mean(abs(Perm_adj_cases-Perm_adj_ctrls))

      # permuted mad test for tom matrix
      PermmadTOM <- mean(abs(Perm_cases-Perm_ctrls))


      # permuted difference between cases and controls correlation vector
      Perm_cor_diff <- Perm_adj_cases - Perm_adj_ctrls

      # permuted dispersion index
      Perm_dis_ind <- sqrt(mean(Perm_cor_diff^2))

      Perm_cor_diff_tom <- Perm_cases - Perm_ctrls

      Perm_dis_ind_tom <- sqrt(mean(Perm_cor_diff_tom^2))

      Perm_cubtest_Adj <- mean(Perm_cor_diff^4)

      Perm_cubtest_Tom <- mean(Perm_cor_diff_tom^4)


      Perm_ranking_Adj <- rank(abs(Perm_cor_diff))

      Perm_ranking_Tom <- rank(abs(Perm_cor_diff_tom))


      Perm_mad_ranking_adj <- mean(abs(Perm_adj_cases-Perm_adj_ctrls)*Perm_ranking_Adj)
      Perm_mad_ranking_tom <- mean(abs(Perm_cases-Perm_ctrls)*Perm_ranking_Tom)

      Perm_dis_ranking_adj <- sqrt(mean((Perm_cor_diff^2)*Perm_ranking_Adj))
      Perm_dis_ranking_tom <- sqrt(mean((Perm_cor_diff_tom^2)*Perm_ranking_Tom))

      Perm_quadtest_ranking_adj <- mean((Perm_cor_diff^4)*Perm_ranking_Adj)
      Perm_quadtest_ranking_tom <- mean((Perm_cor_diff_tom^4)*Perm_ranking_Tom)

      # create a vector to store all above test statistics
      intermres <- c(StatPerm,ttestPerm,AdjPerm,ttestadjPerm,PermmadADJ,PermmadTOM,Perm_dis_ind,Perm_dis_ind_tom,Perm_cubtest_Adj,Perm_cubtest_Tom,
                     Perm_mad_ranking_adj,Perm_mad_ranking_tom,Perm_dis_ranking_adj,Perm_dis_ranking_tom,Perm_quadtest_ranking_adj,
                     Perm_quadtest_ranking_tom)

      # name that vector
      names(intermres) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffcubicAdj","DiffcubicTom",
                            "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

      return(intermres)
    }

    # store as a data frame
    pararesults <- as.data.frame(pararesults)

    T1storeperm [[k]] <- pararesults

    # calculate p value for wilcoxin tom
    pvaluePerm <- ((sum(pararesults$wTOM <= result$statistic))+1)/(permnum+1)
    store[k] <- pvaluePerm

    # calculate p value for t test tom
    pvaluePermttest <- ((sum(abs(pararesults$tTOM) >= abs(ttestresult$statistic)))+1)/(permnum+1)
    storettest[k] <- pvaluePermttest

    # calculate p values for wilcoxin adj
    pvaluePermAdj <- ((sum(pararesults$wADJ <= adjresult$statistic))+1)/(permnum+1)
    storeAdj[k] <- pvaluePermAdj

    # calculate p values for t test adj
    pvaluePermttestadj <- ((sum(abs(pararesults$tADJ) >= abs(adjttestresult$statistic)))+1)/(permnum+1)
    storettestadj[k] <- pvaluePermttestadj

    # calculate p values for mad adj
    pvalueMadADJ <- ((sum(pararesults$MadADJ >= madADJ))+1)/(permnum+1)
    storeMadADJ[k] <- pvalueMadADJ

    # calculate p values for mad tom
    pvalueMadTOM <- ((sum(pararesults$MadTOM >= madTOM))+1)/(permnum+1)
    storeMadTOM[k] <- pvalueMadTOM

    # calculate p values for Dispersion Index
    pvalueDisper <- ((sum(pararesults$DispersionIndexAdj >= obs_dis_ind))+1)/(permnum+1)
    storeDisperInd[k] <- pvalueDisper

    # calculate p values for Dispersion Index
    pvalueDisper_tom <- ((sum(pararesults$DispersionIndexTom >= obs_dis_ind_tom))+1)/(permnum+1)
    storeDisperInd_tom[k] <- pvalueDisper_tom

    pvaluequadAdj <- ((sum(pararesults$DiffcubicAdj >= obs_cubtest_Adj))+1)/(permnum+1)
    storecubAdj[k] <- pvaluequadAdj

    pvaluequadTom <- ((sum(pararesults$DiffcubicTom >= obs_cubtest_Tom))+1)/(permnum+1)
    storecubTom[k] <- pvaluequadTom

    ########### Newly added #####
    pvalue_mad_ranking_adj <- ((sum(pararesults$MadRankingAdj >= obs_mad_ranking_adj))+1)/(permnum+1)
    storeMadRankingAdj[k] <- pvalue_mad_ranking_adj

    pvalue_mad_ranking_tom <- ((sum(pararesults$MadRankingTom >= obs_mad_ranking_tom))+1)/(permnum+1)
    storeMadRankingTom[k] <- pvalue_mad_ranking_tom

    pvalue_dispersion_ranking_adj <- ((sum(pararesults$DispersionRankingAdj >= obs_dis_ranking_adj))+1)/(permnum+1)
    storeDisperRankingAdj[k] <- pvalue_dispersion_ranking_adj

    pvalue_dispersion_ranking_tom <- ((sum(pararesults$DispersionRankingTom >= obs_dis_ranking_tom))+1)/(permnum+1)
    storeDisperRankingTom[k] <- pvalue_dispersion_ranking_tom

    pvalue_quadtest_ranking_adj <- ((sum(pararesults$QuadtestRankingAdj >= obs_quadtest_Ranking_Adj))+1)/(permnum+1)
    storeQuadtestRankingAdj[k] <- pvalue_quadtest_ranking_adj

    pvalue_quadtest_ranking_tom <- ((sum(pararesults$QuandtestRankingTom >= obs_quadtest_Ranking_Tom))+1)/(permnum+1)
    storeQuadtestRankingTom[k] <- pvalue_quadtest_ranking_tom

    testRes <-rbind(pvaluePermAdj,pvaluePerm,pvaluePermttestadj,pvaluePermttest,pvalueMadADJ,pvalueMadTOM,pvalueDisper,pvalueDisper_tom,
                    pvaluequadAdj,pvaluequadTom,pvalue_mad_ranking_adj,pvalue_mad_ranking_tom,pvalue_dispersion_ranking_adj,pvalue_dispersion_ranking_tom,
                    pvalue_quadtest_ranking_adj,pvalue_quadtest_ranking_tom)

    storemodres[,k] <- rbind(storeAdj[k],store[k],storettestadj[k],storettest[k],storeMadADJ[k],storeMadTOM[k],storeDisperInd[k],storeDisperInd_tom[k],
                            storecubAdj[k],storecubTom[k],storeMadRankingAdj[k],storeMadRankingTom[k],storeDisperRankingAdj[k],storeDisperRankingTom[k],
                            storeQuadtestRankingAdj[k],storeQuadtestRankingTom[k])

    testRes <- as.data.frame(testRes)

 #   Resname <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","QuadtestAdj","QuadtestTom",
 #                 "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

    Resname <- c(paste("Wilcoxon_adj",adjpower,sep = "_"),"Wilcoxon_corr", paste("Ttest_adj",adjpower,sep = "_"),"Ttest_corr",paste("Mad_adj",adjpower,sep = "_"),"Mad_corr",
                 paste("DispersionIndex_adj",adjpower,sep = "_"),"DispersionIndex_corr",paste("QuarticTest_adj",adjpower,sep = "_"),"QuarticTest_corr",
                 paste("Mad_Ranking_adj",adjpower,sep = "_"),"Mad_Ranking_corr",paste("Dispersion_Ranking_adj",adjpower,sep = "_"),"Dispersion_Ranking_corr",
                 paste("QuarticTest_Ranking_adj",adjpower,sep = "_"),"Quandtest_Ranking_corr")

    rownames(testRes) <- Resname
    colnames(testRes) <- paste("Module",k,"stats",sep = "_")

    combres[[k]] <- testRes
    }
   if (dim(as.data.frame(Modcases[k]))[2] == 1)
   {
     print("code error session 3")
     warning(paste("Module",k,"only has one gene",sep = "_"))
     testRes <- NA
     combres[[k]] <- testRes
   }
  }

  finalRes <- do.call(cbind,combres)
  }

  if (compare_med == "adjvstom")

  {
    for (k in 1:optimal) {

    if (dim(as.data.frame(Modcases[k]))[2] > 1)
    {
      print("code error session 2")
      #    k=1
      print(k)

      obs_adjMat_cases <- adjacency(as.data.frame(Modcases[k]),
                                    selectCols = NULL,
                                    type = type,
                                    power = adjpower,
                                    corFnc = "cor", corOptions = list(use = "p"),
                                    weights = NULL,
                                    distFnc = "dist", distOptions = "method = 'euclidean'",
                                    weightArgNames = c("weights.x", "weights.y"))


      # check if adj matrix is valid
      checkAdjMat(obs_adjMat_cases, min = 0, max = 1)

      # calculate tom matrix for cases
      obs_Tom_cases <- TOMsimilarity(obs_adjMat_cases,
                                     TOMType = type,
                                     TOMDenom = "min",
                                     suppressTOMForZeroAdjacencies = FALSE,
                                     suppressNegativeTOM = FALSE,
                                     useInternalMatrixAlgebra = FALSE,
                                     verbose = 1,
                                     indent = 0)

      # calcualte adj matrix for controls
      obs_adjMat_ctrls <- adjacency(as.data.frame(Modctrls[k]),
                                    selectCols = NULL,
                                    type = type,
                                    power = adjpower,
                                    corFnc = "cor", corOptions = list(use = "p"),
                                    weights = NULL,
                                    distFnc = "dist", distOptions = "method = 'euclidean'",
                                    weightArgNames = c("weights.x", "weights.y"))

      # check adj matrix for controls
      checkAdjMat(obs_adjMat_ctrls, min = 0, max = 1)

      # calculate tom matrix for controls
      obs_Tom_ctrls <- TOMsimilarity(obs_adjMat_ctrls,
                                     TOMType = type,
                                     TOMDenom = "min",
                                     suppressTOMForZeroAdjacencies = FALSE,
                                     suppressNegativeTOM = FALSE,
                                     useInternalMatrixAlgebra = FALSE,
                                     verbose = 1,
                                     indent = 0)

      # use lower triangle of observed tom cases matrix
      obs_TOM_cases <- obs_Tom_cases[lower.tri(obs_Tom_cases)]

      # use lower triangle of observed tom controls matrix
      obs_TOM_ctrls <- obs_Tom_ctrls[lower.tri(obs_Tom_ctrls)]

      # use lower triangle of observed adj cases matrix
      obs_adj_cases <- obs_adjMat_cases[lower.tri(obs_adjMat_cases)]

      # use lower triangle of observed adj controls matrix
      obs_adj_ctrls <- obs_adjMat_ctrls[lower.tri(obs_adjMat_ctrls)]

      # observed wilconxin test results for tom matrix
      result <- wilcox.test(obs_TOM_cases,obs_TOM_ctrls, paired = TRUE)

      # observed wilconxin test results for adj matrix
      adjresult <-  wilcox.test(obs_adj_cases,obs_adj_ctrls,paired = TRUE)

      # observed t test results for tom matrix
      ttestresult <- t.test(obs_TOM_cases,obs_TOM_ctrls, paired = TRUE)

      # observed t test results for adj matrix
      adjttestresult <- t.test(obs_adj_cases,obs_adj_ctrls,paired = TRUE)


      # observed mad test results for adj matrix
      madADJ <- mean(abs(obs_adj_cases-obs_adj_ctrls))

      # observed mad test results for tom matrix
      madTOM <- mean(abs(obs_TOM_cases-obs_TOM_ctrls))



      # calculate difference between cases correlation vector and controls correlation vector
      cor_diff <- obs_adj_cases - obs_adj_ctrls

      # observed dispersion index
      obs_dis_ind <- sqrt(mean(cor_diff^2))

      # calculate difference between cases correlation vector and controls correlation vector
      cor_diff_tom <- obs_TOM_cases - obs_TOM_ctrls

      # observed dispersion index
      obs_dis_ind_tom <- sqrt(mean(cor_diff_tom^2))

      obs_cubtest_Adj <- mean(cor_diff^4)

      obs_cubtest_Tom <- mean(cor_diff_tom^4)

      #### Newly Added #####
      ranking_adj <- rank(abs(cor_diff))
      ranking_tom <- rank(abs(cor_diff_tom))


      obs_mad_ranking_adj <- mean(abs(obs_adj_cases-obs_adj_ctrls)*ranking_adj)
      obs_mad_ranking_tom <- mean(abs(obs_TOM_cases-obs_TOM_ctrls)*ranking_tom)

      obs_dis_ranking_adj <- sqrt(mean((cor_diff^2)*ranking_adj))
      obs_dis_ranking_tom <- sqrt(mean((cor_diff_tom^2)*ranking_tom))

      obs_quadtest_Ranking_Adj <- mean((cor_diff^4)*ranking_adj)
      obs_quadtest_Ranking_Tom <- mean((cor_diff_tom^4)*ranking_tom)

      # register parallel computing with number of threads (for intel 4th gen and after, one core has two threads)
      #   registerDoParallel(numCores)

      #    permnum <- 1000
      #    G <- dim(eval(as.symbol(paste("cases",k,"module",sep = "_"))))[2]
      G <- dim(as.data.frame(Modcases[k]))[2]

      #colnames(cases_1_module) <- NULL
      #rownames(cases_1_module) <- NULL
      #colnames(ctrls_1_module) <- NULL
      #rownames(ctrls_1_module) <- NULL

      #    tempcases <- eval(as.symbol(paste("cases",k,"module",sep = "_")))
      tempcases <- as.data.frame(Modcases[k])

      #    tempctrls <- eval(as.symbol(paste("ctrls",k,"module",sep = "_")))

      tempctrls <- as.data.frame(Modctrls[k])

      T1obstat <- rbind(result$statistic,ttestresult$statistic,adjresult$statistic,
                        adjttestresult$statistic,madADJ,madTOM,obs_dis_ind,obs_dis_ind_tom,
                        obs_cubtest_Adj,obs_cubtest_Tom,obs_mad_ranking_adj,obs_mad_ranking_tom,obs_dis_ranking_adj,
                        obs_dis_ranking_tom,obs_quadtest_Ranking_Adj,obs_quadtest_Ranking_Tom)

      rownames(T1obstat) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffcubicAdj","DiffcubicTom",
                              "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

      T1storeobs[[k]] <- T1obstat

      pararesults <- foreach(j = 1:permnum,.combine=rbind,.options.snow = opts)%dopar%{
        # set.seed(123)

        library(WGCNA)

        # restore cases module as data frame
        casesPermDF <- as.data.frame(tempcases)
        #     colnames(casesPermDF) <- NULL
        #     rownames(casesPermDF) <- NULL
        # create labels for cases module
        casesPermDF$ind <- "cases"

        # restore controls module as data frame
        ctrlsPermDF <- as.data.frame(tempctrls)
        #      colnames(ctrlsPermDF) <- NULL
        #      rownames(ctrlsPermDF) <- NULL
        # create lables for controls module
        ctrlsPermDF$ind <- "ctrls"

        # rbind cases module and controls module into a single data frame
        totalDF <- rbind(casesPermDF,ctrlsPermDF)

        # restore this data frame
        newdf <- totalDF

        # shuffle the lables to resample cases and controls
        newdf$ind <- sample(newdf$ind,length(newdf$ind),FALSE)

        # new cases data frame
        newcasesdf <- newdf[which(newdf$ind=="cases"),1:G]

        # new controls data frame
        newctrlsdf <- newdf[which(newdf$ind=="ctrls"),1:G]


        # new adj matrix for cases after the resample
        adjMat_cases <- adjacency(newcasesdf,
                                  selectCols = NULL,
                                  type = "signed",
                                  power = adjpower,
                                  corFnc = "cor", corOptions = list(use = "p"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  weightArgNames = c("weights.x", "weights.y"))

        # check new adj matrix for cases
        checkAdjMat(adjMat_cases, min = 0, max = 1)

        # new tom matrix for cases after the resample
        Tom_cases <- TOMsimilarity(adjMat_cases,
                                   TOMType = "signed",
                                   TOMDenom = "min",
                                   suppressTOMForZeroAdjacencies = FALSE,
                                   suppressNegativeTOM = FALSE,
                                   useInternalMatrixAlgebra = FALSE,
                                   verbose = 1,
                                   indent = 0)

        # new adj matrix for controls after the resample
        adjMat_ctrls <- adjacency(newctrlsdf,
                                  selectCols = NULL,
                                  type = "signed",
                                  power = adjpower,
                                  corFnc = "cor", corOptions = list(use = "p"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  weightArgNames = c("weights.x", "weights.y"))

        # check new adj matrix for controls
        checkAdjMat(adjMat_ctrls, min = 0, max = 1)

        # new tom matrix for controls after the resample
        Tom_ctrls <- TOMsimilarity(adjMat_ctrls,
                                   TOMType = "signed",
                                   TOMDenom = "min",
                                   suppressTOMForZeroAdjacencies = FALSE,
                                   suppressNegativeTOM = FALSE,
                                   useInternalMatrixAlgebra = FALSE,
                                   verbose = 1,
                                   indent = 0)

        # use lower triangle of permuted tom matrix for cases
        Perm_cases <- Tom_cases[lower.tri(Tom_cases)]

        # use lower triangle of permuted tom matrix for controls
        Perm_ctrls <- Tom_ctrls[lower.tri(Tom_ctrls)]

        # use lower triangle of permuted adj matrix for cases
        Perm_adj_cases <- adjMat_cases[lower.tri(adjMat_cases)]

        # use lower triangle of permuted adj matrix for controls
        Perm_adj_ctrls <- adjMat_ctrls[lower.tri(adjMat_ctrls)]



        # permuted wilconxin test results for tom matrix
        newresults <- wilcox.test(Perm_cases,Perm_ctrls, paired = TRUE)

        # store wilconxin test statistics
        StatPerm <- newresults$statistic

        # permuted t test results for tom matrix
        newtresults <- t.test(Perm_cases,Perm_ctrls, paired = TRUE)

        # store t test statistics
        ttestPerm <- newtresults$statistic

        # permuted wilconxin test results for adj matrix
        newadjresults <- wilcox.test(Perm_adj_cases,Perm_adj_ctrls,paired = TRUE)

        # store wilconxin test statistics
        AdjPerm <- newadjresults$statistic

        # permuted t test results for adj matrix
        newtresultsadj <- t.test(Perm_adj_cases,Perm_adj_ctrls,paired = TRUE)

        # store t test statistics
        ttestadjPerm <- newtresultsadj$statistic

        # permuted mad test for adj matrix
        PermmadADJ <- mean(abs(Perm_adj_cases-Perm_adj_ctrls))

        # permuted mad test for tom matrix
        PermmadTOM <- mean(abs(Perm_cases-Perm_ctrls))


        # permuted difference between cases and controls correlation vector
        Perm_cor_diff <- Perm_adj_cases - Perm_adj_ctrls

        # permuted dispersion index
        Perm_dis_ind <- sqrt(mean(Perm_cor_diff^2))

        Perm_cor_diff_tom <- Perm_cases - Perm_ctrls

        Perm_dis_ind_tom <- sqrt(mean(Perm_cor_diff_tom^2))

        Perm_cubtest_Adj <- mean(Perm_cor_diff^4)

        Perm_cubtest_Tom <- mean(Perm_cor_diff_tom^4)


        Perm_ranking_Adj <- rank(abs(Perm_cor_diff))

        Perm_ranking_Tom <- rank(abs(Perm_cor_diff_tom))


        Perm_mad_ranking_adj <- mean(abs(Perm_adj_cases-Perm_adj_ctrls)*Perm_ranking_Adj)
        Perm_mad_ranking_tom <- mean(abs(Perm_cases-Perm_ctrls)*Perm_ranking_Tom)

        Perm_dis_ranking_adj <- sqrt(mean((Perm_cor_diff^2)*Perm_ranking_Adj))
        Perm_dis_ranking_tom <- sqrt(mean((Perm_cor_diff_tom^2)*Perm_ranking_Tom))

        Perm_quadtest_ranking_adj <- mean((Perm_cor_diff^4)*Perm_ranking_Adj)
        Perm_quadtest_ranking_tom <- mean((Perm_cor_diff_tom^4)*Perm_ranking_Tom)

        # create a vector to store all above test statistics
        intermres <- c(StatPerm,ttestPerm,AdjPerm,ttestadjPerm,PermmadADJ,PermmadTOM,Perm_dis_ind,Perm_dis_ind_tom,Perm_cubtest_Adj,Perm_cubtest_Tom,
                       Perm_mad_ranking_adj,Perm_mad_ranking_tom,Perm_dis_ranking_adj,Perm_dis_ranking_tom,Perm_quadtest_ranking_adj,
                       Perm_quadtest_ranking_tom)

        # name that vector
        names(intermres) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffcubicAdj","DiffcubicTom",
                              "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

        return(intermres)
      }

      # store as a data frame
      pararesults <- as.data.frame(pararesults)

      T1storeperm [[k]] <- pararesults

      # calculate p value for wilcoxin tom
      pvaluePerm <- ((sum(pararesults$wTOM <= result$statistic))+1)/(permnum+1)
      store[k] <- pvaluePerm

      # calculate p value for t test tom
      pvaluePermttest <- ((sum(abs(pararesults$tTOM) >= abs(ttestresult$statistic)))+1)/(permnum+1)
      storettest[k] <- pvaluePermttest

      # calculate p values for wilcoxin adj
      pvaluePermAdj <- ((sum(pararesults$wADJ <= adjresult$statistic))+1)/(permnum+1)
      storeAdj[k] <- pvaluePermAdj

      # calculate p values for t test adj
      pvaluePermttestadj <- ((sum(abs(pararesults$tADJ) >= abs(adjttestresult$statistic)))+1)/(permnum+1)
      storettestadj[k] <- pvaluePermttestadj

      # calculate p values for mad adj
      pvalueMadADJ <- ((sum(pararesults$MadADJ >= madADJ))+1)/(permnum+1)
      storeMadADJ[k] <- pvalueMadADJ

      # calculate p values for mad tom
      pvalueMadTOM <- ((sum(pararesults$MadTOM >= madTOM))+1)/(permnum+1)
      storeMadTOM[k] <- pvalueMadTOM

      # calculate p values for Dispersion Index
      pvalueDisper <- ((sum(pararesults$DispersionIndexAdj >= obs_dis_ind))+1)/(permnum+1)
      storeDisperInd[k] <- pvalueDisper

      # calculate p values for Dispersion Index
      pvalueDisper_tom <- ((sum(pararesults$DispersionIndexTom >= obs_dis_ind_tom))+1)/(permnum+1)
      storeDisperInd_tom[k] <- pvalueDisper_tom

      pvaluequadAdj <- ((sum(pararesults$DiffcubicAdj >= obs_cubtest_Adj))+1)/(permnum+1)
      storecubAdj[k] <- pvaluequadAdj

      pvaluequadTom <- ((sum(pararesults$DiffcubicTom >= obs_cubtest_Tom))+1)/(permnum+1)
      storecubTom[k] <- pvaluequadTom

      ########### Newly added #####
      pvalue_mad_ranking_adj <- ((sum(pararesults$MadRankingAdj >= obs_mad_ranking_adj))+1)/(permnum+1)
      storeMadRankingAdj[k] <- pvalue_mad_ranking_adj

      pvalue_mad_ranking_tom <- ((sum(pararesults$MadRankingTom >= obs_mad_ranking_tom))+1)/(permnum+1)
      storeMadRankingTom[k] <- pvalue_mad_ranking_tom

      pvalue_dispersion_ranking_adj <- ((sum(pararesults$DispersionRankingAdj >= obs_dis_ranking_adj))+1)/(permnum+1)
      storeDisperRankingAdj[k] <- pvalue_dispersion_ranking_adj

      pvalue_dispersion_ranking_tom <- ((sum(pararesults$DispersionRankingTom >= obs_dis_ranking_tom))+1)/(permnum+1)
      storeDisperRankingTom[k] <- pvalue_dispersion_ranking_tom

      pvalue_quadtest_ranking_adj <- ((sum(pararesults$QuadtestRankingAdj >= obs_quadtest_Ranking_Adj))+1)/(permnum+1)
      storeQuadtestRankingAdj[k] <- pvalue_quadtest_ranking_adj

      pvalue_quadtest_ranking_tom <- ((sum(pararesults$QuandtestRankingTom >= obs_quadtest_Ranking_Tom))+1)/(permnum+1)
      storeQuadtestRankingTom[k] <- pvalue_quadtest_ranking_tom

      testRes <-rbind(pvaluePermAdj,pvaluePerm,pvaluePermttestadj,pvaluePermttest,pvalueMadADJ,pvalueMadTOM,pvalueDisper,pvalueDisper_tom,
                      pvaluequadAdj,pvaluequadTom,pvalue_mad_ranking_adj,pvalue_mad_ranking_tom,pvalue_dispersion_ranking_adj,pvalue_dispersion_ranking_tom,
                      pvalue_quadtest_ranking_adj,pvalue_quadtest_ranking_tom)

      storemodres[,k] <- rbind(storeAdj[k],store[k],storettestadj[k],storettest[k],storeMadADJ[k],storeMadTOM[k],storeDisperInd[k],storeDisperInd_tom[k],
                               storecubAdj[k],storecubTom[k],storeMadRankingAdj[k],storeMadRankingTom[k],storeDisperRankingAdj[k],storeDisperRankingTom[k],
                               storeQuadtestRankingAdj[k],storeQuadtestRankingTom[k])

      testRes <- as.data.frame(testRes)

      #   Resname <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","QuadtestAdj","QuadtestTom",
      #                 "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

      Resname <- c(paste("Wilcoxon_adj",adjpower,sep = "_"),"Wilcoxon_tom", paste("Ttest_adj",adjpower,sep = "_"),"Ttest_tom",paste("Mad_adj",adjpower,sep = "_"),"Mad_tom",
                   paste("DispersionIndex_adj",adjpower,sep = "_"),"DispersionIndex_tom",paste("QuarticTest_adj",adjpower,sep = "_"),"QuarticTest_tom",
                   paste("Mad_Ranking_adj",adjpower,sep = "_"),"Mad_Ranking_tom",paste("Dispersion_Ranking_adj",adjpower,sep = "_"),"Dispersion_Ranking_tom",
                   paste("QuarticTest_Ranking_adj",adjpower,sep = "_"),"Quandtest_Ranking_tom")

      rownames(testRes) <- Resname
      colnames(testRes) <- paste("Module",k,"stats",sep = "_")

      combres[[k]] <- testRes
    }
    if (dim(as.data.frame(Modcases[k]))[2] == 1)
    {
      print("code error session 3")
      warning(paste("Module",k,"only has one gene",sep = "_"))
      testRes <- NA
      combres[[k]] <- testRes
    }
  }

    finalRes <- do.call(cbind,combres)
  }
#    for (p in 1:dim(finalRes)) {
#      assign(colnames(i), )
#    }

#  return(finalRes)
  if(save_res == TRUE) {
  assign("Test_res",finalRes,envir = .GlobalEnv)
  assign("perm_res",T1storeperm,envir = .GlobalEnv)
  assign("obs_res",T1storeobs,envir = .GlobalEnv)
  assign("storep",storemodres,envir = .GlobalEnv)
  }
  invisible(return(list("Test_res"=finalRes,"perm_res"=T1storeperm,"obs_res"=T1storeobs,"storep"=storemodres)))

}
