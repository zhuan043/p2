#' Title
#'
#' @param rho Cov parameter
#' @param G   # of Genes
#' @param N   # of Subjects
#' @param M   # of Modules
#' @param permnum  # of Permutations
#' @param alpha    # Alpha level
#' @param propred  aa
#' @param redpt    bb
#' @param CVstru   cc
#' @param adjpower  dd
#' @param type      ee
#' @param save_res  gg
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
#'
#' @export
#'
#' @examples
#' trial1 <- DMCmod(rho=0.7, G=10, N=50, M=10, permnum=10, alpha=0.05, propred=c(0.1,0.3,0.5,0.7,0.9),redpt=c(0,0.3,0.5),CVstru="AR",adjpower=1,type = "signed", save_res = FALSE)


DMCmod <- function (rho, G, N, M=500, permnum=500, alpha=0.05, propred=c(0.1,0.3,0.5,0.7,0.9),redpt,CVstru,adjpower=1,type = "signed", save_res =TRUE) {
#  numCores <- detectCores()
 # numCores
 # registerDoParallel(numCores)  # use multicore, set to the number of our cores

 # chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

#  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
#    numCores <- 2L
#  } else {
    # use all cores in devtools::test()
#    numCores <- parallel::detectCores()
#  }

  require(WGCNA)

  numCores <- 2

  registerDoParallel(numCores)


  #K=3
  #rho=0.7
  #G=10
  #N=50
  #permnum = 20
  # prop0 = 0.1
  #alpha = 0.1
  #M = 5
  # increment =0.1
  #redpt = -1

  K <- length(propred)

  # store cases modules
  caselist <- list()

  # store controls modules
  controlist <- list()

  ## store test statistics

  # wilcoxin tom
  store <- rep(0,M)

  # wilcoxin adj
  storeAdj <- rep(0,M)

  # t-test tom
  storettest <- rep(0,M)

  # t-test adj
  storettestadj <- rep(0,M)

  # Mad adj
  storeMadADJ <- rep(0,M)

  # Mad tom
  storeMadTOM <- rep(0,M)

  # Dispersion Method
  storeDisperInd <- rep(0,M)

  storeDisperInd_tom <- rep(0,M)

  storecubAdj <- rep(0,M)

  storecubTom <- rep(0,M)
  ##### Newly Added
  storeMadRankingAdj <- rep(0,M)
  storeMadRankingTom <- rep(0,M)
  storeDisperRankingAdj <- rep(0,M)
  storeDisperRankingTom <- rep(0,M)
  storeQuadtestRankingAdj <- rep(0,M)
  storeQuadtestRankingTom <- rep(0,M)

  storeHD <- rep(0,M)
  storeCLX <- rep(0,M)
  storeScott <- rep(0,M)
  storeLC <- rep(0,M)

  Wilcox_count <- rep(0,K)
  Wilcox_adj_count <- rep(0,K)
  Ttest_count <- rep(0,K)
  Ttest_adj_count <- rep(0,K)
  MadADJ_count <- rep(0,K)
  MadTOM_count <- rep(0,K)
  Disper_count <- rep(0,K)
  Disper_tom_count <- rep(0,K)
  Quad_adj_count <- rep(0,K)
  Quad_tom_count <- rep(0,K)

  ##### Newly Added ranking methods#######
  Mad_ranking_adj_count <- rep(0,K)
  Mad_ranking_tom_count <- rep(0,K)
  Dispersion_ranking_adj_count <- rep(0,K)
  Dispersion_ranking_tom_count <- rep(0,K)
  Quadtest_ranking_adj_count <- rep(0,K)
  Quadtest_ranking_tom_count <- rep(0,K)

  ##### Test if two covriance matrix is different####
  HD_count <- rep(0,K)
  CLX_count  <- rep(0,K)
  Scott_count  <- rep(0,K)
  LC_count  <- rep(0,K)

  Wilcox_TOM_Power <- rep(0,K)
  Wilcox_ADJ_Power <- rep(0,K)
  Ttest_TOM_Power <- rep(0,K)
  Ttest_ADJ_Power <- rep(0,K)
  MadADJ_Power <- rep(0,K)
  MadTOM_Power <- rep(0,K)
  Disper_Power <- rep(0,K)
  Disper_Power_tom <- rep(0,K)
  Quad_adj_Power <- rep(0,K)
  Quad_tom_Power <- rep(0,K)

  ##### Newly Added ######
  Mad_ranking_adj_power <- rep(0,K)
  Mad_ranking_tom_power <- rep(0,K)
  Dispersion_ranking_adj_power <- rep(0,K)
  Dispersion_ranking_tom_power <- rep(0,K)
  Quadtest_ranking_adj_power <- rep(0,K)
  Quadtest_ranking_tom_power <- rep(0,K)

  HD_power <- rep(0,K)
  CLX_power  <- rep(0,K)
  Scott_power  <- rep(0,K)
  LC_power  <- rep(0,K)

  plot_store <- rep(0,K)

  methodNum <- 16

  interm_observed_res <- array(rep(0,M*methodNum),dim = c(M,methodNum))

  tdc_observed_res <- array(rep(0,M*methodNum*K),dim = c(M,methodNum,K))

  interm_perm_res <- array(rep(0,M*permnum*methodNum), dim = c(permnum,methodNum,M))

  permstat_res <- array(rep(0,M*permnum*methodNum*K), dim = c(permnum,methodNum,M,K))

#  propred <- c(0.1,0.3,0.5)


  for (k in 1:K) {

 # prop0 <- (0.1+(k-1)*increment)
 # plot_store[k] = prop0
  prop0 <- propred[k]
  plot_store[k] = prop0
  print(prop0)

  for (i in 1:M) {
    #i=1
    print(i)
    # rho=0.7
    #   G=10
    #  N=50
    # permnum = 20
    # prop0 = 0.3
    ## construct cases and controls modules

    # sigma matrix

    if (CVstru =="AR"||CVstru =="CS")
    {
    if (CVstru=="AR")
    {
      times <- 1:G
      #rho <- 0.5
      powerAr <- 1
      ###############
      H <- abs(outer(times, times, "-"))
      sigmaAr <- powerAr * rho^H
      sirow <- nrow(sigmaAr)
      sigmaAr[cbind(1:sirow, 1:sirow)] <- sigmaAr[cbind(1:sirow, 1:sirow)] * powerAr
      sigma <- sigmaAr
      sigmacases <- sigmaAr
    }

    if (CVstru=="CS")
    {
      sigma=matrix(rho,nrow=G,ncol=G)
      diag(sigma)=1
      sigmacases <- sigma
    }


    len <- length(sigma[lower.tri(sigma)])


#    redpt <- c(0,0.3,0.5)
#    prop0 <- 0.5
    dimredpt <- length(redpt)
    num <- rep(0,dimredpt)
    realnum <- rep(0,dimredpt)

    for (countredpt in 1:dimredpt) {

      num[countredpt] <- len*prop0/dimredpt

      if (num[countredpt] %% 2 ==0 ) {
        realnum[countredpt] <- num[countredpt]
      } else
      {realnum[countredpt] <- floor(num[countredpt])}
    }

    d <- sample(1:len, (sum(realnum)), replace = FALSE)

    apply_sigma <- matrix(rep(0,length(d)),ncol = (length(d)/length(realnum)))

    for (count_d in 1:length(realnum)) {
#      count_d <- 2
      apply_sigma[count_d,] <- d[((realnum[count_d]*(count_d-1))+1):(realnum[count_d]*count_d)]
      sigma[lower.tri(sigma)][apply_sigma[count_d,]] <- sigma[lower.tri(sigma)][apply_sigma[count_d]]*redpt[count_d]

    }


    # make sigma positively defined, it has to be positively defined for the rmvnorm function on run.
    sigma[upper.tri(sigma)] = t(sigma)[upper.tri(sigma)]
    sigmactrls <- make.positive.definite(sigma)
    # sigmactrls <- sigma
    #  is.positive.definite(sigmactrls)

    cases_1_module=rmvnorm(N,mean=rep(0,G),sigma=sigmactrls)
    caselist[[i]]<- cases_1_module

    ctrls_1_module=rmvnorm(N,mean=rep(0,G),sigma=sigmacases)
    controlist[[i]] <- ctrls_1_module

  }
   if (CVstru =="SM_CS")
   {
    splitnum <- 2

    initsp <- propred[k]
    sigma=matrix(rho,nrow=G,ncol=G)
    diag(sigma)=1
    sigmacases <- sigma

    ctrls_1_module=rmvnorm(N,mean=rep(0,G),sigma=sigmacases)

    #cases_test_module <- array(rep(0,N*(G/splitnum)*splitnum), dim = c(N,(G/splitnum),splitnum))



    #  interm_perm_res <- array(rep(0,M*permnum*10), dim = c(permnum,10,M))

    cases_test_module1 <- array(rep(0,N*G*initsp),dim = c(N,G*initsp))

    cases_test_module2 <- array(rep(0,N*G*(1-initsp)),dim = c(N,G*(1-initsp)))

    #  cases_test_module <- as.list(cases_test_module1,cases_test_module2)
    # gets <- abind(cases_test_module1,cases_test_module2,along = 3)


    siglen <- c(initsp,1-initsp)
    # rho <- c(0.1,0.7)

    for (counter0 in 1:splitnum)
    {
      #   sigma=matrix(rho,nrow=G,ncol=G)
      #  diag(sigma)=1
      #  sigmacases <- sigma

      # sigdim <- dim(sigma)[1]/splitnum


      #siglen <- (counter0-1)*initsp

      sigdim <- as.numeric(as.character(dim(sigma)[1]*siglen[counter0]))

      sigmasplit <- matrix(rho,nrow = sigdim,ncol = sigdim)
      diag(sigmasplit) <- 1

      sigmactrls <- make.positive.definite(sigmasplit)

      assign(paste("cases_test_module", counter0, sep=""), rmvnorm(N,mean=rep(0,sigdim),sigma=sigmactrls))


      # cases_test_module[,,counter0]=rmvnorm(N,mean=rep(0,sigdim),sigma=sigmactrls)


    }

    cases_1_module <- cbind(cases_test_module1,cases_test_module2)

   }

    if (CVstru =="SM_AR")
    {
      splitnum <- 2
      initsp <- propred[k]

      times <- 1:G
      #rho <- 0.5
      powerAr <- 1
      ###############
      H <- abs(outer(times, times, "-"))
      sigmaAr <- powerAr * rho^H
      sirow <- nrow(sigmaAr)
      sigmaAr[cbind(1:sirow, 1:sirow)] <- sigmaAr[cbind(1:sirow, 1:sirow)] * powerAr
      sigma <- sigmaAr
      sigmacases <- sigmaAr

      ctrls_1_module=rmvnorm(N,mean=rep(0,G),sigma=sigmacases)

      #cases_test_module <- array(rep(0,N*(G/splitnum)*splitnum), dim = c(N,(G/splitnum),splitnum))



      #  interm_perm_res <- array(rep(0,M*permnum*10), dim = c(permnum,10,M))

      cases_test_module1 <- array(rep(0,N*G*initsp),dim = c(N,G*initsp))

      cases_test_module2 <- array(rep(0,N*G*(1-initsp)),dim = c(N,G*(1-initsp)))

      #  cases_test_module <- as.list(cases_test_module1,cases_test_module2)
      # gets <- abind(cases_test_module1,cases_test_module2,along = 3)


      siglen <- c(initsp,1-initsp)
      # rho <- c(0.1,0.7)

      for (counter0 in 1:splitnum)
      {
        #   sigma=matrix(rho,nrow=G,ncol=G)
        #  diag(sigma)=1
        #  sigmacases <- sigma

        # sigdim <- dim(sigma)[1]/splitnum


        #siglen <- (counter0-1)*initsp

        sigdim <- as.numeric(as.character(dim(sigmacases)[1]*siglen[counter0]))

        #     sigmasplit <- matrix(rho,nrow = sigdim,ncol = sigdim)

        times <- 1:sigdim
        #rho <- 0.5
        powerAr <- 1
        ###############
        H <- abs(outer(times, times, "-"))
        sigmaAr <- powerAr * rho^H
        sirow <- nrow(sigmaAr)
        sigmaAr[cbind(1:sirow, 1:sirow)] <- sigmaAr[cbind(1:sirow, 1:sirow)] * powerAr
        sigma <- sigmaAr
        sigmasplit <- sigmaAr


        diag(sigmasplit) <- 1

        sigmactrls <- make.positive.definite(sigmasplit)

        assign(paste("cases_test_module", counter0, sep=""), rmvnorm(N,mean=rep(0,sigdim),sigma=sigmactrls))


        # cases_test_module[,,counter0]=rmvnorm(N,mean=rep(0,sigdim),sigma=sigmactrls)


      }

      cases_1_module <- cbind(cases_test_module1,cases_test_module2)
    }


    # calculate adj matrix for cases
    obs_adjMat_cases <- WGCNA::adjacency(cases_1_module,
                                  selectCols = NULL,
                                  type = type,
                                  power = adjpower,
                                  corFnc = "cor", corOptions = list(use = "p"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  weightArgNames = c("weights.x", "weights.y"))

    # check if adj matrix is valid
 #   checkAdjMat(obs_adjMat_cases, min = 0, max = 1)

    # calcualte adj matrix for controls
    obs_adjMat_ctrls <- WGCNA::adjacency(ctrls_1_module,
                                  selectCols = NULL,
                                  type = type,
                                  power = adjpower,
                                  corFnc = "cor", corOptions = list(use = "p"),
                                  weights = NULL,
                                  distFnc = "dist", distOptions = "method = 'euclidean'",
                                  weightArgNames = c("weights.x", "weights.y"))

    # check adj matrix for controls
 #   checkAdjMat(obs_adjMat_ctrls, min = 0, max = 1)

    obs_TOM_cases <- cor(cases_1_module,method = "pearson")
    obs_TOM_cases <- obs_TOM_cases[lower.tri(obs_TOM_cases)]

    obs_TOM_ctrls <- cor(ctrls_1_module, method = "pearson")
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

    covp <- testCov(cases_1_module,ctrls_1_module, method = "ALL", J = 2500, alpha = 0.05, n.core = 1)
    HDp <- covp$HD$p.value
    CLXp <- covp$CLX$p.value
    Scottp <- covp$Scott$p.value
    LCp <- covp[[4]]$p.value

    names(CLXp) <- NULL
    names(Scottp) <- NULL


    # register parallel computing with number of threads (for intel 4th gen and after, one core has two threads)
 #   registerDoParallel(numCores)



    ## parallel computing for permutations
    pararesults <- foreach (j = 1:permnum, .combine = 'rbind') %dopar% {

      require(WGCNA)

      # restore cases module as data frame
      casesPermDF <- as.data.frame(cases_1_module)

      # create labels for cases module
      casesPermDF$ind <- "cases"

      # restore controls module as data frame
      ctrlsPermDF <- as.data.frame(ctrls_1_module)

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
      adjMat_cases <- WGCNA::adjacency(newcasesdf,
                                selectCols = NULL,
                                type = type,
                                power = adjpower,
                                corFnc = "cor", corOptions = list(use = "p"),
                                weights = NULL,
                                distFnc = "dist", distOptions = "method = 'euclidean'",
                                weightArgNames = c("weights.x", "weights.y"))

      # check new adj matrix for cases
    #  checkAdjMat(adjMat_cases, min = 0, max = 1)


      # new adj matrix for controls after the resample
      adjMat_ctrls <- WGCNA::adjacency(newctrlsdf,
                                selectCols = NULL,
                                type = type,
                                power = adjpower,
                                corFnc = "cor", corOptions = list(use = "p"),
                                weights = NULL,
                                distFnc = "dist", distOptions = "method = 'euclidean'",
                                weightArgNames = c("weights.x", "weights.y"))

      # check new adj matrix for controls
 #     checkAdjMat(adjMat_ctrls, min = 0, max = 1)

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

      # create a vector to store all above test statistics
      # intermres <- c(StatPerm,ttestPerm,AdjPerm,ttestadjPerm,PermmadADJ,PermmadTOM,Perm_dis_ind,Perm_dis_ind_tom,Perm_cubtest_Adj,Perm_cubtest_Tom)
      # name that vector
      #  names(intermres) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffquadAdj","DiffquadTom")


      return(intermres)
    }

    interm_observed_res[i,] <- rbind(result$statistic,ttestresult$statistic,adjresult$statistic,
                                     adjttestresult$statistic,madADJ,madTOM,obs_dis_ind,obs_dis_ind_tom,
                                     obs_cubtest_Adj,obs_cubtest_Tom,obs_mad_ranking_adj,obs_mad_ranking_tom,obs_dis_ranking_adj,
                                     obs_dis_ranking_tom,obs_quadtest_Ranking_Adj,obs_quadtest_Ranking_Tom)

    interm_perm_res[,,i] <- pararesults
    # store as a data frame
    pararesults <- as.data.frame(pararesults)

    # calculate p value for wilcoxin tom
    pvaluePerm <- ((sum(pararesults$wTOM <= result$statistic))+1)/(permnum+1)
    store[i] <- pvaluePerm

    # calculate p value for t test tom
    pvaluePermttest <- ((sum(abs(pararesults$tTOM) >= abs(ttestresult$statistic)))+1)/(permnum+1)
    storettest[i] <- pvaluePermttest

    # calculate p values for wilcoxin adj
    pvaluePermAdj <- ((sum(pararesults$wADJ <= adjresult$statistic))+1)/(permnum+1)
    storeAdj[i] <- pvaluePermAdj

    # calculate p values for t test adj
    pvaluePermttestadj <- ((sum(abs(pararesults$tADJ) >= abs(adjttestresult$statistic)))+1)/(permnum+1)
    storettestadj[i] <- pvaluePermttestadj

    # calculate p values for mad adj
    pvalueMadADJ <- ((sum(pararesults$MadADJ >= madADJ))+1)/(permnum+1)
    storeMadADJ[i] <- pvalueMadADJ

    # calculate p values for mad tom
    pvalueMadTOM <- ((sum(pararesults$MadTOM >= madTOM))+1)/(permnum+1)
    storeMadTOM[i] <- pvalueMadTOM

    # calculate p values for Dispersion Index
    pvalueDisper <- ((sum(pararesults$DispersionIndexAdj >= obs_dis_ind))+1)/(permnum+1)
    storeDisperInd[i] <- pvalueDisper

    # calculate p values for Dispersion Index
    pvalueDisper_tom <- ((sum(pararesults$DispersionIndexTom >= obs_dis_ind_tom))+1)/(permnum+1)
    storeDisperInd_tom[i] <- pvalueDisper_tom

    pvaluequadAdj <- ((sum(pararesults$DiffcubicAdj >= obs_cubtest_Adj))+1)/(permnum+1)
    storecubAdj[i] <- pvaluequadAdj

    pvaluequadTom <- ((sum(pararesults$DiffcubicTom >= obs_cubtest_Tom))+1)/(permnum+1)
    storecubTom[i] <- pvaluequadTom

    ########### Newly added #####
    pvalue_mad_ranking_adj <- ((sum(pararesults$MadRankingAdj >= obs_mad_ranking_adj))+1)/(permnum+1)
    storeMadRankingAdj[i] <- pvalue_mad_ranking_adj

    pvalue_mad_ranking_tom <- ((sum(pararesults$MadRankingTom >= obs_mad_ranking_tom))+1)/(permnum+1)
    storeMadRankingTom[i] <- pvalue_mad_ranking_tom

    pvalue_dispersion_ranking_adj <- ((sum(pararesults$DispersionRankingAdj >= obs_dis_ranking_adj))+1)/(permnum+1)
    storeDisperRankingAdj[i] <- pvalue_dispersion_ranking_adj

    pvalue_dispersion_ranking_tom <- ((sum(pararesults$DispersionRankingTom >= obs_dis_ranking_tom))+1)/(permnum+1)
    storeDisperRankingTom[i] <- pvalue_dispersion_ranking_tom

    pvalue_quadtest_ranking_adj <- ((sum(pararesults$QuadtestRankingAdj >= obs_quadtest_Ranking_Adj))+1)/(permnum+1)
    storeQuadtestRankingAdj[i] <- pvalue_quadtest_ranking_adj

    pvalue_quadtest_ranking_tom <- ((sum(pararesults$QuandtestRankingTom >= obs_quadtest_Ranking_Tom))+1)/(permnum+1)
    storeQuadtestRankingTom[i] <- pvalue_quadtest_ranking_tom

    storeHD[i] <- HDp
    storeCLX[i] <- CLXp
    storeScott[i] <- Scottp
    storeLC[i] <- LCp

  }

  tdc_observed_res[,,k] <- interm_observed_res

  permstat_res[,,,k] <- interm_perm_res

  # count # of significant

  # add <= to alpha

  Wilcox_count[k] <- sum(store < alpha)
  Wilcox_adj_count[k] <- sum(storeAdj < alpha)
  Ttest_count[k] <- sum(storettest < alpha)
  Ttest_adj_count[k] <- sum(storettestadj < alpha)
  MadADJ_count[k] <- sum(storeMadADJ < alpha)
  MadTOM_count[k] <- sum(storeMadTOM < alpha)
  Disper_count[k] <- sum(storeDisperInd < alpha)
  Disper_tom_count[k] <- sum(storeDisperInd_tom < alpha)
  Quad_adj_count[k] <- sum(storecubAdj < alpha)
  Quad_tom_count[k] <- sum(storecubTom < alpha)
  ##### Newly Added #######
  Mad_ranking_adj_count[k] <- sum(storeMadRankingAdj < alpha)
  Mad_ranking_tom_count[k] <- sum(storeMadRankingTom < alpha)
  Dispersion_ranking_adj_count[k] <- sum(storeDisperRankingAdj < alpha)
  Dispersion_ranking_tom_count[k] <- sum(storeDisperRankingTom < alpha)
  Quadtest_ranking_adj_count[k] <- sum(storeQuadtestRankingAdj < alpha)
  Quadtest_ranking_tom_count[k] <- sum(storeQuadtestRankingTom < alpha)

  HD_count[k] <- sum(storeHD < alpha)
  CLX_count[k]  <- sum(storeCLX < alpha)
  Scott_count[k]  <- sum(storeScott < alpha)
  LC_count[k]  <- sum(storeLC < alpha)

  # calculate Power
  Wilcox_TOM_Power <- Wilcox_count/M
  Wilcox_ADJ_Power <- Wilcox_adj_count/M
  Ttest_TOM_Power <- Ttest_count/M
  Ttest_ADJ_Power <- Ttest_adj_count/M
  MadADJ_Power <- MadADJ_count/M
  MadTOM_Power <- MadTOM_count/M
  Disper_Power <- Disper_count/M
  Disper_tom_Power <- Disper_tom_count/M
  Quad_adj_Power <- Quad_adj_count/M
  Quad_tom_Power <- Quad_tom_count/M

  ##### Newly Added ######
  Mad_ranking_adj_power <- Mad_ranking_adj_count/M
  Mad_ranking_tom_power <- Mad_ranking_tom_count/M
  Dispersion_ranking_adj_power <- Dispersion_ranking_adj_count/M
  Dispersion_ranking_tom_power <- Dispersion_ranking_tom_count/M
  Quadtest_ranking_adj_power <- Quadtest_ranking_adj_count/M
  Quadtest_ranking_tom_power <- Quadtest_ranking_tom_count/M

  HD_power <- HD_count/M
  CLX_power  <- CLX_count/M
  Scott_power  <- Scott_count/M
  LC_power  <- LC_count/M

  }

  plot_Wilcox_TOM_Power <- as.data.frame(Wilcox_TOM_Power)
  plot_Wilcox_ADJ_Power <- as.data.frame(Wilcox_ADJ_Power)
  plot_Ttest_TOM_Power <- as.data.frame(Ttest_TOM_Power)
  plot_Ttest_ADJ_Power <- as.data.frame(Ttest_ADJ_Power)
  plot_MadADJ_Power <- as.data.frame(MadADJ_Power)
  plot_MadTOM_Power <- as.data.frame(MadTOM_Power)
  plot_Disper_Power <- as.data.frame(Disper_Power)
  plot_Disper_tom_Power <- as.data.frame(Disper_tom_Power)
  plot_Quad_adj_Power <- as.data.frame(Quad_adj_Power)
  plot_Quad_tom_Power <- as.data.frame(Quad_tom_Power)
  ##### Newly Added ######
  plot_Mad_ranking_adj_power <- as.data.frame(Mad_ranking_adj_power)
  plot_Mad_ranking_tom_power <- as.data.frame(Mad_ranking_tom_power)
  plot_Dispersion_ranking_adj_power <- as.data.frame(Dispersion_ranking_adj_power)
  plot_Dispersion_ranking_tom_power <- as.data.frame(Dispersion_ranking_tom_power)
  plot_Quadtest_ranking_adj_power <- as.data.frame(Quadtest_ranking_adj_power)
  plot_Quadtest_ranking_tom_power <- as.data.frame(Quadtest_ranking_tom_power)

  plot_HD_power <- as.data.frame(HD_power)
  plot_CLX_power <- as.data.frame(CLX_power)
  plot_Scott_power <- as.data.frame(Scott_power)
  plot_LC_power <- as.data.frame(LC_power)

  plot_Wilcox_TOM_Power$dropOut <- plot_store
  plot_Wilcox_ADJ_Power$dropOut <- plot_store
  plot_Ttest_TOM_Power$dropOut <- plot_store
  plot_Ttest_ADJ_Power$dropOut <- plot_store
  plot_MadADJ_Power$dropOut <- plot_store
  plot_MadTOM_Power$dropOut <- plot_store
  plot_Disper_Power$dropOut <- plot_store
  plot_Disper_tom_Power$dropOut <- plot_store
  plot_Quad_adj_Power$dropOut <- plot_store
  plot_Quad_tom_Power$dropOut <- plot_store
  ####### Newly Added ######
  plot_Mad_ranking_adj_power$dropOut  <- plot_store
  plot_Mad_ranking_tom_power$dropOut  <- plot_store
  plot_Dispersion_ranking_adj_power$dropOut  <- plot_store
  plot_Dispersion_ranking_tom_power$dropOut  <- plot_store
  plot_Quadtest_ranking_adj_power$dropOut  <- plot_store
  plot_Quadtest_ranking_tom_power$dropOut  <- plot_store

  plot_HD_power$dropOut <- plot_store
  plot_CLX_power$dropOut <- plot_store
  plot_Scott_power$dropOut <- plot_store
  plot_LC_power$dropOut <- plot_store

  plot_Wilcox_TOM_Power$method <- "Wilcox_TOM"
  plot_Wilcox_ADJ_Power$method <- "Wilcox_ADJ"
  plot_Ttest_TOM_Power$method <- "Ttest_TOM"
  plot_Ttest_ADJ_Power$method <- "Ttest_ADJ"
  plot_MadADJ_Power$method <- "MadADJ_Power"
  plot_MadTOM_Power$method <- "MadTOM_Power"
  plot_Disper_Power$method <- "Disper_ADJ_Power"
  plot_Disper_tom_Power$method <- "Disper_TOM_Power"
  plot_Quad_adj_Power$method <- "Quad_ADJ_Power"
  plot_Quad_tom_Power$method <- "Quad_TOM_Power"
  #### Newly Added #####
  plot_Mad_ranking_adj_power$method <- "Mad_ranking_ADJ_power"
  plot_Mad_ranking_tom_power$method <- "Mad_ranking_TOM_power"
  plot_Dispersion_ranking_adj_power$method <- "Dispersion_ranking_ADJ_power"
  plot_Dispersion_ranking_tom_power$method <- "Dispersion_ranking_TOM_power"
  plot_Quadtest_ranking_adj_power$method <- "Quadtest_ranking_ADJ_power"
  plot_Quadtest_ranking_tom_power$method <- "Quadtest_ranking_TOM_power"

  plot_HD_power$method  <- "HD_power"
  plot_CLX_power$method  <- "CLX_power"
  plot_Scott_power$method  <- "Scott_power"
  plot_LC_power$method  <- "LC_power"

  colnames(plot_Wilcox_TOM_Power)[1] <- "Power"
  colnames(plot_Wilcox_ADJ_Power)[1] <- "Power"
  colnames(plot_Ttest_TOM_Power)[1] <- "Power"
  colnames(plot_Ttest_ADJ_Power)[1] <- "Power"
  colnames(plot_MadADJ_Power)[1] <- "Power"
  colnames(plot_MadTOM_Power)[1] <- "Power"
  colnames(plot_Disper_Power)[1] <- "Power"
  colnames(plot_Disper_tom_Power)[1] <- "Power"
  colnames(plot_Quad_adj_Power)[1] <- "Power"
  colnames(plot_Quad_tom_Power)[1] <- "Power"
  #### Newly Added #####
  colnames(plot_Mad_ranking_adj_power)[1] <- "Power"
  colnames(plot_Mad_ranking_tom_power)[1] <- "Power"
  colnames(plot_Dispersion_ranking_adj_power)[1] <- "Power"
  colnames(plot_Dispersion_ranking_tom_power)[1] <- "Power"
  colnames(plot_Quadtest_ranking_adj_power)[1] <- "Power"
  colnames(plot_Quadtest_ranking_tom_power)[1] <- "Power"

  colnames(plot_HD_power)[1] <- "Power"
  colnames(plot_CLX_power)[1] <- "Power"
  colnames(plot_Scott_power)[1] <- "Power"
  colnames(plot_LC_power)[1] <- "Power"

  plotdf <- rbind(plot_Wilcox_TOM_Power,plot_Wilcox_ADJ_Power,plot_Ttest_TOM_Power,
                  plot_Ttest_ADJ_Power,plot_MadADJ_Power,plot_MadTOM_Power,plot_Disper_Power,
                  plot_Disper_tom_Power,plot_Quad_adj_Power,plot_Quad_tom_Power,
                  plot_Mad_ranking_adj_power,plot_Mad_ranking_tom_power,plot_Dispersion_ranking_adj_power,
                  plot_Dispersion_ranking_tom_power,plot_Quadtest_ranking_adj_power,plot_Quadtest_ranking_tom_power,
                  plot_HD_power,plot_CLX_power,plot_Scott_power,plot_LC_power)

  plotdf$linety <- "solid"
  plotdf$linety[plotdf$method == "Wilcox_TOM"] <- "dashed"
  plotdf$linety[plotdf$method == "Ttest_TOM"] <- "dashed"
  plotdf$linety[plotdf$method == "MadTOM_Power"] <- "dashed"
  plotdf$linety[plotdf$method == "Disper_TOM_Power"] <- "dashed"
  plotdf$linety[plotdf$method == "Quad_TOM_Power"] <- "dashed"
  ### Newly Added
  plotdf$linety[plotdf$method == "Mad_ranking_TOM_power"] <- "dashed"
  plotdf$linety[plotdf$method == "Dispersion_ranking_TOM_power"] <- "dashed"
  plotdf$linety[plotdf$method == "Quadtest_ranking_TOM_power"] <- "dashed"

  plotdf$linety[plotdf$method == "HD_power"] <- "dotted"
  plotdf$linety[plotdf$method == "CLX_power"] <- "dotted"
  plotdf$linety[plotdf$method == "Scott_power"] <- "dotted"
  plotdf$linety[plotdf$method == "LC_power"] <- "dotted"

 # powerPlot <- ggplot(plotdf,aes(y = Power,x = dropOut,color = method)) +
  #  geom_line(aes(linetype=linety)) +
  #  ggtitle("Power plot for comparing different methods")

  xvalue <- plotdf$dropOut

  plotdf$method <- factor(plotdf$method, levels = c("Disper_ADJ_Power","Disper_TOM_Power",
                                                    "MadADJ_Power","MadTOM_Power","Quad_ADJ_Power","Quad_TOM_Power","Ttest_ADJ","Ttest_TOM","Wilcox_ADJ","Wilcox_TOM",
                                                    "Mad_ranking_ADJ_power","Mad_ranking_TOM_power","Dispersion_ranking_ADJ_power","Dispersion_ranking_TOM_power",
                                                    "Quadtest_ranking_ADJ_power","Quadtest_ranking_TOM_power","HD_power","CLX_power","Scott_power",
                                                    "LC_power"))

  #MyDF$MyCat<-factor(MyDF$MyCat, levels=c("Dark","DarkLight","Medium","LightDark","Light"))


  powerPlot <- ggplot(plotdf,aes(y = Power,x = dropOut,color = method,linetype=method)) +
    geom_line() +
    scale_x_discrete(name="Prop. of changed correlations",limits= xvalue,expand = expand_scale(mult = 0.05, add = 0))+
    scale_y_continuous(limits = c(0,1))+
    #geom_point() +
    ggtitle("Power plot for comparing DMC tests") +
    scale_color_manual(name="Method",labels= c("Disper_ADJ_Power","Disper_TOM_Power",
                                               "MadADJ_Power","MadTOM_Power","Quad_ADJ_Power","Quad_TOM_Power","Ttest_ADJ","Ttest_TOM","Wilcox_ADJ","Wilcox_TOM",
                                               "Mad_ranking_ADJ_power","Mad_ranking_TOM_power","Dispersion_ranking_ADJ_power","Dispersion_ranking_TOM_power",
                                               "Quadtest_ranking_ADJ_power","Quadtest_ranking_TOM_power","HD_power","CLX_power","Scott_power",
                                               "LC_power"),
                       values = c("blue","blue", "red","red", "green","green", "orange","orange","black","black","yellow","yellow","pink","pink","purple","purple","burlywood","peachpuff4","grey50","tan2")) +
    scale_linetype_manual(name="Method", labels= c("Disper_ADJ_Power","Disper_TOM_Power",
                                                   "MadADJ_Power","MadTOM_Power","Quad_ADJ_Power","Quad_TOM_Power","Ttest_ADJ","Ttest_TOM","Wilcox_ADJ","Wilcox_TOM",
                                                   "Mad_ranking_ADJ_power","Mad_ranking_TOM_power","Dispersion_ranking_ADJ_power","Dispersion_ranking_TOM_power",
                                                   "Quadtest_ranking_ADJ_power","Quadtest_ranking_TOM_power","HD_power","CLX_power","Scott_power",
                                                   "LC_power"),
                          values = c(rep(c("solid", "dashed"),(methodNum/2)),rep(c("dotted"),4)))

  colnames(interm_observed_res) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffquadAdj","DiffquadTom",
                                     "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

  colnames(tdc_observed_res) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffquadAdj","DiffquadTom",
                                  "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

  colnames(interm_perm_res) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffquadAdj","DiffquadTom",
                                 "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")

  colnames(permstat_res) <- c("wTOM","tTOM","wADJ","tADJ","MadADJ","MadTOM","DispersionIndexAdj","DispersionIndexTom","DiffquadAdj","DiffquadTom",
                              "MadRankingAdj","MadRankingTom","DispersionRankingAdj","DispersionRankingTom","QuadtestRankingAdj","QuandtestRankingTom")


  modtype <- paste("DMC_mod",CVstru,length(redpt),"types of changes",sep = "_")


  modname <- paste(modtype,", rho =",rho,", G =",G,", N =",N,", M =",M,", # of perms = ",permnum,", alpha = ",alpha,", type = ",type,", adj_power = ",adjpower)

  storeintolist <- list(powerPlot,plotdf,modname)

  if (save_res == TRUE) {

  save(interm_observed_res,file = paste(modtype,rho,G,N,M,permnum,alpha,"interm_observed_res.RData",sep = "_"))
  save(tdc_observed_res, file = paste(modtype,rho,G,N,M,permnum,alpha,"tdc_observed_res.RData",sep = "_"))

  save(interm_perm_res, file = paste(modtype,rho,G,N,M,permnum,alpha,"interm_perm_res.RData",sep = "_"))
  save(permstat_res,file = paste(modtype,rho,G,N,M,permnum,alpha,"permstat_res.RData",sep = "_"))

  save(storeintolist,file = paste(modtype,rho,G,N,M,permnum,alpha,"tdc_final_result.RData",sep = "_"))

  }

  assign("plotAndDF",storeintolist,envir = .GlobalEnv)

  invisible(return(list("PowerPlot"=powerPlot,"Plot_dataFrame"=plotdf,"scenario"=modname)))

}
