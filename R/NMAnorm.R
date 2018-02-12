#' Bayesian network meta-analysis
#'
#' This function links in the NMAmodel.f which builds
#' the network meta-regression proposed by Li et al. (2017)
#'
#' @param Y vector of length N for the response variable (aggregate mean)
#' @param SD vector of length N for the standard deviation
#' @param Npt vector of length N for the number of patients
#' @param X design matrix, N by p
#' @param ID_study study ID
#' @param ID_arm arm ID
#' @param ID_model model ID
#' @param Nrep MCMC sample size
#' @param Nthinning thinning size
#' @param Nwarmup warm-up size
#' @param Iseed random seed for the FORTRAN random number generator
#' @param alpha desired HPD credible level
#' @param DIC (Output) a vector of length 5
#' @param logCPO (Output) a vector of length K, K is the number of trials
#' @param betapo (Output) a matrix of nx by 4, nx is the number of covariates
#' @param gampo (Output) a matrix of nT by 4, nT is the number of treatments
#' @param tau2po (Output) a matrix of nT by 4
#' @param eRho (Output) a matrix of nT by nT
#' @param seRho (Output) a matrix of nT by nT 
#' @param Rholow (Output) a matrix of nT by nT
#' @param Rhoupp (Output) a matrix of nT by nT
#' @return Bayesian comparison criteria (DIC and LPML), posterior
#' estimate of coefficients of covariates, treatment effects, variances
#' and correlation matrix of treatment random effects
#' @useDynLib NMAnalysis, .registration = TRUE
#' @export
NMAmodel <- function(ID_model, Nrep=100L, Nthinning=1L, Nwarmup=10L,
                     Iseed=1234567L, alpha=0.05,
                     ID_study=NMAdata[, 1], ID_arm=NMAdata[, 2], 
                     Npt=NMAdata[, 3], Y=NMAdata[, 4], SD=NMAdata[, 5], 
                     X=matrix(as.double(as.matrix(NMAdata[, 6:15])), nrow = 73, ncol = 10), 
                     DIC=rep(0.0,5),
                     logCPO=rep(0.0,29), betapo=matrix(0.0, 10, 4), 
                     gampo=matrix(0.0, 11, 4), 
                     tau2po=matrix(0.0, 11, 4), eRho=matrix(0.0, 11, 11), 
                     seRho=matrix(0.0, 11, 11), Rholow=matrix(0.0, 11, 11), 
                     Rhoupp=matrix(0.0, 11, 11)){
  if(!is.null(Iseed)){set.seed(Iseed)} else return('require a seed value')
  if(ID_model == 1){
    ID_group <- rep(1L, 11)
  } else if(ID_model == 2){
    ID_group <- c(1L,2L,2L,2L,2L,2L,3L,4L,4L,4L,4L)
  } else if(ID_model == 3){
    ID_group <- c(1L,2L,3L,3L,3L,3L,4L,5L,5L,5L,5L)
  } else if(ID_model == 4){
    ID_group <- c(1L,3L,3L,3L,2L,3L,4L,5L,5L,5L,5L)
  } else if(ID_model == 5){
    ID_group <- c(1L,3L,2L,3L,3L,3L,4L,5L,5L,5L,5L)
  } else if(ID_model == 6){
    ID_group <- c(1L,2L,4L,4L,3L,4L,5L,6L,6L,6L,6L)
  } else if(ID_model == 7){
    ID_group <- c(1L,2L,3L,4L,4L,4L,5L,6L,6L,6L,6L)
  } else if(ID_model == 8){
    ID_group <- c(1L,4L,2L,4L,3L,4L,5L,6L,6L,6L,6L)
  } 
  ret <- .Fortran("NMAmodel", y1 = as.double(Y), sd1 = as.double(SD),
                  npt1 = as.integer(Npt),
                  x1 = X, ids1 = as.integer(ID_study), iarm1 = as.integer(ID_arm),
                  igroup1 = ID_group, nrep = as.integer(Nrep),
                  nthin = as.integer(Nthinning),
                  nwarm = as.integer(Nwarmup), iseed = as.integer(Iseed),
                  conf = as.double(alpha),
                  DIC = as.double(DIC), logCPO = as.double(logCPO),
                  Beta = betapo, Gam = gampo, tau2 = tau2po, 
                  eRho1 = eRho, seRho1 = seRho, Rholow1 = Rholow,
                  Rhoupp1 = Rhoupp, PACKAGE = "NMAnalysis")
  dic <- ret$DIC[1]
  pd <- ret$DIC[2]
  bardic <- ret$DIC[3]
  dicbar <- ret$DIC[4]
  lpml <- ret$DIC[5]
  ret$Beta <- as.data.frame(ret$Beta)
  names(ret$Beta) <- c('MEAN', 'SD', 'HPD Lower', 'HPD Upper')
  ret$Gam <- as.data.frame(ret$Gam)
  names(ret$Beta) <- c('MEAN', 'SD', 'HPD Lower', 'HPD Upper')
  ret$tau2 <- as.data.frame(ret$tau2)
  names(ret$Beta) <- c('MEAN', 'SD', 'HPD Lower', 'HPD Upper')
  output <- list(Size_of_Simmulation = Nrep, 
                Size_of_Thinning = Nthinning, 
                Size_of_Warmup = Nwarmup,
                ID_group = ID_group,
                DIC = dic, pD = pd, barDIC = bardic, DICbar = dicbar,
                LPML = lpml, Beta_Posterior = ret$Beta, 
                Gamma_Posterior = ret$Gam, Tau2_Posterior = ret$tau2,
                Rho_MEAN = ret$eRho1, Rho_SD= ret$seRho1, 
                Rho_HPD_Low = ret$Rholow1,
                Rho_HPD_Upp = ret$Rhoupp1)
  return(output)
}





