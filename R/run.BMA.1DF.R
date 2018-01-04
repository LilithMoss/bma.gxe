#' BMA function for 1-degree-of-freedom Test
#'
#' @param Y input A vector of binary outcomes/phenotypes. Defaults to NULL.
#' @param E input A vector of binary exposure status. Defaults to NULL.
#' @param G input A vector of binary genotype. Defaults to NULL.
#' @param Cov of all covariate variables. Defaults to NULL.
#' @param cc input A number. Defaults to 0.5.
#' @param co input A number. Defaults to 0.5.
#' @param phi input A number. Defaults to 1.
#' @param psi input A number. Defaults to 1000.
#' @param formul input A string in the form: "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+... . Defaults to NULL.
#' @return Vector containing BMA interaction effect estimate, standard error, Z-score, and p-value.
#' @import BMA
#' @export

# library(BMA)
run.BMA.1DF <- function(Y=NULL,E=NULL,G=NULL,Cov=NULL,cc=0.05,co=0.05,phi=1,psi=1000,formul=NULL) {
  if(length(Cov)>0){
    cterms <- ncol(Cov)
    covnam <- paste0("C",seq(1,ncol(Cov)))
    names(Cov) <- covnam
    nterms <- 2*cterms
    Sample.complete <- cbind(Y,G,E,Cov)
    Sample.complete[Sample.complete != 0 & Sample.complete != 1] <- NA
    Sample.complete <- Sample.complete[complete.cases(Sample.complete), ]

    if(length(unique(Sample.complete[,1]))>1 & length(unique(Sample.complete[,2]))>1 & length(unique(Sample.complete[,3]))>1){
      Dat <- as.data.frame(ftable(Sample.complete))
      colnames(Dat)[ncol(Dat)] <- "Count"
      pmw <- c(cc,co)
      pmw <- pmw/sum(pmw)
      alt.models<- c(1:2)

      if(is.null(formul)){
        models <- rbind(c(1,1,1,1,1,1,1,rep(1,nterms)), c(1,1,1,0,1,1,1,rep(1,nterms)))
        covariate.terms <- c(paste0("C", 1:cterms),paste0("Y:",covnam))
        base.terms <- c("E","G","E:G","Y","E:Y","G:Y","E:G:Y")
        model.terms <- c(base.terms,covariate.terms)
        cov.model.terms <- paste0("Count~",paste(model.terms,collapse="+"))
        X <- as.data.frame(model.matrix(as.formula(cov.model.terms),data=Dat))[,-1]
        main.terms <- c("E1","G1","Y1","E1:G1","E1:Y1","G1:Y1","E1:G1:Y1")
        extra.terms <- names(X[ , -which(names(X) %in% main.terms)])
        Xr <- X[,c(main.terms,extra.terms)]
      } else {
        cov.model.terms <- formul
        X <- as.data.frame(model.matrix(as.formula(cov.model.terms),data=Dat))[,-1]
        main.terms <- c("E1","G1","Y1","E1:G1","E1:Y1","G1:Y1","E1:G1:Y1")
        extra.terms <- names(X[ , -which(names(X) %in% main.terms)])
        Xr <- X[,c(main.terms,extra.terms)]
        models <- rbind(c(1,1,1,1,1,1,1,rep(1,(ncol(Xr)-7))), c(1,1,1,0,1,1,1,rep(1,(ncol(Xr)-7))))
      }
      p=ncol(Xr)
      r <- BMA::glib(Xr,y=Dat$Count, error="poisson", link = "log", phi=c(1), psi=c(1000),models=models, pmw=pmw, priormean=rep(0, p+1),output.postvar=T)
      Int.est <- r$posterior$mean[match("E1:G1:Y1",names(Xr))]
      Int.sd <- r$posterior$sd[match("E1:G1:Y1",names(Xr))]
      PrM.data <- sum( r$bf$postprob[alt.models])
      PrM <- sum( pmw[alt.models])
      Z.score <- Int.est/Int.sd
      p.value <- 2*pnorm(-abs(Z.score))
      bma.result <- c(Int.est, Int.sd, Z.score, p.value)
      return(bma.result)
    } else {
      bma.result <- c(NA,NA,NA,NA)
      return(bma.result)
    }

  } else {
    if(length(unique(Sample.complete[,1]))>1 & length(unique(Sample.complete[,2]))>1 & length(unique(Sample.complete[,3]))>1){
      Sample.complete <- as.data.frame(cbind(Y,G,E))
      Sample.complete[Sample.complete != 0 & Sample.complete != 1] <- NA
      Sample.complete <- Sample.complete[complete.cases(Sample.complete), ]
      Dat <- as.data.frame(ftable(Sample.complete))
      colnames(Dat)[ncol(Dat)] <- "Count"
      pmw <- c(cc,co)
      pmw <- pmw/sum(pmw)
      alt.models<- c(1:2)
      if(is.null(formul)){
        models <- rbind(c(1,1,1,1,1,1,1), c(1,1,1,0,1,1,1))
        base.terms <- c("E","G","E:G","Y","E:Y","G:Y","E:G:Y")
        model.terms <- c(base.terms)
        cov.model.terms <- paste0("Count~",paste(model.terms,collapse="+"))
        X <- as.data.frame(model.matrix(as.formula(cov.model.terms),data=Dat))[,-1]
        main.terms <- c("E1","G1","Y1","E1:G1","E1:Y1","G1:Y1","E1:G1:Y1")
        extra.terms <- names(X[ , -which(names(X) %in% main.terms)])
        Xr <- X[,c(main.terms,extra.terms)]
      } else {
        cov.model.terms <- formul
        X <- as.data.frame(model.matrix(as.formula(cov.model.terms),data=Dat))[,-1]
        main.terms <- c("E1","G1","Y1","E1:G1","E1:Y1","G1:Y1","E1:G1:Y1")
        extra.terms <- names(X[ , -which(names(X) %in% main.terms)])
        Xr <- X[,c(main.terms,extra.terms)]
        models <- rbind(c(1,1,1,1,1,1,1,rep(1,(ncol(Xr)-7))), c(1,1,1,0,1,1,1,rep(1,(ncol(Xr)-7))))
      }
      p=ncol(Xr)
      r <- BMA::glib(Xr,y=Dat$Count, error="poisson", link = "log", phi=c(1), psi=c(1000),models=models, pmw=pmw, priormean=rep(0, p+1),output.postvar=T)
      Int.est <- r$posterior$mean[match("E1:G1:Y1",names(Xr))]
      Int.sd <- r$posterior$sd[match("E1:G1:Y1",names(Xr))]
      PrM.data <- sum( r$bf$postprob[alt.models])
      PrM <- sum( pmw[alt.models])
      Z.score <- Int.est/Int.sd
      p.value <- 2*pnorm(-abs(Z.score))
      bma.result <- c(Int.est, Int.sd, Z.score, p.value)
      return(bma.result)
    } else {
      return(c(NA,NA,NA,NA))
    }
  }
}
