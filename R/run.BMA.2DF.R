#' BMA function for 2-degree-of-freedom Test
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
#' @return The p-value for a BMA joint test of marginal and interaction effect.
#' @export

# library(BMA)
run.BMA.2DF <- function(Y=NULL,E=NULL,G=NULL,Cov=NULL,cc=0.05,co=0.05,phi=1,psi=1000,formul=NULL) {
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
        X <- as.data.frame(model.matrix(as.formula(cov.model.terms),data=Dat))[,2:14]
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
      Ebeta <- matrix(c(r$posterior$mean[6],r$posterior$mean[7]),ncol=1)
      Gamma1 <- matrix(r$posterior.bymodel$var[[1]],ncol=ncol(r$posterior.bymodel$var[[1]]))
      Gamma2 <- matrix(r$posterior.bymodel$var[[2]],ncol=ncol(r$posterior.bymodel$var[[2]]))
      Gamma1.1 <- Gamma1[7:8,7:8]
      Gamma2.2 <- Gamma2[6:7,6:7]
      beta.hat1 <- matrix(unlist(r$posterior.bymodel$mean[1])[7:8],ncol=1)
      beta.hat2 <- matrix(unlist(r$posterior.bymodel$mean[2])[6:7],ncol=1)
      pr1 <- r$bf$postprob[1]
      pr2 <- r$bf$postprob[2]
      Var1 <- (Gamma1.1+(beta.hat1%*%t(beta.hat1)))*pr1
      Var2 <- (Gamma2.2+(beta.hat2%*%t(beta.hat2)))*pr2
      Gamma <- Var1+Var2 - Ebeta%*%t(Ebeta)
      inv.Gamma <- solve(Gamma)
      W <- t(Ebeta)%*%inv.Gamma%*%Ebeta
      pval <- pchisq(W,2,lower.tail=F)
      return(pval)
    } else {
      pval <- NA
      return(pval)
    }

  } else if(length(Cov)==0){
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
      Ebeta <- matrix(c(r$posterior$mean[6],r$posterior$mean[7]),ncol=1)
      Gamma1 <- matrix(r$posterior.bymodel$var[[1]],ncol=ncol(r$posterior.bymodel$var[[1]]))
      Gamma2 <- matrix(r$posterior.bymodel$var[[2]],ncol=ncol(r$posterior.bymodel$var[[2]]))
      Gamma1.1 <- Gamma1[7:8,7:8]
      Gamma2.2 <- Gamma2[6:7,6:7]
      beta.hat1 <- matrix(unlist(r$posterior.bymodel$mean[1])[7:8],ncol=1)
      beta.hat2 <- matrix(unlist(r$posterior.bymodel$mean[2])[6:7],ncol=1)
      pr1 <- r$bf$postprob[1]
      pr2 <- r$bf$postprob[2]
      Var1 <- (Gamma1.1+(beta.hat1%*%t(beta.hat1)))*pr1
      Var2 <- (Gamma2.2+(beta.hat2%*%t(beta.hat2)))*pr2
      Gamma <- Var1+Var2 - Ebeta%*%t(Ebeta)
      inv.Gamma <- solve(Gamma)
      W <- t(Ebeta)%*%inv.Gamma%*%Ebeta
      pval <- pchisq(W,2,lower.tail=F)
      return(pval)
    } else{
      pval <- NA
      return(pval)
    }
  }
}

