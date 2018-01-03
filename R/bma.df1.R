#' Re-name data and call BMA function for 1-degree-of-freedom Test
#'
#' @param Sample input A data frame or matrix with columns in the order: "Y" "E" "G". Defaults to NULL.
#' @param Cov input A data frame of all covariate variables. Defaults to NULL.
#' @param cc input A number. Defaults to 0.5.
#' @param co input A number. Defaults to 0.5.
#' @param phi input A number. Defaults to 1.
#' @param psi input A number. Defaults to 1000.
#' @param formula input A string in the form: "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+... . Defaults to NULL.
#' @return Vector containing BMA interaction effect estimate, standard error, Z-score, and p-value.
#' @examples
#' data("Cov")
#' data("Sample")
#' formula <- "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+C1+C2+C3+Y:C1+Y:C2+Y:C3+E:C2"
#' bma.df1(Sample=Sample,Cov=Cov,cc=0.5,co=0.5,phi=1,psi=1000,formula=NULL)
#' bma.df2(Sample=Sample,Cov=Cov,cc=0.5,co=0.5,phi=1,psi=1000,formula=NULL)
#' bma.df1(Sample=Sample,Cov=Cov,cc=0.5,co=0.5,phi=1,psi=1000,formula=formula)
#' bma.df2(Sample=Sample,Cov=Cov,cc=0.5,co=0.5,phi=1,psi=1000,formula=formula)
#' @export

bma.df1 <- function(Sample=NULL,Cov=NULL,cc=0.5,co=0.5,phi=1,psi=1000,formula=NULL){
  if(cc+co!=1){
    print("error: cc and co values must add up to 1")
  } else {
    if(length(Cov)>0){
      if(nrow(Sample)!=nrow(Cov)){
        print("error: Sample and Cov must have same number of rows")
      } else {
        if(colnames(Sample)[1]!="Y" | colnames(Sample)[2]!= "E" | colnames(Sample)[3] != "G"){
          print("error: names do not match 'Y,E,G'")
        } else {
          Sample[Sample != 0 & Sample != 1] <- NA
          Sample <- Sample[complete.cases(Sample), ]
          Y<-Sample$Y;E<-Sample$E;G<-Sample$G
          names(Cov) <- paste0("C",seq(1,ncol(Cov)))
          stats <- t(as.matrix(run.BMA.1DF(Y=Y,E=E,G=G,Cov=Cov,cc=cc,co=co,phi=phi,psi=psi,formula=formula)))
          colnames(stats) <- c("Int.est","Int.sd","Z.score","p.value")
          return(stats)
        }
      }

    } else {
        if(colnames(Sample)[1]!="Y" | colnames(Sample)[2]!= "E" | colnames(Sample)[3] != "G"){
          print("error: names do not match 'Y,E,G'")
        } else {
          Sample[Sample != 0 & Sample != 1] <- NA
          Sample <- Sample[complete.cases(Sample), ]
          Y<-Sample$Y;E<-Sample$E;G<-Sample$G
          stats <- t(as.matrix(run.BMA.1DF(Y=Y,E=E,G=G,Cov=NULL,cc=cc,co=co,phi=phi,psi=psi,formula=formula)))
          colnames(stats) <- c("Int.est","Int.sd","Z.score","p.value")
          return(stats)
        }
      }
  }
}

