#' Re-name data and call BMA function for 1-degree-of-freedom Test
#'
#' @param Sample input A data frame or matrix with columns in the order: "Y" "E" "G". Defaults to NULL.
#' @param Covar input A data frame of all covariate variables. Defaults to NULL.
#' @param cc input A number. Defaults to 0.5.
#' @param co input A number. Defaults to 0.5.
#' @param phi input A number. Defaults to 1.
#' @param psi input A number. Defaults to 1000.
#' @param formul input A string in the form: "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+... . Defaults to NULL.
#' @return Vector containing BMA interaction effect estimate, standard error, Z-score, and p-value.
#' @examples
#' data("Covariates")
#' data("Sample")
#' formul <- "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+C1+C2+C3+Y:C1+Y:C2+Y:C3+E:C2"
#' bma.df1(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=NULL)
#' bma.df2(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=NULL)
#' bma.df1(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=formul)
#' bma.df2(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=formul)
#' @export

bma.df1 <- function(Sample=NULL,Covar-NULL,cc=0.5,co=0.5,phi=1,psi=1000,formul=NULL){
  if(cc+co!=1){
    print("error: cc and co values must add up to 1")
  } else {
    # if(exists("Covar")){
    # if(!missing(Covar)){
    if(!missing(Covar) & length(Covar)>0){
      if(nrow(Sample)!=nrow(Covar)){
        print("error: Sample and Covar must have same number of rows")
      } else {
        if(colnames(Sample)[1]!="Y" | colnames(Sample)[2]!= "E" | colnames(Sample)[3] != "G"){
          print("error: names do not match 'Y,E,G'")
        } else {
          Y<-Sample$Y;E<-Sample$E;G<-Sample$G
          names(Covar) <- paste0("C",seq(1,ncol(Covar)))
          stats <- t(as.matrix(run.BMA.1DF(Y=Y,E=E,G=G,Covar=Covar,cc=cc,co=co,phi=phi,psi=psi,formul=formul)))
          colnames(stats) <- c("Int.est","Int.sd","Z.score","p.value")
          return(stats)
        }
      }

    } else {
        if(colnames(Sample)[1]!="Y" | colnames(Sample)[2]!= "E" | colnames(Sample)[3] != "G"){
          print("error: names do not match 'Y,E,G'")
        } else {
          Y<-Sample$Y;E<-Sample$E;G<-Sample$G
          # stats <- t(as.matrix(run.BMA.1DF(Y=Y,E=E,G=G,Covar=NULL,cc=cc,co=co,phi=phi,psi=psi,formul=formul)))
          stats <- t(as.matrix(run.BMA.1DF(Y=Y,E=E,G=G,cc=cc,co=co,phi=phi,psi=psi,formul=formul)))
          colnames(stats) <- c("Int.est","Int.sd","Z.score","p.value")
          return(stats)
        }
      }
  }
}

