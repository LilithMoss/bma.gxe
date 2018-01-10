#' Re-name data and call BMA function for 2-degree-of-freedom Test
#'
#' @param Sample input A data frame or matrix with columns in the order: "Y" "E" "G". Defaults to NULL.
#' @param Covar input A data frame of all covariate variables. Defaults to NULL.
#' @param cc input A number. Defaults to 0.5.
#' @param co input A number. Defaults to 0.5.
#' @param phi input A number. Defaults to 1.
#' @param psi input A number. Defaults to 1000.
#' @param formul input A string in the form: "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+... . Defaults to NULL.
#' @return A number containing the p-value for a BMA joint test of marginal and interaction effect.
#' @examples
#' data("Covariates")
#' data("Sample")
#' formul <- "Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+C1+C2+C3+Y:C1+Y:C2+Y:C3+E:C2"
#' bma.df1(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=NULL)
#' bma.df2(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=NULL)
#' bma.df1(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=formul)
#' bma.df2(Sample=Sample,Covar=Covariates,cc=0.5,co=0.5,phi=1,psi=1000,formul=formul)
#' @export

bma.df2 <- function(Sample=NULL,Covar=NULL,cc=0.5,co=0.5,phi=1,psi=1000,formul=NULL){
  if(cc+co!=1){
    print("error: cc and co values must add up to 1")
  } else {
    if(!missing(Covar) & length(Covar)>0){
      if(nrow(Sample)!=nrow(Covar)){
        print("error: Sample and Covar must have same number of rows")
      } else {
        if(colnames(Sample)[1]!="Y" | colnames(Sample)[2]!= "E" | colnames(Sample)[3] != "G"){
          print("error: names do not match 'Y,E,G'")
        } else {
          Y<-Sample$Y;E<-Sample$E;G<-Sample$G
          names(Covar) <- paste0("C",seq(1,ncol(Covar)))
          pval <- run.BMA.2DF(Y=Y,E=E,G=G,Covar=Covar,cc=cc,co=co,phi=phi,psi=psi,formul=formul)
          return(pval)
        }
      }

    } else {
      if(colnames(Sample)[1]!="Y" | colnames(Sample)[2]!= "E" | colnames(Sample)[3] != "G"){
        print("error: names do not match 'Y,E,G'")
      } else {
        Y<-Sample$Y;E<-Sample$E;G<-Sample$G
        pval <- run.BMA.2DF(Y=Y,E=E,G=G,Covar=NULL,cc=cc,co=co,phi=phi,psi=psi,formul=formul)
        return(pval)
      }
    }
  }
}
