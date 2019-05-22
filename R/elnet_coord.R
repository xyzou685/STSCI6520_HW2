#' elnet_coord
#' ﬁts elastic net to data using coordinate descent algorithm.
#'
#' @param x n*p matrix;explanatory variables
#' @param y n*1 vector;response variable
#' @param a a constant;the l-2 regularization coefficient
#' @param lambda a constant;the l-1 and l-2 regularization coefficient
#' @param betahat p*1 vector;the initial regression coefficient;default is a zero vector
#' @param maxiteration maximum number of iteration allowed
#'
#' @return betahat the estimated regression coefficient
#' @export
#'
#' @examples 
#' x <- matrix(rnorm(500),100,5)
#' y <- x%*%rep(-1,5)+rnorm(100)
#' elnet_coord(x,y,0.5,1)
elnet_coord <- function(x,y,a,lambda,betahat=rep(0,NCOL(x)),maxiteration=100){
  for (i in 1:maxiteration) {
    beta1 <- betahat
    for (j in 1:NCOL(x)) {
      r <- y-x%*%beta1+x[,j]*beta1[j]
      num <- 2*(x[,j]%*%r)/n-a*lambda
      denom <- 2*(sum(x[,j]^2)/n+(1-a)*lambda)
      if (num==0){
        beta1[j]=0
      }
      else{
        #solve beta1[j] by equating derivatives of minimizing function to zero. There're 2 cases:
        #1. beta1>0; 2. beta1<0
        beta1[j]=ifelse(num>0,num/denom,(num+2*a*lambda)/denom)
      }
    }
    if(sum(abs(beta1-betahat))<=10^(-7)){
      break
    }
    betahat <- beta1
  }
  return(betahat)
}
