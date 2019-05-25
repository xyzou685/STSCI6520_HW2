#' Algorithmic Leveraging Function
#' 
#' implements algorithmic leveraging for linear regression using uniform and leverage score based on subsampling of rows
#' @param xi independent variables
#' @param yi response variable
#' @param r the sampling size,default=20% of the input data
#' @return beta,the estimated regression coefficients
#' @export
#' 
#' @examples
#' x <- matrix(rnorm(500),100,5)
#' y <- x%*%rep(-1,5)+rnorm(100)
#' algo_leverage(x,y)
#'
algo_leverage <- function(xi,yi,r=floor(0.2*length(yi))){
  n=length(yi)
  hii <- diag(xi%*%(t(xi)%*%xi)%*%t(xi))

  unif <- sample(n,r,replace = TRUE)
  blev <- sample(n,r,replace = TRUE,prob =hii/sum(hii))
  if(NCOL(xi)>1){
    xunif <- xi[unif,]
    xblev <- xi[blev,]
  }
  else{
    xunif <- xi[unif]
    xblev <- xi[blev]
  }
  yunif <- yi[unif]
  yblev <- yi[blev]
  b_unif= solve(t(xunif)%*%xunif)%*%t(xunif)%*%yunif
  b_blev= solve(t(xblev)%*%xblev)%*%t(xblev)%*%yblev

  b_est <- list(b_unif,b_blev)
  return(b_est)
}
