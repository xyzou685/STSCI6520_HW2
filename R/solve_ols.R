#' solve_ols
#' 
#' Solve linear system by Gauss-Seidel, Jacobi(sequential), or Jacobi(parallel) method
#' @param A n*n matrix of the input system Ax=b
#' @param b n*1 vector of the input system Ax=b
#' @param cores Number of cores used in parallel computing of Jacobi Mathod. If cores=1,use sequential method, otherwise use parallel computing.
#' @param method Method to use. "GS"=Gauss-Seidel,"Jacobi"=Jacobi
#' @param iteration Number of iterations
#'
#' @import doParallel,foreach
#' @return x,n*1 vector, solution of the input system Ax=b
#' @export
#'
#' @examples
#' n=10
#' L <- diag(0, n)
#' L[(row(L) - col(L)) == 1] <- -1
#' U <- diag(0, n)
#' U[(row(U) - col(U)) == -1] <- -1
#' D <- diag(2, n)
#' a <- L+D+U
#' v <- as.matrix(rep(1,10))
#' b=a%*%v
#' solve_ols(a,b,method = "GS",iteration = 100)
#' solve_ols(a,b,cores=1,method = "Jacobi",iteration=100)
#' solve_ols(a,b,cores=2,method = "Jacobi",iteration=100)

solve_ols <- function(A,b,cores=1,method,iteration){

  n=length(b)
  x0 = as.vector(rep(0,n))
  L <-matrix(0,n,n)
  D <-matrix(0,n,n)
  U <-matrix(0,n,n)
  diag(D) <- diag(A)
  L[lower.tri(L)] <- A[lower.tri(A)]
  U[upper.tri(U)] <- A[upper.tri(A)]
  xi=x0

  if (method=="GS"){
    #Gauss-Seidel
    for (i in 1:iteration){
      xi = (solve(L+D))%*%(b-U%*%xi)
      #xi=solve(L+D)%*%b-solve(L+D)%*%U%*%xi
    }
  }
  else if(method=="Jacobi"){
    if(cores==1){
      #Jacobi (sequential)
      for (i in 1:iteration){
        xi = solve(D)%*%(b-((L+U)%*%xi))
      }
    }
    else{
      #Jacobi (parallel)
      library(doParallel)
      library(foreach)
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      outlist=foreach (i=1:iteration) %dopar% {
        xi = solve(D)%*%(b-(L+U)%*%xi)
      }
      xi=outlist[iteration]
    }}
  return(xi)}
