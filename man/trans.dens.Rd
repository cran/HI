\name{trans.dens}
\alias{trans.dens}
\alias{trans.up}
\title{Functions for transdimensional MCMC}
\description{
Computes the value of the 'g' density at a given point and, optionally, 
returns the backtransformed point and the model to which the point
belongs.
}
\usage{
trans.dens(y, ldens.list, which.models, ..., back.transform=F)
trans.up(x, ldens.list, which.models, ...)
}
\arguments{
  \item{y}{Vector or matrix of points (by row) at which the 
  density of the absolutely continuous auxiliary distribution 
  has to be evaluated}
  \item{x}{Vector or matrix of points corresponding to \code{y} (see examples)}
  \item{ldens.list}{List of densities (of submodels)}
  \item{which.models}{List of integers, in increasing order,
    giving the number of components to be dropped when evaluating
    the density in \code{which.models} in the corresponding position.
    A first element equal 0 (full model) is added if not already present}
  \item{back.transform}{Logical that determines the output}
  \item{...}{Other arguments passed to the functions in \code{ldens.list}}
}
\details{See the reference for details}
\value{
   If \code{back.transform=F}, \code{trans.dens} returns the density
   of the absolutely continuous auxiliary distribution evaluated 
   at the point(s) \code{y}.\\
   If \code{back.transform=T}, \code{trans.dens} returns in addition 
   the point \code{x} corresponding to \code{y} in the original space 
   and the index of the subspace to which \code{x} belongs.\\
   \code{trans.up} is a (stochastic) right inverse of the correspondence
   between \code{y} and \code{x}
}
\references{Petris \& Tardella, A geometric approach to transdimensional Markov
  chain Monte Carlo.
   \emph{The Canadian Journal of Statistics}, vol.31, n.4, (2003) (?). }
\author{Giovanni Petris \email{GPetris@uark.edu}, Luca Tardella}

\examples{
#### ==> Warning: running the examples may take a few minutes! <== ####    
### Generate a sample from a mixture of 0,1,2-dim standard normals
ldens.list <- list(f0 = function(x) sum(dnorm(x,log=TRUE)),
                   f1 = function(x) dnorm(x,log=TRUE),
                   f2 = function() 0)
trans.mix <- function(y) {
    trans.dens(y, ldens.list=ldens.list, which.models=0:2)
}

trans.rmix <- arms(c(0,0), trans.mix, function(x) crossprod(x)<1e4, 500)
rmix <- trans.dens(y=trans.rmix, ldens.list=ldens.list,
                      which.models=0:2, back.transform = TRUE)
table(rmix[,2])/nrow(rmix) # should be about equally distributed
plot(trans.rmix,col=rmix[,2]+3,asp=1, xlab="y.1", ylab="y.2",
     main="A sample from the auxiliary continuous distribution")
x <- rmix[,-(1:2)]
plot(x, col=rmix[,2]+3, asp=1,
     main="The sample transformed back to the original space")
### trans.up as a right inverse of trans.dens
y <- trans.up(x, ldens.list, 0:2)
stopifnot(all.equal(x, trans.dens(y, ldens.list, 0:2, back.transform=TRUE)[,-(1:2)]))

### More trans.up
z <- trans.up(matrix(0,1000,2), ldens.list, 0:2)
plot(z,asp=1,col=5) # should look uniform in a circle corresponding to model 2
z <- trans.up(cbind(runif(1000,-3,3),0), ldens.list, 0:2)
plot(z,asp=1,col=4) # should look uniform in a region corresponding to model 1
}
\keyword{distribution}

