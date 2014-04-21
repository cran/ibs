\name{ibs}
\alias{ibs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrated B-spline
}
\description{
Calculate the integral of a B-spline function.
}
\usage{
ibs(x, knots, ord, coef = rep(1, length(knots) - ord))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Numerical value or vector. The value(s) at which to evaluate the integral of the B-spline; must
    be in the bewteen \code{min(knots)} and \code{max(knots)}. 
}
  \item{knots}{
    Numerical vector. The knot positions/sites of the B-spline function to be integrated.
}
  \item{ord}{
    An integer >=1. The order of the B-spline integrand function to be integrated. Equals degree plus 1.
}
  \item{coef}{
   A numerical vector. The coefficients (de Boor points) defining the B-spline integrand function.
}
}
\details{
  The function returns the integral(s) of the B-spline function
  specified  by knots \code{knots}, order \code{ord}, and coefficients
  \code{coef}, from the minimum knot position to each \code{x}
  value. The evaluation is based on a closed form expression of the
  integral in terms of higher order B-splines, given on page 128 of de
  Boor (2001).  
}
\value{
  A numerical equal to the integral(s).
}
\references{
  de Boor, C (2001) A Practical Guide to Splines. Revised Edition. Springer: New York.
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
% \note{
%   %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{bspline}}
}
\examples{
kns <- c(rep(0,4),1:4*0.2,rep(1,4))
co <- rnorm(length(kns)-3)
integrate(bspline,knots=kns,ord=3,coef=co,0,0.75)
require(splines)
integrate(function(x)splineDesign(kns,x,3) \%*\% co,0,0.75)
ibs(0.75,kns,3,co)

set.seed(2);
nk <- 20; o <- 4;
kns <- sort(rnorm(nk))
co <- rnorm(nk-o)
op <- par(mfrow=c(2,1))
curve(ibs(x,kns,o,co),kns[1],kns[nk],n=1001,xlab="x",ylab="ibs")
bspl <- curve(bspline(x,kns,o,co),kns[1],kns[nk],n=1001,
              xlab="x",ylab="B-spline")
abline(h=0,col="gray",lty=2)
par(op)
op <- par(new=TRUE)
plot(0,0,xlim=range(kns),xaxt="n",yaxt="n",type="n",xlab="",ylab="",bty="n")
abline(v=bspl$x[abs(bspl$y)<2e-2],col="gray")
par(op)

}
\keyword{math}
