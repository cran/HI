\name{rballunif}
\alias{rballunif}
\title{Function to generate a random vector uniformly distributed in a sphere} 
\description{
Generates a random vector from a uniform distribution in a sphere
centered at the origin
}
\usage{
rballunif(n,d)
}
\arguments{
   \item{n}{Dimension of the random vector to be generated}
   \item{d}{Radius of the sphere}
}
\details{
   The function is not vectorized: it is intended to generate
   one random vector at a time.
}
\value{
   An pseudo-random vector from a uniform distribution in the
   \code{n}-dimensional sphere centered at the origin and having radius 
   \code{d}
}
\author{Giovanni Petris \email{GPetris@uark.edu}} 
\keyword{distribution}

