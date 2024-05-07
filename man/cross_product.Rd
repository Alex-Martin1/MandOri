\name{cross_product}
\alias{cross_product}
\title{
Computes the result of the cross product between two specific vectors
}
\description{
This function computes the result of the cross product between two input vectors. These vectors must be present in the R environment for the function to execute correctly
}
\usage{
cross_product(x, y)
}
\arguments{
  \item{x}{
x must be an R vector consisting of three distinct numerical values, each separated from the others by a space.
}
  \item{y}{
y must be another R vector composed of three distinct numerical values, each separated from the others by a space.
}
}
\value{
Cross product of the two input vectors.
}
\author{
Alejandro Martín-Sánchez
}
\note{
As a low-level user of R, while I grasp the logic of the language and handle its fundamental principles, I have turned to AI to generate the code for the entire package in order to save hours of work. A seasoned developer or even someone proficient in the field may have easily noticed this, as the code may not be as efficient as it could be. For this reason, any suggestions on how to improve it are more than welcome :)
}
\examples{
require(MandOri)
u_vector <- "5 8 4"
v_vector <- "2 5 8"
u <- text2vector(u_vector)
v <- text2vector(v_vector)
cross_product(u, v)
{
  }
}
