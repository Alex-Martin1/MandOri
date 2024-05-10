\name{t2m}
\alias{t2m}
\title{
*Before using, please see details - Transform text inputs into a 4x3 matrix where each row are the coordinates of a specific point
}
\description{
Enables the generation of a 4x3 matrix comprising three-dimensional points (each column representing the X, Y, and Z coordinates of each point) from a text input containing the values of these coordinates. It accepts both scientific notation and decimal notation. The order in which the coordinates of the points are introduced is crucial, as the labels of the points are always generated in the same order, irrespective of whether the coordinates introduced truly match those of that specific point.
}
\usage{
t2m(x)
}
\arguments{
  \item{x}{
must be a text string containing all the coordinates of the 4 points to be included in the matrix. Each value must be separated by a space from the next, and the order must follow: X coordinate of point 1, Y coordinate of point 1, Z coordinate of point 1, X coordinate of point 2, Y coordinate of point 2...
}
}
\details{
*Do not use this unless you have coleected your landmarks coordinates following the exact same method of Godinho et al., 2020 (DOI: 10.1002/ajpa.24095)
}
\value{
A 5x3 matrix with each of the points introduced labelled. 
}
\author{
Alejandro Martín-Sánchez
}
\note{
As a low-level user of R, while I grasp the logic of the language and handle its fundamental principles, I have turned to AI to generate the code for the entire package in order to save hours of work. A seasoned developer or even someone proficient in the field may have easily noticed this, as the code may not be as efficient as it could be. For this reason, any suggestions on how to improve it are more than welcome :)
}
\examples{
require(MandOri)
text <- "2 2 8 5 8 4 6 9 7 6 3 2"
result_matrix <- t2m(text)
result_matrix
}