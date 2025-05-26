\name{getmatrix}
\alias{getmatrix}
\title{
Convert string text inputs into an R×3 numeric matrix.
}
\description{
Enables the generation of a Rx3 matrix, where R is computed by dividing the total number of space-separated numeric values by three. This structure assumes that the input encodes a sequence of three-dimensional coordinates, with each row representing the x, y, and z values of a single point in 3D space. 
It accepts both scientific notation and decimal notation. 
The order in which the coordinates of the points are introduced is crucial, as the labels of the points are always generated in the same order, irrespective of whether the coordinates introduced truly match those of that specific point. To determine the correct order in which these coordinates should be entered to generate the matrix correctly, please refer to the link provided in the "references" section.   
}
\usage{
getmatrix(text)
}
\arguments{
  \item{text}{
must be a text string containing all the coordinates of the points to be included in the matrix. Each value must be separated by a space from the next, and the order must follow: X coordinate of point 1, Y coordinate of point 1, Z coordinate of point 1, X coordinate of point 2, Y coordinate of point 2...
}
}
\details{
}
\value{
A Rx3 matrix with each of the points introduced labelled. 
}
\references{
Protocols.io - to be published
}
\author{
Alejandro Martín-Sánchez
}
\note{
As a low-level user of R, while I grasp the logic of the language and handle its fundamental principles, I have turned to AI to generate the code for the entire package in order to save hours of work. A seasoned developer or even someone proficient in the field may have easily noticed this, as the code may not be as efficient as it could be. For this reason, any suggestions on how to improve it are more than welcome :)
}
\examples{
require(MandOri)
text <- "2 2 8 5 8 4 6 9 7 6 3 2 5 9 8 5 4 8 5 9 4 5 4 7 4 5 8"
result_matrix <- getmatrix(text)
result_matrix
}