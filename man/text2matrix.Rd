\name{text2matrix}
\alias{text2matrix}
\title{
Transform text inputs into matrices
}
\description{
Enables generating a matrix of K three-dimensional points (Kx3 matrix) from a text input containing the values of these coordinates. It accepts both scientific notation and decimal notation.
}
\usage{
text2matrix(x)
}
\arguments{
  \item{x}{
x must be a text string containing all the coordinates of the points to be included in the matrix. Each value must be separated by a space from the next, and the order must follow: X coordinate of point 1, Y coordinate of point 1, Z coordinate of point 1, X coordinate of point 2, Y coordinate of point 2...
}
}
\value{
A Kx3 matrix, where K (the number of rows of the matrix) is the number of three-dimensional points included ([number of individual values sepated by spaces] divided by 3)
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
text <- "2 2 8 5 8 4 6 9 7 6 3 2"
result_matrix <- text2matrix(text)
result_matrix
}
