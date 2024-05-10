\name{mirror_A}
\alias{mirror_A}
\title{
Mirroing point A respect to the P plane.
}
\description{
This function computes the X, Y and Z coordinates of a new point (A' or A_Line) symmetric to the first introduced point (A) with respect to a plane defined by the following three points that were introduced as an input (B, C, D). Thus, the order in which the coordinates are entered is crucial for the subsequent result. Once the function is executed, it will automatically update the existing "result_matrix" object by adding the new point A_Line as the last row.
}
\usage{
mirror_A(x)
}
\arguments{
  \item{x}{
must be a matrix where the points A, B, C, and D are arranged from top to bottom in that same order.
}
}
\details{
This function can only be executed with matrices containing three columns; i.e., with three-dimensional points where the coordinates for X, Y, and Z are each located in a separate column)
}
\value{
Three-dimentional coordinates of the new calculated point (A_Line).
Besides, the function will update the existing "result_matrix" object by adding the new point A_Line as the last row.
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
require(MandiOri)
text <- "2 2 8 5 8 4 6 9 7 6 3 2 5 9 8 5 4 8 5 9 4 5 4 7 4 5 8"
result_matrix <- getmatrix(text)
mirror_A(result_matrix)
View(result_matrix_A_Prime)
}
