\name{getmatrix}
\alias{getmatrix}
\title{
Transform text inputs into a 9x3 matrix where each row are the coordinates of a specific point
}
\description{
Enables the generation of a 9x3 matrix comprising three-dimensional points (each column representing the X, Y, and Z coordinates of each point) from a text input containing the values of these coordinates. It accepts both scientific notation and decimal notation. The order in which the coordinates of the points are introduced is crucial, as the labels of the points are always generated in the same order, irrespective of whether the coordinates introduced truly match those of that specific point. To determine the correct order in which these coordinates should be entered to generate the matrix correctly, please refer to the link provided in the "references" section.   
}
\usage{
getmatrix(x)
}
\arguments{
  \item{x}{
must be a text string containing all the coordinates of the 9 points to be included in the matrix. Each value must be separated by a space from the next, and the order must follow: X coordinate of point 1, Y coordinate of point 1, Z coordinate of point 1, X coordinate of point 2, Y coordinate of point 2...
}
}
\details{
The reason for including two functions that are so similar (this one and t2m) in the package is that when designing the protocol, a function should be included to allow following the described workflow in the standard protocol (protocols.io). However, in our personal research case, the workflow was slightly different, so we had to desing a new function that addressed these circumstances (t2m).
}
\value{
A 9x3 matrix with each of the points introduced labelled. 
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