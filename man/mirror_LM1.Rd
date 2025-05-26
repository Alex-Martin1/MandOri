\name{mirror_LM1}
\alias{mirror_LM1}
\title{
Mirroing landmark 1 respect to the plane P.
}
\description{
This function computes the X, Y and Z coordinates of a new landmark (1' or 1_Line) symmetric to the first introduced landmark (1) with respect to a plane defined by the following three landmarks that were introduced as an input (2, 3, 4). Thus, the order in which the coordinates are entered is crucial for the subsequent result. Once the function is executed, it will automatically update the existing "result_matrix" object by adding the new landmark 1_Line as the last row.
}
\usage{
mirror_LM1(result_matrix)
}
\arguments{
  \item{result_matrix}{
must be a matrix where the landmarks 1, 2, 3, and 4 are arranged from top to bottom in that same order.
}
}
\details{
This function can only be executed with matrices containing three columns; i.e., with three-dimensional points where the coordinates for X, Y, and Z are each located in a separate column)
}
\value{
Three-dimentional coordinates of the new calculated point (1_Line).
Besides, the function will update the existing "result_matrix" object by adding the new landmark 1_Line as the last row.
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
text <- "-3.176e+01 -7.557e+00 -1.441e+02 -1.478e+00 6.264e+00 -1.046e+02 -5.153e-02 -2.313e+01 -9.379e+01 -5.700e-01 -1.662e+01 -1.083e+02 -2.924e+01 -5.551e+00 -1.335e+02 -2.191e+01 -5.099e+00 -1.348e+02 -2.557e+01 -1.035e+00 -1.221e+02 -1.980e+01 -6.242e-01 -1.243e+02 -4.616e+01 -3.496e+01 -1.652e+02 -1.790e+01 2.789e+00 -1.106e+02 -4.694e+01 1.235e+01 -2.013e+02"
result_matrix <- getmatrix(text)
mirror_LM1(result_matrix)
View(result_matrix)
}
