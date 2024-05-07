\name{PointDistancemm}
\alias{PointDistancemm}
\title{
Calculate the distance in millimeters between two points
}
\description{
This function calculate the distance in millimeters between two specific points. These points must be present in the R environment for the function to execute correctly
}
\usage{
PointDistancemm(x, y)
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
Existing distance in milimeters between the two input points.
}
\author{
Alejandro Martín-Sánchez
}
\note{
As a low-level user of R, while I grasp the logic of the language and handle its fundamental principles, I have turned to AI to generate the code for the entire package in order to save hours of work. A seasoned developer or even someone proficient in the field may have easily noticed this, as the code may not be as efficient as it could be. For this reason, any suggestions on how to improve it are more than welcome :)
}
\examples{
require(MandOri)
Point_A <- "1 2 3"
Point_B <- "5 5 3"
point_A <- text2vector(Point_A)
point_B <- text2vector(Point_B)
PointDistancemm(point_A, point_B)
}
