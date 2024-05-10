\name{distR}
\alias{distR}
\title{
*Before using, please see details - Calculate linear distances in mm from a given set of three-dimentional points
}
\description{
Calculate linear distances (in millimetres) that will subsequently be used for size standardization processes.The difference with CalEstDistances is that this function utilizes a specific set of points designed by Godinho et al., 2020 (DOI: 10.1002/ajpa.24095), which has a specific order and size (21 points entered as text input plus another one that is subsequently added from the R environment, as it has been calculated in previous steps).
}
\usage{
distR(x)
}
\arguments{
  \item{x}{
must be a text string containing all the coordinates of the points to be included in the calculations. Each value must be separated by a space from the next, and the order must follow: X coordinate of point 1, Y coordinate of point 1, Z coordinate of point 1, X coordinate of point 2, Y coordinate of point 2...
}
}
\details{
*Do not use this unless you have coleected your landmarks coordinates following the exact same method of Godinho et al., 2020 (DOI: 10.1002/ajpa.24095)
The reason for including two functions that are so similar in the package is that when designing the protocol, a function should be included to allow following the described workflow in the stard protocol (CalcEstDistances). However, in our personal research case, the workflow was slightly different because we already had the coordinates of the points stored in a different way, so we had to desing a new function that addressed these circumstances (CalcEstDistancesR).
}
\value{
Distance in milimeters of:
- Mandibular corpus length
- Mandibular ramuss lateral length
- Mandibular ramus width
- Dental arch breadth
For further details on the measurements and the reason for their calculation, see references (see below).
}
\references{
Godinho RM, O'Higgins P, Gonçalves C. Assessing the reliability of virtual reconstruction of mandibles. Am J Phys Anthropol. 2020; 172: 723–734. https://doi.org/10.1002/ajpa.24095
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
mirror_A(result_matrix)
pMV(result_matrix_A_Prime)
est_points_text <- "24 17 91 53 82 64 35 49 70 92 18 58 76 29 83 45 67 51 22 11 72 13 55 40 20 95 34 61 87 77 43 28 19 39 86 31 75 52 68 50 16 98 25 63 12 57 14 21 84 47 38 96 30 73 65 36 80 59 48 88 46 90 33"
distR(est_points_text)
}
