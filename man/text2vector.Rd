\name{text2vector}
\alias{text2vector}
\title{
Allows transforming a text input (string) containing numeric characters into a R vector (a list of items of the same type)
}
\usage{
text2vector(x)
}
\arguments{
  \item{x}{
x should be a text string containing values separated by a space from the nexts
}
}
\value{
R vector (A list of items that R conceives as the same type)
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
vector <- text2vector(text)
}
