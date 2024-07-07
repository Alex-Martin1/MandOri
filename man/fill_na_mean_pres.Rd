\name{fill_na_mean_pres}
\alias{fill_na_mean_pres}
\title{
Function that fills NA cells in a column preserving the mean of the previusly existing values and generates a new dataset with NA cells filled
}
\description{
This function calculates the mean ignoring the NA values and then generates X values, where X is the number of NA cells that exist within the selected column. Those newly generated values will always be between the maximun and minimun pre-existing values. Finally, it generates a new dataset with NA cells filled.    
}
\usage{
fill_na_mean_pres(df, "column")
}
\arguments{
  \item{df}{
must be a previuosly uploaded tidy dataframe
}
\item{column}{
must be a column of the uploaded dataframe
}
}
\author{
Alejandro Martín-Sánchez
}
\note{
As a low-level user of R, while I grasp the logic of the language and handle its fundamental principles, I have turned to AI to generate the code for the entire package in order to save hours of work. A seasoned developer or even someone proficient in the field may have easily noticed this, as the code may not be as efficient as it could be. For this reason, any suggestions on how to improve it are more than welcome :)
}



