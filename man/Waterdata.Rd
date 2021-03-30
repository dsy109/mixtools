\name{Waterdata}
\docType{data}
\title{Water-Level Task Data Set}
\alias{Waterdata}
\alias{WaterdataFull}
\usage{
data(Waterdata)
}
\description{
  This data set arises from the water-level task proposed by the Swiss 
  psychologist Jean Piaget to assess children's understanding of the physical 
  world.  This involves presenting a child with a rectangular vessel with a cap,
  affixed to a wall, that can be tilted (like the minute hand of a clock) to point
  in any direction.  A separate disk with a water line indicated on it, which
  can similarly be spun so that the water line may assume any desired angle with
  the horizontal, is positioned so that by spinning this disk, the child subject 
  may make the hypothetical surface of water inside the vessel assume any desired
  orientation.  For each of eight different orientations of the vessel, corresponding
  to the clock angles at 1:00, 2:00, 4:00, 5:00, 7:00, 8:00, 10:00, and 11:00,
  the child subject is asked to position the water level as it would appear in reality
  if water were in the vessel.  The measurement is the acute angle with 
  the horizontal, in degrees,
  assumed by the water line after it is positioned by the child.  A sign is attached
  to the measurement to indicate whether the line slopes up 
  (positive) or down (negative) from left to right.  Thus, each 
  child has 8 repeated measurements, one for each vessel angle, and the range of
  possible values are from -90 to 90.
  
  The setup of the experiment, along with a photograph of the testing apparatus,
  is given by Thomas and Jamison (1975).  A more detailed analysis using a 
  subset of 405 of the original 579 subjects is given
  by Thomas and Lohaus (1993); further analyses using the functions in
  \code{mixtools} are given by Benaglia et al (2008) and Levine et al (2011), 
  among others.
  
  There are two versions of the dataset included in \code{mixtools}.  The full 
  dataset, called \code{WaterdataFull}, has 579 individuals.  The dataset
  called \code{Waterdata} is a subset of 405 individuals, comprising all children
  aged 11 years or more and omitting any individuals with any observations equal
  to 100, which in this context indicates a missing value (since all of the degree
  measurements should be in the range from -90 to +90, 100 is not a possible value).
}
\format{These data frames consist of 405 or 579 rows, one row for each child. 
  There are ten columns:
  The age (in years) and sex (where 1=male and 0=female) are given for each 
  individual along with the 
  degree of deviation from the horizontal
  for 8 specified clock-hour orientations 
  (11, 4, 2, 7, 10, 5, 1, and 8 o'clock, in order).
}
\source{
  Benaglia, T., Chauveau, D., and Hunter, D.R. (2009), 
  An EM-Like Algorithm for Semi- and Non-Parametric Estimation in 
  Multivariate Mixtures, 
  \emph{Journal of Computational and Graphical Statistics}, 
  18: 505-526.

  Levine, M., Hunter, D.R., and Chauveau, D. (2011), 
  Maximum Smoothed Likelihood for Multivariate Mixtures, 
  \emph{Biometrika}, 
  98(2): 403-416.

  Thomas, H. and Jamison, W. (1975),
  On the Acquisition of Understanding that Still Water is Horizontal,
  \emph{Merrill-Palmer Quarterly of Behavior and Development},
  21(1): 31-44.

  Thomas, H. and Lohaus, A.  (1993),
  \emph{Modeling Growth and Individual Differences in Spatial Tasks},
  University of Chicago Press, Chicago, available on JSTOR.  
} 
\keyword{datasets}
