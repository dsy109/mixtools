\name{RTdata2}
\docType{data}
\title{Reaction Time (RT) Data Set (No. 2)} 
\alias{RTdata2}
\usage{
data(RTdata2)
}
\description{
  This data set involves normally developing children 9 years of age presented 
  visual simuli on a computer monitor.  There are three different experimental
  conditions, according to the length of the delay after which the stimulus 
  was displayed on the screen.  Each subject experienced each condition eight
  times, and these 24 trials were given in random order.
  These data give the 82 children for whom
  there are complete measurements among over 200 total subjects.
}
\format{This data frame consists of 82 children (the rows) and their 24 responses (the columns)
to the stimulus presented.  The response is recorded in milliseconds.  The columns are
not in the order in which the stimuli were presented to the children;  rather, they
are arranged into three blocks of eight columns each so that each eight-column
block contains only trials from one of the three conditions.
}
\references{
  Miller, C. A., Kail, R., Leonard, L. B. and Tomblin, J. B. (2001) Speed of Processing in Children with
  Specific Language Impairment, \emph{Journal of Speech, Language, and Hearing Research} \bold{44(2)}, 416--433.
}
\seealso{
\code{\link{RTdata}}
}
\keyword{datasets}
