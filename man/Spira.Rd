\name{Spira}
\alias{Spira}
\docType{data}
\title{ Microarray experiment with smokers}
\description{
Microarray experiment lung tissue for 3 groups: smokers (C: 34), former smokers (F: 18)
and never smokers (N: 23), n=75, p=9968.
}
\usage{data(Spira)}
\format{
A list with three named components Patients, X, and genes.
The component Patients is a dataframe containing covariates about the 75 patients
who particpated. The second component X contains the microarray data corresponding
to the 75 patients.  The matrix X has 75 rows and 9968 columns. 
The vector genes indicates the 487 genes (columns) in X that
were selected using SCRDA. 
}
\details{
The first component of Spira is a dataframe which is named Patients. The components are:
\itemize{
    \item{\code{PATIENT_ID}}{a numeric vector}
    \item{\code{STATUS}}{a factor with levels \code{C} \code{F} \code{N}
    corresponding to current smoker, former smoker and never smoker
    }
    \item{\code{START_YEAR}}{a numeric vector}
    \item{\code{END_YEAR}}{a numeric vector}
    \item{\code{PACK_YEARS}}{a numeric vector}
    \item{\code{AGE}}{a numeric vector}
    \item{\code{SEX}}{a factor with levels \code{F} \code{M}}
    \item{\code{RACIAL}}{a factor with levels \code{AFA} \code{AMI} \code{ASI} \code{CAU} \code{HIS} \code{OTH}}
    \item{\code{y}}{a factor with levels \code{smoker} \code{nonsmoker}}
}
}

\source{
 \url{http://pulm.bumc.bu.edu/aged/index.html}
}

\references{
Avrum Spira, Jennifer Beane, Vishal Shah, Gang Liu, Frank Schembri, Xuemei Yang, John Palma, and Jerome S. Brody
Effects of cigarette smoke on the human airway epithelial cell transcriptome
PNAS, Jul 2004; 101: 10143 - 10148. 

\url{http://www.pnas.org/cgi/citemap?id=pnas;101/27/10143}

}


\examples{
#compute mean expression for each gene
data(Spira)
attach(Spira)
apply(X, 2, mean)
#list names of covariates
names(Patients)
}
\keyword{datasets}


