% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/21_cardioModel.R
\name{cardioModel}
\alias{cardioModel}
\title{Cardiovascular safety exposure-response modeling in early-phase clinical studies}
\usage{
cardioModel(x, ID = "ID", NTAFD = "NTAFD", TOD = "TOD", CLOCK = "CLOCK",
  RESPONSE = "RESPONSE", EXPOSURE = "EXPOSURE", PERIOD = "PERIOD",
  SEQUENCE = "SEQUENCE", COHORT = "COHORT", AIC.k = 2,
  delay.confidence = 0.99, study.name = "unknown study",
  drug.name = "unknown drug")
}
\arguments{
\item{x}{A dataframe with one measure of cardiovascular response (e.g., QT interval, heart rate, pulse rate, blood pressure), one measure of drug exposure (e.g., plasma concentration).}

\item{ID}{The name of the column (in x) containing the subject unique identification number.}

\item{NTAFD}{The name of the column containing the nominal time after first dose.}

\item{TOD}{(required if CLOCK is not supplied) The name of the column containing the time of the day, as a continuous variable with range [0, 24).  If both TOD and CLOCK are provided, then only TOD will be used in the model fitting.}

\item{CLOCK}{(required if TOD is not supplied) The name of the column containing the time of the day, as 24 hour clock time [0:00, 24:00).  If both TOD and CLOCK are provided, then only TOD will be used in the model fitting.}

\item{RESPONSE}{The name of the column containing a measure of cardiovascular response (e.g., QT interval, heart rate, pulse rate, blood pressure).}

\item{EXPOSURE}{The name of the column containing a measure of drug exposure (e.g., plasma concentration).}

\item{PERIOD}{(optional) The name of the column specifying the study period for studies where the subjects receive different doses at different occasions.  If a column with this name does not exist, then it is assumed that there is only a single occasion.  In this case, only models without between occasion variability will be fit.}

\item{SEQUENCE}{(optional) The name of the column specifying the sequence when multiple measurements of the response are taken with the same recorded value of NTAFD.  If a column with this name does not exist, then it is assumed that multiple measurements are not recorded for the same subject at the same time.}

\item{COHORT}{(optional)    The name of the column specifying the study cohort when period numbers are nested within cohorts. If a column with this name does not exist, then it is assumed that the period numbers are not nested within cohorts.}

\item{AIC.k}{numeric, the penalty per parameter to be used for computing the AIC; the default here is 2, which is the classical AIC.}

\item{delay.confidence}{numeric in the range (0, 1), the one-sided confidence level threshold for the linear regression analysis performed in the diagnosis for the anticlockwise hysteresis. The default here is 0.99.}

\item{study.name}{(optional) The name of the study.  This is used to label the output.}

\item{drug.name}{(optional) The name of the drug.  This is used to label the output.}
}
\value{
A list with two elements, "summary" and "models". \cr\cr
The first element in the list ("summary") is a data frame with a row for each of the tested models, ordered by ascending values of the Akaike Information Criterion. \cr\cr
The first column of the data frame is named "MODEL", and it provides a brief description of the model formula, while the second column ("AIC") provides Akaike's Information Criterion for the fit. \cr\cr
The next 11 columns ("SLOPE", "VAR.SLOPE", "EMAX", "VAR.EMAX", "EC50", "VAR.EC50", "COV.EMAX.EC50", "HILL", "VAR.HILL", "COV.EMAX.HILL", "COV.EC50.HILL") contain the parameter estimates and associated uncertainties.
If a parameter is not present in a particular model, the corresponding value in the data frame will be NA. For models that include a linear drug effect, the "SLOPE" column shows the estimate for the slope, and the "VAR.SLOPE" column shows the estimated variance as computed by the stats::vcov() function.
Likewise, for models that include an Emax or sigmoidal Emax drug effect, the associated variable estimates are in the columns "EMAX", "EC50", and "HILL".  The components of the estimated covariance matrix are in the columns "VAR.EMAX", "VAR.EC50", "VAR.HILL", "COV.EMAX.EC50", "COV.EMAX.HILL", "and COV.EC50.HILL".  \cr\cr
The second last column of the data frame (DRUG.EFFECT.DELAY) is a diagnosis for anticlockwise hysteresis which is performed to evaluate the appropriateness of the direct-link exposure-response analysis. The first derivative of the individual drug concentration is numerically calculated with respect to time. The derivative (x-axis) is plotted versus standardized (Pearson) residuals from the nlme (y-axis).
A linear regression analysis is performed and the upper confidence bound for the slope is calculated. If the higher confidence bound for the slope is less than zero, then the presence of a drug effect delay is inferred, and the value will be "yes"; otherwise, the value will be "no". \cr\cr
The last column of the data frame (ERROR.MESSAGE) stores the error (if any) that occurred while fitting the model.  If no error occurred, then the value will be "none". \cr\cr
The second element in the list ("models") is a list of nlme objects corresponding to each of the models.  The nlme objects are named with the same names as in the "MODEL" column of "summary" data frame.  The order of the nlme objects is
the same as in the "summary" data frame. \cr\cr
}
\description{
Includes 110 mixed-effects model structures describing the relationship between drug concentration and QT interval, heart rate/pulse rate or blood pressure. Given an exposure-response dataset, the tool fits each model structure to the observed data.
}
\note{
A complete description of each model is provided in Conrado et al. listed in the references.  The MODEL column of the output contains only a very brief description of the model.  For each model, the description starts with a number from 1 to 110 to label the model.  Next, the description shows the form of the dose-response relationship.  Note that models 101 through 110 do not have any drug effect.  In addition, the description lists the terms allowed to have between-subject variability (BSV) and between-occasion variability (BOV).  \cr\cr
Each model is fitted using the the nlme() function.  Thus, any error listed in the ERROR.MESSAGE column originates from the nlme() function or from the subsequent extraction of the fitted parameters. \cr\cr
When fitting the nlme models, the initial parameter estimates are set as follows.
For models with a linear component, the intercept is initially set to the mean of the response for the pooled observations.
The slope is initially set to 0.  For models with a cosine component, the
amplitude is initially set to equal the standard deviation of the response for the pooled observations.
These models also set the initial phase shift of the cosine terms to 0.
For models that include an Emax drug effect, Emax is initially set to equal the standard deviation of the response for the pooled observations,
while the EC50 is initially set to equal the mean of the pooled exposures.
For models that include a Hill coefficient, this is initially set to be 1.\cr\cr
}
\examples{
head(sim_QTcf)
result <- cardioModel(sim_QTcf)
names(result)
mySummary <- result[[1]]
head(mySummary)
}
\references{
Conrado DJ, Chen D, Denney WS. Cardiovascular Safety Assessment in Early-Phase Clinical Studies: A Meta-Analytical Comparison of Exposure-Response Models. [Submitted] \cr\cr
Conrado DJ, Hather GJ, Chen D, Denney WS. Facilitating Exposure-Response Analysis of Cardiovascular Safety Markers in Early-Phase Clinical Studies with the cardioModel Package for R. American Conference on Pharmacometrics (ACoP6), 2015. \url{http://metrumrg.com/assets/pubs/Conrado_ACoP_2015.pdf}
}

