#' @name cardioModel
#' @title Cardiovascular safety exposure-response modeling in early-phase clinical studies
#' @description Includes 100 mixed-effects model structures describing the relationship between drug concentration and QT interval, heart rate/pulse rate or blood pressure. Given an exposure-response dataset, the tool fits each model structure to the observed data.
#' @param x A dataframe with one measure of cardiovascular response (e.g., QT interval, heart rate, pulse rate, blood pressure), one measure of drug exposure (e.g., plasma concentration).
#'
#' @param ID The name of the column (in x) containing the subject unique identification number.
#' @param NTAFD The name of the column containing the nominal time after first dose.
#' @param TOD The name of the column containing the time of the day, as a continuous variable with range [0, 24).
#' @param RESPONSE The name of the column containing a measure of cardiovascular response (e.g., QT interval, heart rate, pulse rate, blood pressure).
#' @param EXPOSURE The name of the column containing a measure of drug exposure (e.g., plasma concentration).
#'
#' @param PERIOD (optional) The name of the column specifying the study period for studies where the subjects receive different doses at different occasions.  If a column with this name does not exist, then it is assumed that there is only a single occasion.  In this case, only models without between occasion variability will be fit.
#' @param SEQUENCE (optional) The name of the column specifying the sequence when multiple measurements of the response are taken with the same recorded value of NTAFD.  If a column with this name does not exist, then it is assumed that multiple measurements are not recorded for the same subject at the same time.
#' @param COHORT (optional)	The name of the column specifying the study cohort when period numbers are nested within cohorts. If a column with this name does not exist, then it is assumed that the period numbers are not nested within cohorts.
#' @param AIC.k numeric, the penalty per parameter to be used for computing the AIC; the default here is 2, which is the classical AIC.
#' @param study.name (optional) The name of the study.  This is used to label the output.
#' @param drug.name (optional) The name of the drug.  This is used to label the output.
#'
#' @return A data frame with a row for each of the tested models, ordered by ascending values of the Akaike Information Criterion. \cr\cr
#' The first column of the data frame is named "MODEL", and it provides a brief description of the model formula, while the second column ("AIC") provides Akaike's Information Criterion for the fit. \cr\cr
#' The next 11 columns ("SLOPE", "VAR.SLOPE", "EMAX", "VAR.EMAX", "EC50", "VAR.EC50", "COV.EMAX.EC50", "HILL", "VAR.HILL", "COV.EMAX.HILL", "COV.EC50.HILL") contain the parameter estimates and associated uncertainties.
#' If a parameter is not present in a particular model, the corresponding value in the data frame will be NA. For models that include a linear drug effect, the "SLOPE" column shows the estimate for the slope, and the "VAR.SLOPE" column shows the estimated variance as computed by the stats::vcov() function.
#' Likewise, for models that include an Emax or sigmoidal Emax drug effect, the associated variable estimates are in the columns "EMAX", "EC50", and "HILL".  The components of the estimated covariance matrix are in the columns "VAR.EMAX", "VAR.EC50", "VAR.HILL", "COV.EMAX.EC50", "COV.EMAX.HILL", "and COV.EC50.HILL".  \cr\cr
#' The second last column of the data frame (DRUG.EFFECT.DELAY) is a diagnosis for anticlockwise hysteresis which is performed to evaluate the appropriateness of the direct-link exposure-response analysis. The first derivative of the individual drug concentration is numerically calculated with respect to time. The derivative (x-axis) is plotted versus standardized (Pearson) residuals from the nlme (y-axis).
#' A linear regression analysis is performed and the 99 percent lower confidence bound for the slope is calculated. If the lower confidence bound for the slope is greater than zero, then the presence of a drug effect delay is inferred, and the value will be "yes"; otherwise, the value will be "no". \cr\cr
#' The last column of the data frame (ERROR.MESSAGE) stores the error (if any) that occurred while fitting the model.  If no error occurred, then the value will be "none". \cr\cr
#'
#' @note A complete description of each model is provided in Conrado et al. listed in the references.  The MODEL column of the output contains only a very brief description of the model.  For each model, the description starts with a number from 1 to 100 to label the model.  Next, the description shows the form of the dose-response relationship.  In addition, the description lists the terms allowed to have between-subject variability (BSV) and between-occasion variability (BOV).  \cr\cr
#' Each model is fitted using the the nlme() function.  Thus, any error listed in the ERROR.MESSAGE column originates from the nlme() function or from the subsequent extraction of the fitted parameters. \cr\cr
#' When fitting the nlme models, the initial parameter estimates are set as follows.
#' For models with a linear component, the intercept is initially set to the mean of the response for the pooled observations.
#' The slope is initially set to 0.  For models with a cosine component, the
#' amplitude is initially set to equal the standard deviation of the response for the pooled observations.
#' These models also set the initial phase shift of the cosine terms to 0.
#' For models that include an Emax drug effect, Emax is initially set to equal the standard deviation of the response for the pooled observations,
#' while the EC50 is initially set to equal the mean of the pooled exposures.
#' For models that include a Hill coefficient, this is initially set to be 1.\cr\cr
#'
#' @examples
#' head(sim_QTcf)
#' result <- cardioModel(sim_QTcf)
#' head(result)
#'
#' @references Conrado DJ, Chen D, Denney WS. Cardiovascular Safety Assessment in Early-Phase Clinical Studies: A Meta-Analytical Comparison of Exposure-Response Models. [To be submitted] \cr\cr
#' Conrado DJ, Hather GJ, Chen D, Denney WS. Facilitating Exposure-Response Analysis of Cardiovascular Safety Markers in Early-Phase Clinical Studies with the cardioModel Package for R. American Conference of Pharmacometrics (ACoP6), 2015. \url{http://metrumrg.com/assets/pubs/Conrado_ACoP_2015.pdf}
#' @export
#' @import nlme plyr
#' @importFrom stats AIC coef fitted lm na.omit residuals sd vcov

cardioModel <- function(x, ID = "ID", NTAFD = "NTAFD", TOD = "TOD", RESPONSE = "RESPONSE", EXPOSURE = "EXPOSURE", PERIOD = "PERIOD", SEQUENCE = "SEQUENCE", COHORT = "COHORT", AIC.k = 2, study.name = "unknown study", drug.name = "unknown drug"){

  # remove extra columns
  x <- x[, colnames(x) %in% c(ID, NTAFD, TOD, RESPONSE, EXPOSURE, PERIOD, SEQUENCE, COHORT)]

  # rename required columns
  x <- rename_column(x, ID, "ID", required = TRUE)
  x <- rename_column(x, NTAFD, "NTAFD", required = TRUE)
  x <- rename_column(x, TOD, "TOD", required = TRUE)
  x <- rename_column(x, RESPONSE, "RESPONSE", required = TRUE)
  x <- rename_column(x, EXPOSURE, "EXPOSURE", required = TRUE)

  # rename optional columns
  x <- rename_column(x, PERIOD, "PERIOD", required = FALSE)
  x <- rename_column(x, SEQUENCE, "SEQUENCE", required = FALSE)
  x <- rename_column(x, COHORT, "COHORT", required = FALSE)

  myData <- x

  # assign classes
  myData  <- myPremodeling(myData)
  
  # workaround for bug in nlme
  extra_e = environment()
  with(extra_e, {
  
    # no BOV
    results10 <- myModeling10(myData, AIC.k)
    results11 <- myModeling11(myData, AIC.k)
    results12 <- myModeling12(myData, AIC.k)
    results13 <- myModeling13(myData, AIC.k)
    
    if("PERIOD" %in% colnames(myData)){
      # with BOV
      results14 <- myModeling14(myData, AIC.k)
      results15 <- myModeling15(myData, AIC.k)
      results16 <- myModeling16(myData, AIC.k)
      results17 <- myModeling17(myData, AIC.k)
      out <- mySummary(results=list(results10, results11, results12, results13,
                                    results14, results15, results16, results17))
    } else {
      out <- mySummary(results=list(results10, results11, results12, results13))
    }
    
  # workaround for bug in nlme
  })
  
  # return result
  out$DRUG <- drug.name
  out$STUDY <- study.name
  out <- out[, c(ncol(out), ncol(out)-1, 1:(ncol(out)-2))]
  return(out)
}
