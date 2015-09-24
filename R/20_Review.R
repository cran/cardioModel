mySummary <- function(results){

  # rbind all dataframes
  results <- ldply(results, data.frame)

  # eliminate duplicates for etas-slope which repeat accross Eta-EC50/EMAX
  results <- unique(results)

  # give all models tested order from the best (lower AIC) to worse (higher AIC)
  out <- results[with(results, order(AIC)),]

  # remove first part of ERROR.MESSAGE so that function names don't show
  out$ERROR.MESSAGE <- gsub(".* : \\n  ", "", out$ERROR.MESSAGE)
  out$ERROR.MESSAGE <- gsub("\\n$", "", out$ERROR.MESSAGE)
  return(out)
}
