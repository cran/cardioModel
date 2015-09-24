rename_column <- function(x, old_name, new_name, required){
  is_present <- old_name %in% colnames(x)
  if (required & !is_present){
    stop(paste("Expected column", old_name, "not found!"))
  }
  if (is_present){
    colnames(x)[colnames(x) == old_name] <- new_name
  }
  x
}

myPremodeling <- function(myData){
  myData$fNTAFD <- as.factor(myData$NTAFD)
  myData$ID <- as.factor(myData$ID)
  myData$RESPONSE <- as.numeric(myData$RESPONSE)
  myData$EXPOSURE <- as.numeric(myData$EXPOSURE)
  if("PERIOD" %in% colnames(myData)){
    myData$PERIOD <- as.factor(myData$PERIOD)
  }
  return(myData)
}
