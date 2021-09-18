
# check and categorize location info based on patterns
catPattern <- function(locIn){
  
  out <- rep(NA, length = length(locIn))
  
  out[grep(
    "[P|p]rov|[P|p]ro.", 
    locIn, ignore.case = TRUE)] <- "province"
  
  out[grep(
    "[C|c]ity|[M|m]unici", 
    locIn, 
    ignore.case = TRUE)] <- "municipality"
  
  out[grep(
    "[B|b]rg|[B|b]arangay|[B|b]arrio", 
    locIn, ignore.case = TRUE)] <- "barangay"
  
  return(out)
  
}


# summarize final admin category of each entry, MULTIPLE if more than one match
sumLocTable <- function(tableIn){
 
  out <- apply(tableIn[ , -1], 1, function(x){
    
    if(sum(is.na(x)) == length(x)){
      
      return(NA)
    }else if(length(unique(x[!is.na(x)])) == 1){
      
      return(unique(x[!is.na(x)]))
      
      
    }else{
      
      return("MULTIPLE")
    }
    
  })
  
  return(out)
  
  
}
