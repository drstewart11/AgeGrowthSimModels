#' Assigns individual length data to a length bin
#'
#' @param len = length data
#' @param by = sets the width of the length bin
#' @param dir = round "up" or round "down"
#' @return length bin
#' @export
lengthbins<-function(len,by,dir=c("up","down")){
  if(dir=="up"){
    by=abs(by)
    by*(len%/%by+as.logical(len%%by))
  }else if(dir=="down"){
    by=-1*abs(by)
    by*(len%/%by+as.logical(len%%by))
  }
}
