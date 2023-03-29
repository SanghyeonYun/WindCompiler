#' Extraction Function (Month)
#'
#' This function works based on the data extracted by the extncdf4 function.
#' Extracts a specific month value for a specific year.
#'
#' @param data Name of the input data. The format of the data follows the output value of the 'extncdf4' function in the same package.
#' @param year Years to be extracted
#' @param degree Interval between measurements of NCEP data (e.g. 2.5 x 2.5 -> 2.5)
#' @param month Month to be extracted (e.g. Jan -> 1, Aug -> 8)
#' @return A matrix of the infile
#' @examples Jan_Uwind = extmonth(Uwind,1985:2020,2.5,1)
#' @export
extmonth = function(data,year,degree,month){
  imsi = matrix(rep(NA,length(year)*dim(data)[1]*dim(data)[2]),ncol=(length(year)))
  for(i in 1:length(year)){imsi[,i] = data[,,month+(12*(i-1))]}
  if(degree == 2.5){imsi = data.frame(rep(seq(0,357.5,2.5),73),rep(seq(90,-90,-2.5),each=144),imsi)
  names(imsi) = c('lon','lat',year)}
  if(degree == 0.25){imsi = data.frame(rep(seq(0,357.5,0.25),1431),rep(seq(90,-90,-0.25),each=721),imsi)
  names(imsi) = c('lon','lat',year)}
  return(imsi)
}