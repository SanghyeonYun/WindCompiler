#' Extraction Function (Seasonal)
#'
#' This function works based on the data extracted by the extncdf4 function.
#' Extracts a seasonal value for a specific year.
#'
#' @param data Name of the input data. The format of the data follows the output value of the 'extncdf4' function in the same package.
#' @param year Years to be extracted
#' @param degree Interval between measurements of NCEP data (e.g. 2.5 x 2.5 -> 2.5)
#' @param start Start of season to be extracted (e.g. June to November -> 6)
#' @param end End of season to be extracted (e.g. June to November -> 11)
#' @return A matrix of the infile
#' @examples JJASON_Uwind = extmonth(Uwind,1985:2020,2.5,6,11)
#' @export
extseason = function(data,year,degree,start,end){
  beforemonth = start-1 ; monthlyrange = length(start:end)
  vec = as.vector(data) ; EndMon = seq(beforemonth,(length(year)*12),12) ; imsi = NULL
  for(i in EndMon){imsi = c(imsi,vec[(10512*i+1):(10512*(i+monthlyrange))])}
  imsimat = matrix(imsi,nrow=(dim(data)[1]*dim(data)[2]))
  extraction.imsi = matrix(rep(NA,length(year)*dim(data)[1]*dim(data)[2]),ncol=length(year))
  for(j in 1:(dim(data)[1]*dim(data)[2])){for(k in 1:length(year)){extraction.imsi[j,k] = mean(imsimat[j,monthlyrange*(k-1)+c(1:monthlyrange)])}}
  if(degree == 2.5){extraction.imsi = data.frame(rep(seq(0,357.5,2.5),73),rep(seq(90,-90,-2.5),each=144),extraction.imsi)
  names(extraction.imsi) = c('lon','lat',year)}
  if(degree == 0.25){extraction.imsi = data.frame(rep(seq(0,357.5,0.25),1431),rep(seq(90,-90,-0.25),each=721),extraction.imsi)
  names(extraction.imsi) = c('lon','lat',year)}
  return(extraction.imsi)
}