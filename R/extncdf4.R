#' Extraction Function
#'
#' This function loads the ncdf4 file of NCEP as a matrix.
#'
#' @param filename Name of the input ncdf4 file (e.g. 'uwnd.mon.mean.nc')
#' @param level Barometric pressure to be output (e.g. 925hPa -> 925)
#' @param dsy Start year of original data (e.g. Jan,1948 -> 1948)
#' @param year Range to be output (e.g. 1985 to 2020 -> 1985:2020)
#' @return A matrix of the infile
#' @examples Uwind = extncdf4('uwnd.mon.mean.nc',850,1948,1985:2020)
#' @export
extncdf4 = function(filename,level,dsy,year){
  ncin = nc_open(paste(getwd(),"/",filename,sep=""))
  value = ncvar_get(ncin,strsplit(filename,"[.]")[[1]][1])
  return(value[,,which(ncvar_get(ncin,"level") == level),((year[1]-dsy)*12+1):((rev(year)[1]-dsy)*12+12)])
}
