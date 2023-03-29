#' Distance Calculating Function
#'
#' This function calculates the distance of coordinates on the ellipsoid.
#'
#' @param data Name of the input data. The format of the data follows the output value of Extraction function ('extmonth' or 'extseason') in the same package.
#' @param degree Interval between measurements of NCEP data (e.g. 2.5 x 2.5 -> 2.5)
#' @return A matrix of the infile
#' @examples distance = partial(JJASON_Uwind,2.5)
#' @export
partial = function(data,degree){
  lonlev = levels(as.factor(data[,1])) ; latlev = levels(as.factor(data[,2]))
  
  dx_p1 = cbind(rep(0,length(latlev)),seq(latlev[length(latlev)],latlev[1],-degree))
  dx_p2 = cbind(rep(2.5,length(latlev)),seq(latlev[length(latlev)],latlev[1],-degree))
  dx = distVincentyEllipsoid(dx_p1,dx_p2,a=6378137, b=6356752.3142, f=1/298.257223563) ; rm(dx_p1,dx_p2)
  
  dy_imsi = as.numeric(c(latlev[length(latlev)],seq(latlev[length(latlev)],latlev[1],-degree),latlev[1]))
  n = length(dy_imsi)
  p1 = cbind(rep(0,(n-2)),dy_imsi[1:(n-2)]) ; p2 = cbind(rep(0,(n-2)),dy_imsi[3:n])
  dy=(1/2)*distVincentyEllipsoid(p1,p2,a=6378137, b=6356752.3142, f=1/298.257223563)
  
  return(cbind(dx,dy))
}

