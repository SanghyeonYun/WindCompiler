#' Wind Component Calculating Function
#'
#' This function calculates the Component of Wind.
#' Winds can be decomposed into rotational and divergent components. These are outputted as vorticity and divergence, respectively.
#'
#' @param uwnd Name of the meridional wind. The format of the data follows the output value of Extraction function ('extmonth' or 'extseason') in the same package.
#' @param vwnd Name of the zonal wind. The format of the data follows the output value of Extraction function ('extmonth' or 'extseason') in the same package.
#' @param degree Interval between measurements of NCEP data (e.g. 2.5 x 2.5 -> 2.5)
#' @param year Range to be analysis (e.g. 1985 to 2020 -> 1985:2020)
#' @param component The wind component to be calculated. This option has a choice of 'rotational' and 'divergent'.
#' @return A matrix of the infile
#' @examples vorticity = windcomp(JJASON_Uwind,2.5,1985:2020,'rotational')
#' @export
windcomp = function(uwnd,vwnd,degree,year,component){
  lonlev = levels(as.factor(uwnd[,1])) ; latlev = levels(as.factor(uwnd[,2]))
  
  dx_p1 = cbind(rep(0,length(latlev)),seq(latlev[length(latlev)],latlev[1],-degree))
  dx_p2 = cbind(rep(2.5,length(latlev)),seq(latlev[length(latlev)],latlev[1],-degree))
  dx = distVincentyEllipsoid(dx_p1,dx_p2,a=6378137, b=6356752.3142, f=1/298.257223563) ; rm(dx_p1,dx_p2)
  
  dy_imsi = as.numeric(c(latlev[length(latlev)],seq(latlev[length(latlev)],latlev[1],-degree),latlev[1]))
  n = length(dy_imsi)
  p1 = cbind(rep(0,(n-2)),dy_imsi[1:(n-2)]) ; p2 = cbind(rep(0,(n-2)),dy_imsi[3:n])
  dy=(1/2)*distVincentyEllipsoid(p1,p2,a=6378137, b=6356752.3142, f=1/298.257223563)
  
  U = as.matrix(uwnd[c(order(uwnd[,1])),][,3:(length(year)+2)])
  V = as.matrix(vwnd[c(order(vwnd[,1])),][,3:(length(year)+2)])
  Lon = as.numeric(levels(as.factor(uwnd[,1]))) ; ilen = length(Lon)
  Lat = rev(as.numeric(levels(as.factor(uwnd[,2])))) ; jlen = length(Lat)
  
  if(component == 'rotational')
  {
    dVdx = NULL ; dUdy = NULL
    for(k in 1:length(year))
    {
      dVdx_ann=NULL ; dUdy_ann=NULL
      U_ann=matrix(U[,k],ncol=ilen) ; V_ann=matrix(V[,k],ncol=ilen)
      for(i in 1:ilen) {d=diff(c(U_ann[1,i], U_ann[,i], U_ann[jlen,i])) ; dUdy_ann=c(dUdy_ann,(d[2:(jlen+1)]+d[1:jlen])/2*(1/dy))}
      for(j in 1:jlen) {d=diff(c(V_ann[j,ilen], V_ann[j,], V_ann[j,1])) ; dVdx_ann=c(dVdx_ann,(d[2:(ilen+1)]+d[1:ilen])/2*(1/dx[j]))}
      dVdx=cbind(dVdx,dVdx_ann) ; dUdy=cbind(dUdy,dUdy_ann)
    }
    dVdx = cbind(rep(Lon,length(Lat)),rep(Lat,each=length(Lon)),dVdx)
    dUdy = cbind(rep(Lon,each=length(Lat)),rep(Lat,length(Lon)),dUdy)
    dUdy = dUdy[c(order(-dUdy[,2])),] ; dUdy[,3:(length(year)+2)] = -1*dUdy[,3:(length(year)+2)]
    vorticity = dVdx[,3:(length(year)+2)] - dUdy[,3:(length(year)+2)]
    
    vorticity = cbind(rep(seq(Lon[1],Lon[length(Lon)],degree),length(Lat)),
                      rep(seq(Lat[1],Lat[length(Lat)],-degree),each=length(Lon)),
                      vorticity)
    colnames(vorticity) = c('lon','lat',year)
    return(vorticity)
  }
  
  if(component == 'divergent')
  {
    dUdx = NULL ; dVdy = NULL
    for(k in 1:length(year))
    {
      dUdx_ann=NULL ; dVdy_ann=NULL
      U_ann=matrix(U[,k],ncol=ilen) ; V_ann=matrix(V[,k],ncol=ilen)
      for(j in 1:jlen) {d=diff(c(U_ann[j,ilen], U_ann[j,], U_ann[j,1])) ; dUdx_ann=c(dUdx_ann,(d[2:(ilen+1)]+d[1:ilen])/2*(1/dx[j]))}
      for(i in 1:ilen) {d=diff(c(V_ann[1,i], V_ann[,i], V_ann[jlen,i])) ; dVdy_ann=c(dVdy_ann,(d[2:(jlen+1)]+d[1:jlen])/2*(1/dy))}
      dUdx=cbind(dUdx,dUdx_ann) ; dVdy=cbind(dVdy,dVdy_ann)
    }
    
    dUdx = cbind(rep(Lon,length(Lat)),rep(Lat,each=length(Lon)),dUdx)
    dVdy = cbind(rep(Lon,each=length(Lat)),rep(Lat,length(Lon)),dVdy)
    dVdy = dVdy[c(order(-dVdy[,2])),] ; dVdy[,3:(length(year)+2)] = -1*dVdy[,3:(length(year)+2)]
    divergence = dUdx[,3:(length(year)+2)] + dVdy[,3:(length(year)+2)]
    
    divergence = cbind(rep(seq(Lon[1],Lon[length(Lon)],degree),length(Lat)),
                       rep(seq(Lat[1],Lat[length(Lat)],-degree),each=length(Lon)),
                       divergence)
    colnames(divergence) = c('lon','lat',year)
    return(divergence)
  }
}

