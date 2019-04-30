function km = lat2km(lat1,lat2)
%function km = lat2km([lat1[,lat2]])
  lonin=-80;
  if ( exist('lat2','var') )
    km = distance_wgs84(lat1,lonin,lat2,lonin);
  elseif ( exist('lat1','var') )
    km = distance_wgs84(lat1,lonin,lat1+1,lonin);
  else
    km = distance_wgs84(25,lonin,1,lonin);
  end;
return;
