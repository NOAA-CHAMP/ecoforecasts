function deg=km2lat(km1,km2,az)
%function deg=km2lat([km1[,km2[,az]]])
  latin=25; lonin=-80;
  if (~exist('az','var')||isempty(az)); az=0; end;
  if exist('km2','var')
    [latout,lonout] = reckon_wgs84(latin,lonin,km2-km1,az);
  elseif exist('km1','var')
    [latout,lonout] = reckon_wgs84(latin,lonin,km1,az);
  else
    [latout,lonout] = reckon_wgs84(latin,lonin,1,az);
  end;
  deg=abs(latout-latin);
return;
