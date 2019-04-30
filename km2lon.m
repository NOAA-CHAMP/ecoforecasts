function deg=km2lon(km1,km2,latin,az)
%function deg=km2lon([km1[,km2[,latin[,az]]]])
  if (~exist('lonin','var')||isempty(lonin)); lonin=-80; end;
  if (~exist('latin','var')||isempty(latin)); latin=+25; end;
  if (~exist('az','var')||isempty(az)); az=90; end;

  if exist('km2','var') && ~isempty(km2)
    [latout,lonout] = reckon_wgs84(latin,lonin,km2-km1,az);
  elseif exist('km1','var') && ~isempty(km1)
    [latout,lonout] = reckon_wgs84(latin,lonin,km1,az);
  else
    [latout,lonout] = reckon_wgs84(latin,lonin,1,az);
  end;
  deg=abs(lonout-lonin);
return;
