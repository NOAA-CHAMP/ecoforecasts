function [LAT,LON] = utm2latlon(X,Y,zn)
%function [LAT,LON] = utm2latlon(X,Y,zn)
%
% Using MAPPING TOOLBOX to convert arrays of UTM easting/northing X and Y for
% UTM zone ZN (DEFAULT: "17R" for areas of Florida), into latitude/longitudes.
% 
% Other zones of interest for coral reefs: 55P Saipan/CNMI, 2L Anerican Samoa,
% 19Q for western Puerto Rico and Dominican Republic, and 20Q for USVI.
% 
% Last Saved Time-stamp: <Tue 2019-02-05 16:14:03 EST lew.gramer>

  persistent ellipsoid
  persistent utmstruct

  if ( ~exist('zn','var') || isempty(zn) )
    zn = '17R';
    warning('Using default UTM zone %s!',zn);
  end;

  if ( ~isfield(utmstruct,'zone') || ~strcmpi(utmstruct.zone,zn) )
    utmstruct=[]; clear utmstruct;
    % Obtain the suggested ellipsoid vector and name for this zone
    [ellipsoid,estr] = utmgeoid(zn);
    % Set up the UTM coordinate system based on this information
    utmstruct = defaultm('utm'); 
    utmstruct.zone = zn;
    utmstruct.geoid = ellipsoid; 
    utmstruct = defaultm(utmstruct);
  else
    %DEBUG:    disp(['Zone is ',utmstruct.zone]);
  end;

  % %[X,Y] = mfwdtran(utmstruct,latlim,lonlim)
  % if ( ~isvector(X) && ~isvector(Y) )
  %   %[LAT,LON] = utm2deg(949622.38,589331.52,'17 R')
  %   %[LAT,LON] = minvtran(utmstruct,949622.38,589331.52)
  %   [LAT,LON] = minvtran(utmstruct,X,Y);
  % else
  %   % MINVTRAN relies on MESHGRID format
  %   warning('Ecoforecasts:UTMgrid','Input were vectors, output are matrices');
  %   [X,Y] = meshgrid(X,Y);
  %   [LAT,LON] = minvtran(utmstruct,X,Y);
  % end;

  % MINVTRAN USED TO rely on MESHGRID format: Now it apparently doesn't??
  [LAT,LON] = minvtran(utmstruct,X,Y);

return;
