function [X,Y,zn] = latlon2utm(LAT,LON,zn)
%function [X,Y,zn] = latlon2utm(LAT,LON,zn)
%
% Using MAPPING TOOLBOX to convert arrays of latitude/longitude to UTM
% easting/northing UTM zone ZN. If ZN is not specified or empty, and function
% UTMZONE (Map Toolbox) exists, then ZN=UTMZONE(LAT(1),LON(1)). Otherwise,
% choose DEFAULT ZN="17R" which corresponds to large parts of Florida.
% 
% Last Saved Time-stamp: <Wed 2019-02-06 12:11:50 EST lew.gramer>

  persistent ellipsoid
  persistent utmstruct

  if ( ~exist('zn','var') || isempty(zn) )
    if ( exist('utmzone') > 1 )
      zn = utmzone(LAT(1),LON(1));
    else
      zn = '17R';
    end;
  end;

  if ( ~isfield(utmstruct,'zone') || ~strcmpi(utmstruct.zone,zn) )
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

  % if ( ~isvector(LAT) && ~isvector(LON) )
  %   [X,Y] = mfwdtran(utmstruct,LAT,LON);
  % else
  %   % MFWDTRAN relies on MESHGRID format
  %   warning('Ecoforecasts:UTMgrid','Input were vectors, output are matrices');
  %   [LAT,LON] = meshgrid(LAT,LON);
  %   [X,Y] = mfwdtran(utmstruct,LAT,LON);
  % end;

  % MFWDTRAN USED TO rely on MESHGRID format: Now it apparently doesn't??
  [X,Y] = mfwdtran(utmstruct,LAT,LON);


return;
