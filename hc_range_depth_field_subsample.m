1;
%% SCRIPT hc_range_depth_field_subsample.m
%
% Subsample seafloor depth (H), slope (BET), slope angle (ANG, if available),
% Horizontal Convection Range (HCRNG), and Horizontal Convection Depth (HCH)
% from HC_RANGE_DEPTH .MAT file named based on BASEMATFNAME (DEFAULT: 'FRT_depth_and_beta_92m'),
% and combine with data from bathymetry subsample .MAT file SUBMATFNAME (DEFAULT:
% [BASEMATFNAME,... '_NANMEAN_3_3_4']) into new .MAT file [SUBMATFNAME,'_hc_range'].
%
% Last Saved Time-stamp: <Sat 2018-11-17 15:52:26 Eastern Standard Time gramer>

if ( ~exist('basematfname','var') || isempty(basematfname) )
  basematfname = 'FRT_depth_and_beta_92m';
end;
if ( ~exist('interpmethod','var') )
  interpmethod = {@nanmean,3,3,4};
end;
if ( ~exist('submatfname','var') || isempty(submatfname) )
  submatfname = [basematfname,'_',interpmethod_to_text(interpmethod)];
end;

disp(['Loading ',basematfname,'.mat']);
fullbath = load([basematfname,'_hc_range_depth.mat']);

disp(['Loading ',submatfname,'.mat']);
load([submatfname,'.mat']);

[LON,LAT] = meshgrid(lon,lat);

if ( exist('ang','var') )
  mth = 'nearest';
  disp(['Interpolating ANG using INTERPMETHOD ',interpmethod_to_text(mth)]);
  ang = interp_field(fullbath.lat,fullbath.lon,fullbath.ang,LAT(:),LON(:),mth);
  ang = reshape(ang,size(LAT)); ang = ang';
end;

mth = interpmethod; mth{1} = @nanmax;
disp(['Interpolating HCRNG using INTERPMETHOD ',interpmethod_to_text(mth)]);
hcrng = interp_field(fullbath.lat,fullbath.lon,fullbath.hcrng,LAT(:),LON(:),mth);
hcrng = reshape(hcrng,size(LAT)); hcrng = hcrng';

mth = interpmethod; mth{1} = @nanmin;
disp(['Interpolating HCH using INTERPMETHOD ',interpmethod_to_text(mth)]);
hch = interp_field(fullbath.lat,fullbath.lon,fullbath.hch,LAT(:),LON(:),mth);
hch = reshape(hch,size(LAT)); hch = hch';

resmatfname = [submatfname,'_hc_range_depth.mat'];
disp(['Saving ',resmatfname]);
if ( exist('ang','var') )
  save(resmatfname,'lon','lat','h','bet','ang','hcrng','hch');
else
  save(resmatfname,'lon','lat','h','bet','hcrng','hch');
end;
