1;
%% SCRIPT hc_range_field_subsample.m
%
% Subsample seafloor slope angle ANG and Horizontal Convection Range HCRNG
% from .MAT file BASEMATFNAME (DEFAULT: 'FRT_depth_and_beta_92m'), and save
% together with data from .MAT file SUBMATFNAME (DEFAULT: [BASEMATFNAME,...
% '_NANMEAN_3_3_4']) into new .MAT file [SUBMATFNAME,'_hc_range'].
%
% Last Saved Time-stamp: <Thu 2017-05-04 20:59:45 Eastern Daylight Time gramer>

if ( ~exist('basematfname','var') || isempty(basematfname) )
  basematfname = 'FRT_depth_and_beta_92m';
end;
if ( ~exist('submatfname','var') || isempty(submatfname) )
  submatfname = [basematfname,'_NANMEAN_3_3_4'];
end;

disp(['Loading ',basematfname,'.mat']);
fullbath = load([basematfname,'_hc_range.mat']);

disp(['Loading ',submatfname,'.mat']);
load([submatfname,'.mat']);

[LON,LAT] = meshgrid(lon,lat);

ang = interp_field(fullbath.lat,fullbath.lon,fullbath.ang,LAT(:),LON(:),'nearest');
ang = reshape(ang,size(LAT)); ang = ang';

disp('Interpolating HCRNG using {@NANMAX,3,3,4}');
hcrng = interp_field(fullbath.lat,fullbath.lon,fullbath.hcrng,LAT(:),LON(:),{@nanmax,3,3,4});
hcrng = reshape(hcrng,size(LAT)); hcrng = hcrng';

disp(['Saving ',[submatfname,'_hc_range.mat']]);
save([submatfname,'_hc_range.mat'],'lon','lat','h','bet','ang','hcrng');
