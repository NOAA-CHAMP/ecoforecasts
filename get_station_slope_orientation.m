function stn = get_station_slope_orientation(stn,N,method,rad,useHighestRes,bathfld)
%function stn = get_station_slope_orientation(stn[,N[,method[,rad[,useHighestRes[,bathfld]]]]])
%
% Return STN with fields .slope, .slope_orientation, .isobath_orientation.
% Calls READ_HIRES_BATHYMETRY if needed, then FIND_NGDC_SLOPE. Optional args
% N and METHOD are passed to the latter, RAD and USEHIGHESTRES to the former.
%
% Last Saved Time-stamp: <Thu 2016-08-18 00:02:36 Eastern Daylight Time gramer>

  if ( ~exist('N','var') )
    N = [];
  end;
  if ( ~exist('method','var') )
    method = [];
  end;
  if ( ~exist('rad','var') )
    rad = [];
  end;
  if ( ~exist('useHighestRes','var') )
    useHighestRes = [];
  end;

  bathfld = 'ngdc_hires_bathy';

  if ( ~isfield(stn,bathfld) )
    stn = read_hires_bathymetry(stn,rad,[],useHighestRes);
  end;
  [stn.slope,stn.slope_orientation,stn.isobath_orientation,stn.(bathfld)] = ...
      find_ngdc_slope(stn.(bathfld),stn.lon,stn.lat,N,method)

return;
