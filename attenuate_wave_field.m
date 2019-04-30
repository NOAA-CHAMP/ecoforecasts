function newhs = attenuate_wave_field(lon,lat,bat,hs,dp)
%function newhs = attenuate_wave_field(lon,lat,bat,hs,dp)
%
% Continuously attenuate wave heights (3-D time series field HS with dominant
% direction 3-D field DP) using bathymetry field BAT. NOTE: LON,LAT,BAT,HS
% must describe the same coordinate field, that is, we need SIZE(BAT,1) ==
% SIZE(HS,2), SIZE(BAT,2) == SIZE(HS,3), ALL(SIZE(HS) == SIZE(DP)),
% LENGTH(LAT) == SIZE(BAT,1), and LENGTH(LON) == SIZE(BAT,2). Attenuated wave
% height time series field NEWHS (ALL(SIZE(NEWHS)==SIZE(HS))) is returned.
%
% Last Saved Time-stamp: <Tue 2016-05-24 10:38:53 Eastern Daylight Time gramer>

  ctrlon = median(lon(:));
  ctrlat = median(lat(:));

  [ig,ctrlonix] = min(abs(lon-ctrlon));
  [ig,ctrlatix] = min(abs(lat-ctrlat));

  [ig,minlatix] = min(lat);
  [ig,maxlatix] = max(lat);
  [ig,minlonix] = min(lon);
  [ig,maxlonix] = max(lon);

  for rix = 1:size(bat,1)
    for cix = 1:size(bat,2)
      minz_ne(rix,cix) = nanmin(nanmin(bat(rix:end,cix:end)));
      minz_nw(rix,cix) = nanmin(nanmin(bat(rix:end,1:cix)));
      minz_se(rix,cix) = nanmin(nanmin(bat(1:rix,cix:end)));
      minz_sw(rix,cix) = nanmin(nanmin(bat(1:rix,1:cix)));

      newhs(
    end;
  end;

%{
  if ( minlatix < maxlatix )
    Nixes = [0:numel(lat)-1];

    N_latixen = minlatix:ctrlatix;
    S_latixen = ctrlatix:maxlatix;
  else
    N_latixen = ctrlatix:maxlatix;
    S_latixen = minlatix:ctrlatix;
  end;
  if ( minlonix < ctrlonix )
    W_lonixen = minlonix:ctrlonix;
    E_lonixen = ctrlonix:maxlonix;
  else
    W_lonixen = ctrlonix:maxlonix;
    E_lonixen = minlonix:ctrlonix;
  end;
%}

return;
