function stn = station_partition_periods(stn,fld,pers)
%function stn = station_partition_periods(stn,fld[,pers])
%
% Partition energy in time series STN.(FLD) into per-hour passbands given by
% the Nx2 vector PERS (DEFAULT gives 3-11h, 11-16h, 17-28h, 2-10d).
%
% Returns STN with new TS fields, e.g., STN.([FLD,'_3_h_hp_8_h_lp']).
%
% Last Saved Time-stamp: <Mon 2018-02-26 22:09:49 Eastern Standard Time gramer>

  if ( ~exist('pers','var') || isempty(pers) )
    pers = [ 3,11 ; 11,16 ; 17,28 ; 2*24,10*24 ];
  end;

  for perix = 1:size(pers,1)
    perfld = sprintf('%s_%g_h_hp_%g_h_lp',fld,pers(perix,1),pers(perix,2));
    %DEBUG:    disp(perfld);
    try,
      stn = verify_variable(stn,perfld);
      stn.(perfld) = ts_nanify_gaps(stn.(perfld),pers(end)/24);
    catch,
      %catchwarn;
    end;
  end;

return;
