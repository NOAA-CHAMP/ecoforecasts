function [newdts,newts,newprof] = window_func_prof(dts,ts,prof,varargin)
%function [newdts,newts,newprof] = window_func_prof(dts,ts,prof,func,nhrs,samplehrs,wintyp,plotResults)
%
% Call WINDOW_FUNC repeatedly on each COLUMN of time series profile PROF.
% Similarly, call WINDOW_FUNC on DTS and (vector) TS.
%
% Last Saved Time-stamp: <Mon 2018-02-26 12:49:04 EST lew.gramer>

  [newdts,newts] = window_func(dts,ts,varargin{:});
  newprof = repmat(nan,[numel(newdts),size(prof,2)]);
  for ix = 1:size(prof,2);
    if ( any(~isnan(prof(:,ix))) )
      [newdts_prof,newprof_prof] = window_func(dts,prof(:,ix),varargin{:});
      [dtix,pdtix] = intersect_dates(newdts,newdts_prof);
      newprof(dtix,ix) = newprof_prof(pdtix,1);
      clear newdts_prof newprof_prof
    end;
  end;

return;
