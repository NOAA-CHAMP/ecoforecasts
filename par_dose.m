function dose = par_dose(par,nhrs)
%function dose = par_dose(par,nhrs)
%
% Convert hourly mean PAR time series (STRUCT with .date and .data expressed
% in [micro-mol/m^2/s]) into PAR "dose" in [mol/m^2/N*hrs]. Result DOSE is an
% *HOURLY*  time series (STRUCT with .date and .data) of cumulative dose for
% the past NHRS hours (DEFAULT: 24 hours for daily dose). Note PAR may also
% be a time series of UV, insolation, or other field with units [W/m^2]: in
% this case, resulting time series DOSE is in MJ/m^2/day.
%
% Last Saved Time-stamp: <Wed 2011-07-20 12:45:09  Lew.Gramer>

  if ( ~exist('nhrs','var') || isempty(nhrs) )
    nhrs = 24;
  end;
  if ( ~isnumeric(nhrs) || ~isscalar(nhrs) || floor(nhrs) ~= nhrs )
    error('Second arg must be an integral scalar number of HOURS for dose');
  end;

  sumfld = ['p_' num2str(nhrs) '_hour_sum'];

  x.p = par;
  x = verify_variable(x,sumfld);

  dose.date = x.(sumfld).date;
  % Convert 24*[micro-mol/m^2/s] -> [mol/m^2/day]
  % (Also converts 24*[W/m^2] -> [MJ/m^2/day])
  dose.data = x.(sumfld).data.*3600./1e6;
  x=[]; clear x;

return;
