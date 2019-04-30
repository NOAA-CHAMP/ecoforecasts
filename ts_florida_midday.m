function idx = ts_florida_midday(ts,dayhours)
%function idx = ts_florida_midday(ts,dayhours)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during (or seasonally close to) midday UTC hours
% for Florida: ISMEMBER(GET_HOUR(TS.date),DAYHOURS) (DEFAULT: [15:18]).

  if ( ~exist('dayhours','var') || isempty(dayhours) )
    dayhours = [15:18];
  end;
  idx = find(ismember(get_hour(ts.date),dayhours));

return;
