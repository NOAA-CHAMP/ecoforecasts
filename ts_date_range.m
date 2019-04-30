function idx = ts_date_range(ts,date1,date2)
%function idx = ts_date_range(ts,date1,date2)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full between DATENUMs DATE1 and DATE2, excluding those
% at exactly DATE2. With just two args, DATE1 must be a 2-vector of DATENUM.

  if ( numel(date1) == 2 )
    date2 = date1(2);
    date1 = date1(1);
  end;
  idx = find(date1 <= ts.date & ts.date < date2);

return;
