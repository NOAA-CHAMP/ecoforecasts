function new_data = gap_fill(dates, data, method)
%function new_data = gap_fill(dates, data, method)
%
% Replace all NaNs in 'data' with interpolated (e.g., spline-fit) values
%
% Last Saved Time-stamp: <Thu 2009-03-19 07:16:13 Eastern Daylight Time gramer>

  if ( ~exist('method', 'var') || isempty(method) )
    method = 'cubic';
  end

  idx = ~isnan(data);
  new_data = interp1(dates(idx), data(idx), dates, method);

return;
