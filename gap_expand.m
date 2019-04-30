function [new_indices, new_data] = gap_expand(indices, data, dt)
%function [new_indices, new_data] = gap_expand(indices, data, dt)
%
% Expand a time series 'data' (e.g., sea temperatures), presumed to be indexed
% by another time series 'indices' (e.g., measurement dates), so that it has
% no indicial gaps: any gaps from the original data are filled with NaNs. The
% 'new_indices' contains a continuous series of equally-spaced datenums from
% the beginning to the end of the vector 'indices' - including datenums for
% all times where a NaN was added to 'new_data'.
%
% Optional DT is time series resolution (DEFAULT: MEDIAN(DIFF(INDICES))).
%
% Last Saved Time-stamp: <Tue 2017-06-06 16:30:42 Eastern Daylight Time lew.gramer>

  % Data may contain nans, but dates ('indices') may not!
  if ( any(~isfinite(indices)) )
    error('Indices must not contain NaN or Inf! Please review data...');
  end;
  if ( min(diff(indices)) < 0 )
      error('Indices must increase monotonically! Please review data...');
  end;
  % Time resolution of time series
  if ( ~exist('dt','var') || isempty(dt) )
    dt = median(diff(indices));
  end;


  % Always deal with column vectors
  indices = indices(:);
  data = data(:);

  % NOTE: 'indices' may come in as a float vector beginning with any value
  idx = round((indices - indices(1)) ./ dt) + 1;

  % Create a new TS big enough to hold all the filled gaps
  new_n = floor(idx(end) - idx(1) + 1);
  % Using LINSPACE will be problematic for time series whose reporting time
  % changes during their record - e.g., hourly data that are reported at
  % 0:56:00 for part of the record, then 0:00:00 for the rest (e.g., LPPR1).
  new_indices = linspace(indices(1), indices(end), new_n)';

  new_data = repmat(nan, [new_n size(data,2)]);
  new_data(idx) = data;

return;
