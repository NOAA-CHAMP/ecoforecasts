function plot_data_and_ranges(stn, fld, rngs)
%function plot_data_and_ranges(stn, fld, rngs)
%
% Plot raw station data from the field named "fld" in the station data struct
% "stn". Annotate the plot with the value ranges specified in struct "rngs".
% Struct "rngs" is OPTIONAL, but if specified, it must contain at least two
% fields, "ranges" and "rangeNames". It is normally one of the result structs
% returned by a call to, e.g., the M-Function "seasonal_data_ranges" (q.v.).
%
% Also plots histogram of raw data values, broken up by the values in rngs.
%
% Last Saved Time-stamp: <Thu 2008-09-11 12:32:28 Eastern Daylight Time gramer>

  if ( ~isstruct(stn) || ~isfield(stn, fld) )
    error('First arg must be a struct, second arg a field name in struct.');
  end;
  if ( ~exist('rngs', 'var') )
    rngs.ranges = [];
  elseif ( ~isstruct(rngs) || ~isfield(rngs, 'ranges') || ~isfield(rngs, 'rangeNames') )
    error(['Third arg "rngs", if given, must be a struct having a vector field ' ...
           ' of numbers "ranges", and a cell array of strings field "rangeNames".']);
  end;


  figure;
  plot(stn.(fld).date, stn.(fld).data);
  datetick;

  % Ignore the -Inf and +Inf on either end
  for idx = 2:(length(rngs.ranges)-1)
    annotline([], rngs.ranges(idx), rngs.rangeNames{idx}, 'red');
  end;


  figure;
  hist(stn.(fld).data, rngs.ranges(2:end-1));

return;
