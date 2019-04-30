function trueOrFalse = is_ts(ts)
%function trueOrFalse = is_ts(ts)
%
% A time series TS is a struct with two fields, TS.date and TS.data, having
% equal numbers of elements, and with TS.date numeric.
% 
% Last Saved Time-stamp: <Sat 2011-04-16 14:03:26  Lew.Gramer>

  trueOrFalse = ( isfield(ts,'date') && isfield(ts,'data') && ...
                  isnumeric(ts.date) && (numel(ts.date) == numel(ts.data)) );

return;
