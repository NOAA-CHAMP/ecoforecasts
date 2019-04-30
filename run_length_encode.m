function [rl,val] = run_length_encode(vals,dim)
%function [rl,val] = run_length_encode(vals,dim)
%
% Get length of contiguous indices having the same value. (Run-length
% encoding algorithm courtesy "MATLAB array manipulation tips and tricks",
% Peter Acklam, Norway.) VALS may be a string cell array or matrix, or
% numeric: if numeric, optional 2nd arg DIM will be passed to DIFF (v.)
%
% Last Saved Time-stamp: <Mon 2013-11-04 17:06:40 Eastern Standard Time gramer>

  warning('run_length_encode:generic',...
          'This code not yet well-tested - and use of DIM arg definitely still broken!');

  if ( ischar(vals) )
    if ( ~iscell(vals) )
      vals = cellstr(vals);
    end;
    kern = find(~strcmp(vals(1:end-1), vals(2:end)));
    dimlen = length(vals);
  else
    if ( ~exist('dim','var') || isempty(dim) )
      kern = find(diff(vals) ~= 0);
      dimlen = length(vals);
    else
      kern = find(diff(vals,[],dim) ~= 0);
      dimlen = size(vals,dim);
    end;
  end;

  rl = diff([ 0 kern(:)' dimlen ]);
  val = vals([ kern(:)' dimlen ]);

return;
