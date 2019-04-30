function newts = subset_ts(ts,idx)
%function newts = subset_ts(ts,idx)
%
% Return a new time series struct NEWTS containing only those indices of time
% series TS selected by IDX, a logical or index vector, or a function_handle
% for a function accepting a time series and returning a vector of indices.
%
% As of 2017, now also subsets .prof and .rawprof (for, e.g., ADCP profiles),
% and .field (for, e.g., 2-D SST fields) based on IDX. However, verifies
% first that SIZE(...,1) of each of these fields matches NUMEL(TS.date).
%
% Last Saved Time-stamp: <Mon 2017-03-13 14:46:35 Eastern Daylight Time gramer>

  if ( isa(idx,'function_handle') )
    idx = idx(ts);
  elseif ( (~isnumeric(idx) && ~islogical(idx)) || ~isvector(idx) )
    error('Second arg IDX must be a function handle, numeric- or logical-vector!');
  end;

  newts.date = ts.date(idx);
  % Some time series have ONLY .prof or .field, etc.
  if ( isfield(ts,'data') )
    newts.data = ts.data(idx);
  end;
  if ( isfield(ts,'prof') )
    if ( ndims(ts.prof) == 2 && size(ts.prof,1) == numel(ts.date) )
      newts.prof = ts.prof(idx,:);
    end;
  end;
  if ( isfield(ts,'rawprof') )
    if ( ndims(ts.rawprof) == 2 && size(ts.rawprof,1) == numel(ts.date) )
      newts.rawprof = ts.rawprof(idx,:);
    end;
  end;
  if ( isfield(ts,'field') )
    if ( ndims(ts.field) == 3 && size(ts.field,1) == numel(ts.date) )
      newts.field = ts.field(idx,:,:);
    end;
  end;

return;
