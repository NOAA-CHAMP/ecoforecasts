function newts = date_range_ts(ts,dt1,dt2)
%function newts = date_range_ts(ts,dt1[,dt2])
%
% Return a new time series struct NEWTS containing only those indices of time
% series TS that fall between DATENUMs DT1 and DT2. See also TS_DATE_RANGE.
% If DT1 is a 2-vector and DT2 is missing, use date range [DT1(1),DT1(2)].
%
% Subsets .data, .prof and .rawprof (for, e.g., ADCP profiles), and .field
% (for, e.g., 2-D SST fields) based on DT1 and DT2.
%
% Last Saved Time-stamp: <Tue 2017-08-22 17:02:38 Eastern Daylight Time gramer>

  if ( ~exist('dt2','var') )
    if ( numel(dt1) ~= 2 || ~isnumeric(dt1) )
      error('If no second arg, first arg must be a 2-vector of DATENUM');
    else
      dt2 = dt1(2);
      dt1 = dt1(1);
    end;
  end;

  idx = ts_date_range(ts,dt1,dt2);
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
