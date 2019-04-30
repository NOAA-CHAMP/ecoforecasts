function result = ts_op(ts1,ts2,op,ix1,ix2,doPlot,tol)
%function result = ts_op(ts1,ts2,op,ix1,ix2,doPlot,tol)
%
% Apply OP (a function-handle, name, or punctuation operator) to coincident
% values from two time series. Calls INTERSECT_DATES(TS1.date,TS2.date) to
% find coincident data points. If optional IX1 or IX2 are *non-empty*, only
% use corresponding points from those sub-indices of TS1 or TS2, resp. IX1
% and IX2 can also be time series-indexing function_handles, e.g., TS_JFM.
% If optional arg DOPLOT is True, plot the result on a new figure. Either
% TS1 or TS2 may also be scalar, for distributive scalar-vector operations.
% Optional argument TOL is passed through to INTERSECT_DATES (v.).
%
% NOTE: MATRIX operators are treated as ARRAY operators, e.g., specifying an
% OP (third argument) of '*' is the same as specifying '.*'. If the user
% really wants a matrix operator on time series for some reason, they must
% specify it by function name/handle, e.g., >> TS3 = TS_OP(TS1,TS2,'MTIMES');
%
% As of 2017, now also operators on .prof and .rawprof (for, e.g., ADCP
% profiles), and .field (for, e.g., 2-D SST fields). However, verifies
% first that SIZE(...,1) of each of these fields matches NUMEL(TS.date).
%
% Last Saved Time-stamp: <Thu 2017-07-06 23:08:44 Eastern Daylight Time gramer>

  if ( isscalar(ts1) && isnumeric(ts1) && is_valid_ts(ts2) )
    val = ts1;
    ts1 = ts2;
    if isfield(ts1,'data');	ts1.data(:) = val;	end;
    if isfield(ts1,'prof');	ts1.prof(:) = val;	end;
    if isfield(ts1,'rawprof');	ts1.rawprof(:) = val;	end;
    if isfield(ts1,'field');	ts1.field(:) = val;	end;
  elseif ( isscalar(ts2) && isnumeric(ts2) && is_valid_ts(ts1) )
    val = ts2;
    ts2 = ts1;
    if isfield(ts2,'data');	ts2.data(:) = val;	end;
    if isfield(ts2,'prof');	ts2.prof(:) = val;	end;
    if isfield(ts2,'rawprof');	ts2.rawprof(:) = val;	end;
    if isfield(ts2,'field');	ts2.field(:) = val;	end;
  end;
  if ( ~isfield(ts1,'date') || ~isfield(ts2,'date') )
    error('At least one of TS1 and TS2 must be a time series STRUCT!');
  end;

  if ( ~exist('ix1','var') || isempty(ix1) )
    ix1 = 1:length(ts1.date);
  end;
  if ( ~exist('ix2','var') || isempty(ix2) )
    ix2 = 1:length(ts2.date);
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    tol = [];
  end;

  dodata = (isfield(ts1,'data') && isfield(ts2,'data'));
  doprof = (isfield(ts1,'prof') && isfield(ts2,'prof'));
  dorawprof = (isfield(ts1,'rawprof') && isfield(ts2,'rawprof'));
  dofield = (isfield(ts1,'field') && isfield(ts2,'field'));

  % Validate assumptions about time series profiles and/or 2-D fields
  if ( dodata )
    if ( size(ts1.data,1) ~= numel(ts1.date) || size(ts2.data,1) ~= numel(ts2.date) )
      error('Data field dimension 1 does not correspond to timestamps!');
    end;
  end;
  if ( doprof )
    if ( ndims(ts1.prof) ~= 2 || ndims(ts2.prof) ~= 2 ...
         || size(ts1.prof,1) ~= numel(ts1.date) || size(ts2.prof,1) ~= numel(ts2.date) )
      error('Profile dimension 1 does not correspond to timestamps!');
    end;
  end;
  if ( dorawprof )
    if ( ndims(ts1.rawprof) ~= 2 || ndims(ts2.rawprof) ~= 2 ...
         || size(ts1.rawprof,1) ~= numel(ts1.date) || size(ts2.rawprof,1) ~= numel(ts2.date) )
      error('Raw profile dimension 1 does not correspond to timestamps!');
    end;
  end;
  if ( dofield )
    if ( ndims(ts1.field) ~= 3 || ndims(ts2.field) ~= 3 ...
         || size(ts1.field,1) ~= numel(ts1.date) || size(ts2.field,1) ~= numel(ts2.date) )
      error('Field dimension 1 does not correspond to timestamps!');
    end;
  end;


  % Subset time series based on "IX1" and "IX2" arguments
  if ( isa(ix1,'function_handle') )
    ix1 = feval(ix1,ts1);
  end;
  if ( isa(ix2,'function_handle') )
    ix2 = feval(ix2,ts2);
  end;

  ts1.date = ts1.date(ix1);
  if dodata;	ts1.data = ts1.data(ix1,:);		end;
  if doprof;	ts1.prof = ts1.prof(ix1,:);		end;
  if dorawprof;	ts1.rawprof = ts1.rawprof(ix1,:);	end;
  if dofield;	ts1.field = ts1.field(ix1,:,:);		end;

  ts2.date = ts2.date(ix2);
  if dodata;	ts2.data = ts2.data(ix2,:);		end;
  if doprof;	ts2.prof = ts2.prof(ix2,:);		end;
  if dorawprof;	ts2.rawprof = ts2.rawprof(ix2,:);	end;
  if dofield;	ts2.field = ts2.field(ix2,:,:);		end;


  % Intersect time series to only matching timestamps
  [ix1,ix2] = intersect_dates(ts1.date,ts2.date,tol);

  if ( numel(unique(ts1.date(ix1))) ~= numel(unique(ts2.date(ix2))) )
    error('INTERSECT_DATES results mismatch! Try adjusting TOL?');
  end;

  ts1.date = ts1.date(ix1);
  if dodata;	ts1.data = ts1.data(ix1,:);		end;
  if doprof;	ts1.prof = ts1.prof(ix1,:);		end;
  if dorawprof;	ts1.rawprof = ts1.rawprof(ix1,:);	end;
  if dofield;	ts1.field = ts1.field(ix1,:,:);		end;
  ts2.date = ts2.date(ix2);
  if dodata;	ts2.data = ts2.data(ix2,:);		end;
  if doprof;	ts2.prof = ts2.prof(ix2,:);		end;
  if dorawprof;	ts2.rawprof = ts2.rawprof(ix2,:);	end;
  if dofield;	ts2.field = ts2.field(ix2,:,:);		end;

  opstr = 'NONSTRING';
  if ( ischar(op) )
    opstr = lower(op);
  end;

  % Produce result time series
  result.date = ts1.date;

  switch ( opstr ),

   case {'+','plus'},
    if dodata;		result.data = ts1.data + ts2.data;		end;
    if doprof;		result.prof = ts1.prof + ts2.prof;		end;
    if dorawprof;	result.rawprof = ts1.rawprof + ts2.rawprof;	end;
    if dofield;		result.field = ts1.field + ts2.field;		end;

   case {'-','minus'},
    if dodata;		result.data = ts1.data - ts2.data;		end;
    if doprof;		result.prof = ts1.prof - ts2.prof;		end;
    if dorawprof;	result.rawprof = ts1.rawprof - ts2.rawprof;	end;
    if dofield;		result.field = ts1.field - ts2.field;		end;

   case {'*','.*','times'},
    % NOTE: MATRIX operators treated as ARRAY operators!
    if dodata;		result.data = ts1.data .* ts2.data;		end;
    if doprof;		result.prof = ts1.prof .* ts2.prof;		end;
    if dorawprof;	result.rawprof = ts1.rawprof .* ts2.rawprof;	end;
    if dofield;		result.field = ts1.field .* ts2.field;		end;

   case {'/','./','rdivide','divide'},
    % NOTE: MATRIX operators treated as ARRAY operators!
    if dodata;		result.data = ts1.data ./ ts2.data;		end;
    if doprof;		result.prof = ts1.prof ./ ts2.prof;		end;
    if dorawprof;	result.rawprof = ts1.rawprof ./ ts2.rawprof;	end;
    if dofield;		result.field = ts1.field ./ ts2.field;		end;

   case {'^','.^','power'},
    % NOTE: MATRIX operators treated as ARRAY operators!
    if dodata;		result.data = ts1.data .^ ts2.data;		end;
    if doprof;		result.prof = ts1.prof .^ ts2.prof;		end;
    if dorawprof;	result.rawprof = ts1.rawprof .^ ts2.rawprof;	end;
    if dofield;		result.field = ts1.field .^ ts2.field;		end;

   otherwise,
    % NOTE: If the user really wants MATRIX operators on time series (but why
    % would they??), they must specify by name, e.g., TS_OP(TS1,TS2,'MTIMES')
    if dodata;		result.data = feval(op,ts1.data,ts2.data);		end;
    if doprof;		result.prof = feval(op,ts1.prof,ts2.prof);		end;
    if dorawprof;	result.rawprof = feval(op,ts1.rawprof,ts2.rawprof);	end;
    if dofield;		result.field = feval(op,ts1.field,ts2.field);		end;
  end;

  if ( doPlot && dodata )
    figure;
    maxigraph;
    plot(result.date,result.data);
    datetick3;
  end;

return;
