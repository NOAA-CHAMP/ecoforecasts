function varargout = intersect_tses(varargin)
%function varargout = intersect_tses([tol,]ts1,ts2,...)
%
% Return vector of time series structs with only those datapoints that
% coincide among EACH of the input time series specified in VARARGIN: a
% "matching" timestamp is defined as for INTERSECT_ALL_DATES (v.). TOL
% tolerance (in days) may be specified as an optional first arg. Otherwise,
% must have length(varargout) == length(TSES), unless output is a single cell
% array variable: then return a single cell array with all TS structs in it.
% In that case, length(varargout{1})== length(TSES).
%
% As of 2017, now also subsets .prof and .rawprof (for, e.g., ADCP profiles),
% and .field (for, e.g., 2-D SST fields) based on intersection. But first,
% verifies that SIZE(...,1) of each of these fields matches NUMEL(TS.date).
%
% CALLS: INTERSECT_ALL_DATES.
%
% Last Saved Time-stamp: <Mon 2017-07-03 15:56:33 Eastern Daylight Time gramer>

  argix = 1;
  if ( isnumeric(varargin{1}) )
    tol = varargin{1};
    argix = argix + 1;
  else
    tol = [];
  end;

  % Extract time series structs to intersect dates for
  tses = {};
  ntses = 0;
  while ( argix <= nargin )
    arg = varargin{argix};
    if ( isstruct(arg) && isfield(arg,'date') && isfield(arg,'data') )
      % Gracefully deal with vectors OR MATRICES of structs passed as a single argument
      for strcix = 1:numel(arg)
        strc = arg(strcix);
        if ( isempty(strc.date) || ~isfield(strc,'date') || ...
             ~isfield(strc,'data') || any(size(strc.date) ~= size(strc.data)) )
          error('Struct arguments must be non-empty time series!');
        end;
        ntses = ntses + 1;
        tses{ntses} = strc;
        dts{ntses} = strc.date;
      end;
    else
      error('All args after optional TOL must be time series structs or arrays of TS structs!');
    end;
    argix = argix + 1;
  end;

  ixes = intersect_all_dates(tol,dts{:});
  for ix = 1:numel(ixes)
    if ( isfield(tses{ix},'prof') )
      if ( ndims(tses{ix}.prof) == 2 && size(tses{ix}.prof,1) == numel(tses{ix}.date) )
        tses{ix}.prof = tses{ix}.prof(ixes{ix},:);
      end;
    end;
    if ( isfield(tses{ix},'rawprof') )
      if ( ndims(tses{ix}.rawprof) == 2 && size(tses{ix}.rawprof,1) == numel(tses{ix}.date) )
        tses{ix}.rawprof = tses{ix}.rawprof(ixes{ix},:);
      end;
    end;
    if ( isfield(tses{ix},'field') )
      if ( ndims(tses{ix}.field) == 3 && size(tses{ix}.field,1) == numel(tses{ix}.date) )
        tses{ix}.field = tses{ix}.field(ixes{ix},:,:);
      end;
    end;
    if ( isfield(tses{ix},'data') )
      tses{ix}.data = tses{ix}.data(ixes{ix});
    end;
    tses{ix}.date = tses{ix}.date(ixes{ix});
  end;

  if ( nargout == 0 )
  elseif ( nargout == 1 )
    varargout{1} = tses;
  elseif ( nargout == numel(tses) )
    for ix = 1:nargout
      varargout{ix} = tses{ix};
    end;
  else
    error('Need one output arg per input TS, or a single cell array');
  end;

return;
