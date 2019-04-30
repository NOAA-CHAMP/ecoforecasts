function varargout = intersect_all_dates(tol, varargin)
%function varargout = intersect_all_dates(tol, varargin)
%
% Return indices of those elements in EACH of the datenum vectors specified
% in VARARGIN, that match among ALL those vectors: a "matching" timestamp is
% defined as one where the difference between the elements of the two series
% is less than some tolerance level TOL. If TOL empty, use DEFAULT of 30 mins
% and a hair: (30+(0.005/60))/(24*60). Length(varargout) == length(varargin),
% unless output is a single variable: then return a single cell array with
% all indices. In that case, length(varargout{1})== length(varargin).
%
% CALLS: INTERSECT_DATES.
%
% Last Saved Time-stamp: <Thu 2011-10-20 16:48:35  lew.gramer>

  if ( isempty(tol) )
    tol = (30.0+(0.005/60.0))/(24.0*60.0);
  end;
  if ( ~isnumeric(tol) || ~isscalar(tol) )
    error('First arg should be a numeric TOLERANCE');
  end;

  dts = varargin;
  ixs = {};

  origdts1 = dts{1};
  for ix = 2:length(dts)
    [ixs{1},ig] = intersect_dates(dts{1},dts{ix},tol);
    dts{1} = dts{1}(ixs{1});
    %DEBUG:    disp({'1 vs ix',numel(ixs{1})});
  end;
  ixs{1} = find(ismember(origdts1,dts{1}))';
  %DEBUG:  disp({'1',numel(ixs{1})});
  for ix = 2:length(dts)
    [ig,ixs{ix}] = intersect_dates(dts{1},dts{ix},tol);
    %DEBUG:    disp({ix,numel(ixs{ix})});
  end;

  % Sanity check
  for ix = 2:length(dts)
    if ( length(ixs{1}) ~= length(ixs{ix}) )
      error('Algorithm problem: Not all index results the same length!');
    end;
  end;

  if ( nargout == 1 )
    varargout(1) = {ixs};
  else
    for ix = 1:nargout
      varargout(ix) = {ixs{ix}};
    end;
  end;

return;
