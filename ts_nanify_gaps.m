function newtses = ts_nanify_gaps(tses,maxgap)
%function newtses = ts_nanify_gaps(tses,maxgap)
%
% Construct new time series STRUCT(s) NEWTSES containing all the data in each
% time series of TSES, but whenever FIND_DATE_RANGES (v.) finds a gap longer
% than MAXGAP days in any TS, insert a NaN at the beginning of that gap.
%
% This function is useful for example when calling PLOT_TS (v.) with linespec
% such as '.-':, since the insertion of NaNs will cause the removal of the
% line that appears to join points in the plot on either side of large gaps.
%
% As of 2018 Oct 01, invalid time series are PASSED THROUGH without change.
%
% DEFAULT MAXGAP: 3*MEDIAN(DIFF(TS.date)) calculated for each TS in TSES
%
% Last Saved Time-stamp: <Tue 2018-10-02 07:25:03 Eastern Daylight Time gramer>

  if ( ~isstruct(tses) )
    error('First arg must be a STRUCT or STRUCT matrix');
  end;

  if ( ~exist('maxgap','var') || isempty(maxgap) )
    default_maxgap = [];
  else
    default_maxgap = maxgap;
  end;

  for ix = 1:numel(tses)
    ts = tses(ix);

    newts = ts;
    
    if ( ~is_valid_ts(ts) )
      warning('Ecoforecasts:Gap:InvalidTS',...
              'TSES(%d) is not a valid time series STRUCT',ix);
      % PASS THROUGH invalid time series STRUCTs
      
    else
      
      % "Normal" time step in time series
      dt = median(diff(ts.date));
      
      if ( isempty(default_maxgap) )
        maxgap = 3*dt;
      else
        maxgap = default_maxgap;
      end;
      if ( maxgap < 2*dt )
        error('MAXGAP is too short for sampling frequency of time series TSES(%d)!',ix);
      end;
      
      doProf = false;
      if ( isfield(newts,'prof') && numel(newts.date) == size(newts.prof,1) )
        doProf = true;
      end;
      
      [rgs,ixes,allixes] = find_date_ranges(ts.date,maxgap);
      
      for ixix=size(ixes,2):-1:2
        newts.date(ixes(1,ixix)+1:end+1) = newts.date(ixes(1,ixix):end);
        newts.data(ixes(1,ixix)+1:end+1) = newts.data(ixes(1,ixix):end);
        if ( doProf )
          newts.prof(ixes(1,ixix)+1:end+1,:) = newts.prof(ixes(1,ixix):end,:);
        end;
        
        newts.date(ixes(1,ixix)) = newts.date(ixes(1,ixix))-dt;
        newts.data(ixes(1,ixix)) = nan;
        if ( doProf )
          newts.prof(ixes(1,ixix),:) = nan;
        end;
      end; %for ixix

    end; %if ( ~is_valid_ts(ts) ) else

    newtses(ix) = newts;
    ts=[]; newts=[]; clear ts newts
  end; %for ix

return;
