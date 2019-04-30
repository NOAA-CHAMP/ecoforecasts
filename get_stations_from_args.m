function stns = get_stations_from_args(varargin)
%function stns = get_stations_from_args([STN|STNS(:)|FIELD|COORDS(1:2,:)],extra_fldnms)
%
% Extract matrix of station STRUCTs with fields .lon,.lat, from the first
% arg, which should be a matrix of STRUCT with .lon,.lat or .lons,.lats, or a
% 1xN cell matrix of srings (or STRING vector), or else a CxN numeric or cell
% matrix (where C=2 + the number of names given in arg 2, EXTRA_FLDNMS).
%
% Last Saved Time-stamp: <Mon 2019-02-18 12:29:30 Eastern Standard Time gramer>

  args = varargin;
  
  stns=[];
  
  if ( isstruct(args{1}) )
    stns = args{1};
    
  elseif ( ischar(args{1}) )
    for ix=1:size(args{1},1)
      stns(ix) = get_station_from_station_name(args{1}(ix,:));
    end;
    
  elseif ( isstring(args{1}) )
    for ix=1:numel(args{1})
      stns(ix) = get_station_from_station_name(args{1}{ix});
    end;
    
  elseif ( iscellstr(args{1}) )
    for ix=1:size(args{1},1)
      stns(ix) = get_station_from_station_name(args{1}{ix});
    end;
    
  elseif ( isnumeric(args{1}) )
    for ix=1:size(args{1},2)
      stns(ix).lon = args{1}(1,ix);
      stns(ix).lat = args{1}(2,ix);
      if ( numel(args) > 1 )
        for fldix=1:numel(args{2})
          fldnm = args{2}{fldix};
          stns(ix).(fldnm) = args{1}(fldix+2,ix);
        end;
      end;
    end;
    
  elseif ( iscell(args{1}) )
    for ix=1:numel(args{1}{1})
      stns(ix).lon = args{1}{1}(ix);
      stns(ix).lat = args{1}{2}(ix);
      if ( numel(args) > 1 )
        for fldix=1:numel(args{2})
          fldnm = args{2}{fldix};
          stns(ix).(fldnm) = args{1}{fldix+2}{ix};
        end;
      end;
    end;
    
  end; %if ( isnumeric(args{1}) ) elseif ( isstruct(args{1}) ) else

return;
