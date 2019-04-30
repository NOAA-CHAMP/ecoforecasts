function h=rosewgtd(dirts,spdts,drs,cumfun)
%function h=rosewgtd(dirts,spdts,drs,cumfun)
%
% "Wind-rose" plot function for direction time series DIRTS. Calls ROSE (v.)
% and plots result with POLAR (v.), but ROSEWGTD accepts a time series of
% direction in degrees True (0o north, 90o east, rather than of a vector of
% east-origin directions in radians), and plots speed-weighted distributions
% of direction using speed time series SPDTS. Uses INTERSECT_TSES (v.) to
% combine the two time series. For purposes of plotting, ROSEWGTD bins dirs
% using HISTC (v.) and a vector of bin "edges" DRS (DEFAULT: [0:1:360]).  If
% SPDTS is not a valid time series, works very similarly to ROSE. If optional
% CUMFUN is a function handle, it is called in place of @NANMEAN.
%
% Last Saved Time-stamp: <Sun 2019-04-14 16:03:51 Eastern Daylight Time gramer>

  if ( ~is_ts(dirts) )
    error('First arg must be a valid time-series struct');
  end;
  if ( ~exist('spdts','var') || ~is_ts(spdts) )
    disp('No speed given: Binning by incidence');
    haveSpeed = false;
    dts = dirts;
  else 
    haveSpeed = true;
    [dts,sts] = intersect_tses(dirts,spdts);
  end;
  if ( ~exist('drs','var') || isempty(drs) )
    drs=0:1:359;
  end;
  if ( ~exist('cumfun','var') || isempty(cumfun) )
    cumfun = @nanmean;
  end;

  % Define the bin edges
  [N,B]=histc(dts.data,drs);
  % [N,B]=histcounts(dts.data,drs);

  % Preallocate direction bins for speed of execution

  if ( strcmpi(char(cumfun),'NANMEAN') )
    wd = repmat(nan,[1,ceil(nanmax(dts.data)*numel(drs))]);
  else
    % NANMEAN while more natural may lead to rounding errors
    wd = repmat(nan,[1,ceil(cumfun(dts.data)*numel(drs))]);
  end;

  endix = 0;
  for dr=1:numel(drs)
    n(dr) = 0;
    ix = find(B==dr);
    if ( numel(ix)>0 )
      if ( haveSpeed )
        n(dr) = round(cumfun(sts.data(ix)));
      else
        n(dr) = round(numel(ix));
      end;
    end;
    wd(endix+1:endix+n(dr)) = repmat(drs(dr),[1,n(dr)]);
    endix = endix + n(dr);
  end;

  % Remove excess NaNs when done
  wd(endix+1:end) = [];

  % [T,R]=rose(deg2rad(wd));
  % h=polar(T,R);
  % view(90,-90);

  [T,R]=rose(deg2rad(wd));
  h=polarplot(T,R);
  set(gca, 'ThetaDir','clockwise', 'ThetaZeroLocation','top');  

return;
