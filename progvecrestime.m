function [res_t,exit_locs,intermed_data] = progvecrestime(uts,vts,bath,varargin)
%function [res_t,exit_locs,intermed_data] = progvecrestime(uts,vts,bath,varargin)
%
% Use Progressive Vector Diagram method (e.g., Fiechter et al. 2008) to
% estimate residence time from point current measurement time series UTS and
% VTS, using the bathymetry data field in BATH.lon, BATH.lat, BATH.field.
%
% Bounds of the domain for residence time estimates is defined by optional
% MAXISO isobath (DEFAULT: 30 m), the shallows inshore (MINISO, DEFAULT: 0
% m), and by a maximum distance from the center location BATH.lon,BATH.lat
% (DEFAULT: edges of domain in BATH). Optional arg METHOD==1 (DEFAULT) means
% use a fast, memory-intensive matrix-arithmetic method, or 2 use a slower,
% more straightforward method of looping over the time points in UTS/VTS. 
%
% Last Saved Time-stamp: <Tue 2017-03-21 17:34:42 Eastern Daylight Time gramer>

  if ( nargin < 2 || ~is_ts(uts) || ~is_ts(vts) )
    error('First two args should be time series STRUCTs (v. IS_TS)');
  end;
  if ( nargin < 3 || ~isfield(bath,'lon') || ~isfield(bath,'lat') || ~isfield(bath,'field') )
    error('Third arg should be a bathymetry STRUCT (v. READ_HIRES_BATHYMETRY)');
  end;

  args = varargin;
  [start_lon,args]      = getarg(args,'start_lon','default',mean(bath.lon));
  [start_lat,args]      = getarg(args,'start_lat','default',mean(bath.lat));
  [method,args]         = getarg(args,'method','default',1);
  [maxiso,args]         = getarg(args,'maxiso','default',30);
  [miniso,args]         = getarg(args,'miniso','default',0);
  [tgtpts,args]         = getarg(args,'tgtpts','default',7000);
  [maxres_t,args]       = getarg(args,'maxres_t','default',60);

  res_t=[];
  res_t.date = uts.date;
  res_t.data = repmat(nan,size(uts.data));

  if ( nargout > 1 )
    exit_locs = repmat(nan,[numel(uts.data),2]);
  end;
  if ( nargout > 2 )
    intermed_data=[];
  end;

  % Median number of points needed to encompass maximum expected residence time
  maxres_t_pts = ceil(maxres_t / nanmedian(diff(uts.date)));
  %DEBUG:
  disp(['Max Res_t ',num2str(maxres_t),' d := ',num2str(maxres_t_pts),' points']);

  cmperlat = 1e2 * 1e3 * distance_wgs84(start_lat,start_lon,start_lat+1,start_lon); %[cm / deg Lat]
  cmperlon = 1e2 * 1e3 * distance_wgs84(start_lat,start_lon,start_lat,start_lon+1); %[cm / deg Lon]

  % Time between data points [s]
  delt = 24*3600*diff(uts.date); delt(end+1) = median(delt);

  %% METHOD 1: Memory intensive but faster (with caveats)
  % Reduce memory and thus improve speed by splitting task into sections
  if ( method == 1 )
      all_uts = uts;
      all_vts = vts;

      all_res_t.date = all_uts.date;
      all_res_t.data = repmat(nan,[numel(all_uts.date),1]);

      nsecs = ceil(numel(all_uts.date)/tgtpts);
      for secix=1:nsecs;
        %DEBUG:
        disp([num2str(secix),' of ',num2str(nsecs),' sections']);

        begix = ((secix-1)*tgtpts) + 1;
        endix = begix + tgtpts + maxres_t_pts;
        % On the last section, count backward from the end
        if ( secix == nsecs )
          endix = numel(all_uts.data);
          %begix = min([begix,endix-maxres_t_pts-1]);
          begix = endix - tgtpts - maxres_t_pts - 1;
        end;
        %DEBUG:
        disp([num2str(begix),':',num2str(endix)]);
        begix = max([begix,1]);
        if ( endix >= numel(all_uts.data) )
          endix = numel(all_uts.data);
          secix = nsecs;
        end;
        %DEBUG:
        disp([num2str(begix),':',num2str(endix)]);
        secixen = begix:endix;

        uts = subset_ts(all_uts,secixen);
        vts = subset_ts(all_vts,secixen);
        delt = 24*3600*diff(uts.date); delt(end+1) = median(delt);

        UTS=[]; VTS=[]; DELT=[]; INTEGX=[]; INTEGY=[]; INTEGZ=[]; IX=[]; JX=[];
        clear UTS VTS DELT INTEGX INTEGY INTEGZ IX JX

        % Rows: t for u(t), Cols: tau for res_t(tau)
        UTS = triu(repmat(uts.data',[numel(uts.data)-1,1]));
        VTS = triu(repmat(vts.data',[numel(vts.data)-1,1]));

        DELT = repmat(delt',[numel(uts.data)-1,1]);
        DX = UTS.*DELT;
        DY = VTS.*DELT;
        DELT=[]; UTS=[]; VTS=[]; clear DELT UTS VTS

        INTEGX = start_lon + ((cumsum(DX,2))./cmperlon);
        INTEGY = start_lat + ((cumsum(DY,2))./cmperlat);
        DX=[]; DY=[]; clear DX DY

        INTEGZ = interp2(bath.lon,bath.lat,bath.field,INTEGX,INTEGY,'nearest');

        [IX,JX] = find(-maxiso>INTEGZ | INTEGZ>=miniso | isnan(INTEGZ));
        if ( nargout > 2 )
          intermed_data.INTEGX = INTEGX;
          intermed_data.INTEGY = INTEGY;
          intermed_data.INTEGZ = INTEGZ;
        end;
        if ( nargout <= 1 )
          INTEGX=[]; INTEGX=[]; INTEGZ=[]; clear INTEGX INTEGY INTEGZ
        end;

        for ix=unique(IX)'
          jx = find(IX==ix,1);
          jx = JX(jx);
          if ( nargout > 1 )
            exit_locs(ix,1) = INTEGX(jx);
            exit_locs(ix,2) = INTEGY(ix);
          end;
          res_t.data(ix) = uts.date(jx) - uts.date(ix);
        end;
        INTEGX=[]; INTEGX=[]; INTEGZ=[]; clear INTEGX INTEGY INTEGZ
        if ( nargout > 2 )
          intermed_data.IX = IX;
          intermed_data.JX = JX;
        end;
        IX=[]; JX=[]; clear IX JX ix jx

        matchix = find(ismember(all_res_t.date,res_t.date));
        all_res_t.data(matchix) = res_t.data;

        if ( secix == nsecs )
          break;
        end;
      end; %for secix...

      res_t = all_res_t;
      all_res_t=[]; clear all_res_t


    %% METHOD 2: Simple loop: slower, but more straightforward and less memory
    else

      %DEBUG:
      disp(numel(uts.data));
      if ( nargout > 2 )
        intermed_data.integx = repmat(nan,size(uts.data));
        intermed_data.integy = repmat(nan,size(uts.data));
        intermed_data.integz = repmat(nan,size(uts.data));
        intermed_data.deepix = repmat(nan,size(uts.data));
      end;

      for ix = 1:numel(uts.data)
        integx = start_lon+(cumsum(delt(ix:end).*uts.data(ix:end))./cmperlon);
        integy = start_lat+(cumsum(delt(ix:end).*vts.data(ix:end))./cmperlat);
        integz = interp2(bath.lon,bath.lat,bath.field,integx,integy);
        deepix = find( (-maxiso>integz | integz>=miniso | isnan(integz)), 1) + ix - 1;
        % Should only be empty in the final days/weeks of a deployment
        if ( ~isempty(deepix) )
          res_t.data(ix) = uts.date(deepix) - uts.date(ix);
        end;
        if ( nargout > 2 )
          intermed_data.integx(ix) = integx;
          intermed_data.integy(ix) = integy;
          intermed_data.integz(ix) = integz;
          intermed_data.deepix(ix) = deepix;
        end;
        integx=[]; integy=[]; integz=[]; deepix=[]; clear integx integy integz deepix ix
      end; %for ix = 1:numel(uts.data)

    end;

return;
