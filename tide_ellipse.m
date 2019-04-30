1;
%% SCRIPT TIDE_ELLIPSE.m: Plot tide ellipses as calculated by TMD tide model toolkit
%  CALLS: ELLIPSE (TMD), ELLIPSE1 (Map Toolbox)
%  PARAMETERS: STNM - station name (DEFAULT 'mlrf1'); BATHRAD - radius of
%   bathymetric map (DEFAULT [40e3,40e3]); MDL - TMD tide model (DEFAULT
%   'Mex'); CMPS - tide components to plot (DEFAULT {'M2','K1','S2','O1'}).
% 
% CALLS: PLOT_CURRENT_ELLIPSE (Ecoforecasts)
% 
% Last Saved Time-stamp: <Fri 2017-03-10 16:33:41 Eastern Standard Time gramer>

    if ( ~exist('stnm','var') || isempty(stnm) )
      stnm = 'mlrf1';
    end;
    if ( ~exist('stn','var') || ~isfield(stn,'lon') || ...
         ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,stnm) )
      stn=[]; clear stn
      stn = get_station_from_station_name(stnm);
    end;

    if ( ~exist('bathrad','var') || isempty(bathrad) )
      bathrad = [40e3,40e3];
    end;
    if ( ~exist('bathcntrs','var') || isempty(bathcntrs) )
      bathcntrs = -[0:4:80];
    end;
    if ( ~exist('mdl','var') || isempty(mdl) )
      mdl = 'Mex';
      % mdl = 'tpxo7.2';
      % mdl = 'AO';
    end;
    if ( ~exist('cmps','var') || isempty(cmps) )
      cmps = {'M2','K1','S2','O1'};
    end;
    % % Primary components
    % cmp = 'M2';
    % cmp = 'K1';
    % % Secondary components
    % cmp = 'S2';
    % cmp = 'O1';

    %stn = plot_hires_bathymetry(stn,-[0:4:80],bathrad,[],[],false);
    if ( prod(bathrad(:)) > (50e3*50e3) )
      % For big domains, don't overtax the poor computer
      useHighestRes = false;
    else
      useHighestRes = true;
    end;
    stn = plot_hires_bathymetry(stn,bathcntrs,bathrad,true,@contour,[],[],useHighestRes);

    % TIDE_PRED only works easily if RUN IN the tmd_toolbox directory
    tide_mpath = which('ellipse');
    if ( isempty(tide_mpath) )
      error('No path found to ELLIPSE function!');
    end;
    [tide_mpath,ig] = fileparts(tide_mpath);
    %DEBUG:  dir(tide_mpath),

    clrs = {'k','b',[0,.5,0],'r','m','c','y',[.5,.5,.5]};
    nclrs = numel(clrs);

    for ix = 1:numel(cmps)
      cmp = cmps{ix};

      mydir = pwd;
      cd(tide_mpath);
      [umj(ix),umn(ix),uph(ix),uinc(ix)] = ellipse(fullfile('DATA',['Model_',mdl]),stn.lat,stn.lon,cmp);
      cd(mydir);

      % % [LAT,LON] = ellipse1(LAT0,LON0,ELLIPSE) computes ellipse(s) with
      % % center(s) at LAT0, LON0.  ELLIPSE must have the form [semimajor axis,
      % % eccentricity].  LAT0 and LON0 can be scalar or column vectors. The
      % % input and output latitudes and longitudes are in units of degrees.
      % % ELLIPSE must have the same number of rows as the input LAT0 and LON0.
      % % The semimajor axis (ELLIPSE(1)) is in degrees of arc length on a
      % % sphere. All ellipses are oriented so that their major axis runs
      % % north-south.
      % %
      % % [LAT,LON] = ellipse1(LAT0,LON0,ELLIPSE,OFFSET) computes the ellipses
      % % where the major axis is rotated from due north by an azimuth OFFSET.
      % % The offset angle is measure clockwise from due north.
      % [lats{ix},lons{ix}] = ellipse1(stn.lat,stn.lon,[umj(ix)/111,axes2ecc(umj(ix),abs(umn(ix)))],uinc(ix));

      [lats{ix},lons{ix}] = calc_current_ellipse(stn.lon,stn.lat,umj(ix),umn(ix),uinc(ix),uph(ix));
      plot(lons{ix},lats{ix},'.' ,'Color',clrs{mod(ix-1,nclrs)+1});
    end;

    legend('Bathymetry','Site','Coast',cmps{:});
    titlename(['Station ',upper(textize(stnm)),' model ',textize(mdl),' tide ellipses']);

    %stn=[]; clear ans stn stnm bathrad mdl clh clrs cmp cmps cs fh h ig ix lats lons mydir nclrs tide_mpath uinc umj umn uph
