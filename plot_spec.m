function [Pxxes,W_days,fh,lhs,Pxxces] = plot_spec(stn,var_or_vars,specfunc,plotfunc,xlm,ylm,printFmt,doNorm,presVar,ixes,confInt,bandAvg)
%function [Pxxes,W_days,fh,lhs,Pxxces] = plot_spec(stn,[var|vars],specfunc,plotfunc,xlm,ylm,printFmt,doNorm,presVar,ixes,confInt,bandAvg)
% 
% Calculate, plot and return a Power Spectral Density estimate (spectrum) for
% a time series named STN.(VAR). Fields STN.(VAR).date and .data must exist.
% VARS may also be a CELLSTR of field names, each plotted in its own color.
% Uses function SPECFUNC (DEFAULT: @pmtm) to calculate the PSD, and  function
% PLOTFUNC (DEFAULT: @loglog) to plot it. If SPECFUNC is a function handle or
% a cell with first elt. function handle, that function must take one arg and
% return two values: PXX (vector of power per unit frequency estimates) and W
% (frequencies), both which are returned by PLOT_SPEC: see, e.g., PWELCH,
% PMTM, PERIODOGRAM, etc. The returned W_DAY is the vector of periods in
% *DAYS*. FH is a new FIGURE handle, LHS a column vector of new lineseries
% handles (v. PLOT): if PLOTFUNC is 'none', then no plot is generated; FH and
% LHS will then be empty.
%
% Optional XLM and YLM specify lower and upper bounds for period and energy,
% respectively. If either is 'auto', use natural data limits (DEFAULTS:
% [2/24,365*25] and [1e-7,1e7], resp.) Optional PRINTFMT specifies format to
% PRINT figure in (e.g.,'tiff'); if empty (DEFAULT), no printout is created.
%
% If DONORM (DEFAULT: false), subtract mean of all variables before PSD; if
% PRESVAR (DEFAULT: false), plot variance-preserving spectra (W*Pxx vs logW).
% In this latter case, the DEFAULT value of YLM becomes 'auto' (see above).
% If VARS is a cellstr with more than one element, PXXES and W_DAYS are cell
% arrays of numerical vectors containing the same number of elements as VARS.
% Otherwise, PXXES and W_DAYS are simple numerical vectors. If optional IXES
% is an index vector, only plot spectrum using those indices in STN.(VAR).
% If it is a FUNCTION_HANDLE, plot indices returned by IXES(STN.(VAR)); if
% VARS is a CELLSTR and IXES is given, NUMEL(IXES) must equal NUMEL(VARS).
% Optional output PXXCES is a confidence interval (see PMTM) for spectrum.
% If optional CONFINT is True, then plot that confidence interval. If CONFINT
% is numeric, use that Confidence Interval. BANDAVG if an integer > 1, asks
% that spectrum be band-averaged: each point on the plot will contain the
% average of BANDAVG points from the natural spectrum.
%
% SAMPLE CALLS below both plot 68% CONFINT; the first also passes a sampling
% frequency and a PSD combination method into spectral function PMTM (v.);
% the second band-averages the PSD plot every 3 periods. Both subtract the
% mean before calculating, and plot a variance-preserving spectrum:
%
% >> plot_spec(stn,'ndbc_sea_t',{@pmtm,[],[],2*pi,'adapt',0.68},[],[],[],[],true,true,[],true);
% >> plot_spec(stn,'ndbc_sea_t',[],[],[],[],[],true,true,[],0.68,3);
% 
% Last Saved Time-stamp: <Thu 2018-08-02 12:54:53 Eastern Daylight Time gramer>

  Pxxes = [];
  W_days = [];
  W_hours = [];
  fh = [];
  lhs = [];

  figspath = get_ecoforecasts_path('figs');

  if ( iscellstr(var_or_vars) )
    vars = var_or_vars;
  elseif ( ischar(var_or_vars) )
    vars = {var_or_vars};
  else
    error('Second arg must be a valid field name or cellstr of them');
  end;
  clear var_or_vars;

  if ( ~exist('specfunc','var') || isempty(specfunc) )
    specfunc = @pmtm;
  end;
  specargs = {};
  if ( iscell(specfunc) )
    if ( numel(specfunc) > 1 )
      specargs = specfunc(2:end);
    end;
    specfunc = specfunc{1};
  end;
  if ( ~exist('presVar','var') || isempty(presVar) )
    presVar = false;
  end;

  xscl = '';
  yscl = '';
  if ( ~exist('plotfunc','var') || isempty(plotfunc) )
    if ( presVar )
      plotfunc = @semilogx;
      xscl = 'log';
      yscl = 'linear';
    else
      plotfunc = @loglog;
      xscl = 'log';
      yscl = 'log';
    end;
  end;
  if ( ~exist('xlm','var') || isempty(xlm) )
    xlm = [(2/24),(365*25)];
  end;
  if ( ~exist('ylm','var') || isempty(ylm) )
    if ( presVar )
      ylm = 'auto';
    else
      ylm = [1e-7,1e+7];
    end;
  end;
  if ( ~exist('printFmt','var') )
    printFmt = '';
  end;
  if ( ~exist('doNorm','var') || isempty(doNorm) )
    doNorm = false;
  end;

  if ( ~exist('ixes','var') || isempty(ixes) )
    ixes = [];
  end;
  if ( ~iscell(ixes) )
    ixes = {ixes};
  end;
  if ( numel(ixes) == 1 )
    ixes = repmat(ixes,[1,numel(vars)]);
  elseif ( numel(ixes) ~= numel(vars) )
    error('If IXES is given, NUMEL(IXES) must equal NUMEL(VARS)');
  end;

  if ( ~exist('confInt','var') || isempty(confInt) )
    confInt = false;
  end;

  annotStr = '';

  if ( confInt )
    % If user wants confidence interval other than 95%, supply arguments to PMTM if needed
    if ( isempty(specargs) && strcmpi(char(specfunc),'pmtm') && isnumeric(confInt) && confInt > 0 )
      specargs = {[],[],2*pi,'adapt',confInt};
      annotStr = num2str(confInt*100,' CI %.1f%%');
    elseif ( numel(specargs) >= 5 )
      annotStr = num2str(specargs{5}*100,' CI %.1f%%');
    else
      annotStr = ' CI 95.0%';
    end;
  end;

  if ( ~exist('bandAvg','var') || isempty(bandAvg) )
    bandAvg = 1;
  end;
  if ( bandAvg < 1 || floor(bandAvg) ~= bandAvg )
    error('Optional arg BANDAVG must be an integer >= 1!');
  end;

  if ( isfield(stn,'station_name') )
    stnm = stn.station_name;
  elseif ( isfield(stn,'name') )
    stnm = stn.name;
  elseif ( length(inputname(1)) > 0 )
    stnm = inputname(1);
  else
    stnm = 'stn';
  end

  % Just for fun, let's output the local inertial period if we calculate it
  if ( ~isfield(stn,'lat') || ~exist('inertial_period') > 1 )
    IP_hr = [];
    IP_day = [];
    warning('Inertial Period not calculated');
  else
    %% Where the heck did these formulae come from?? They're wrong!
    %IP_hr = 6/(1.4584e-4*sind(stn.lat)*3600);
    %IP_day = 6/(1.4584e-4*sind(stn.lat)*24*3600);

    [IP_day,IP_hr,IP_sec] = inertial_period(stn.lat);
    disp(['IP(' num2str(stn.lat) ') = ' num2str(IP_hr) ' hours']);
    disp(['IP(' num2str(stn.lat) ') = ' num2str(IP_day) ' days']);
  end;


  if ( ~strcmpi(plotfunc,'none') )
    fh = figure;
    maxigraph;
    % hold on;
    ax = gca;
    set(gca,'FontSize',20);  % For publication-ready fonts when printing
    co = get(ax,'ColorOrder');
    ncos = size(co,1);
  end;


  min_W_day = +Inf;
  max_W_day = -Inf;
  min_Pxx = +Inf;
  max_Pxx = -Inf;


  varnms = '';
  nvars = numel(vars);
  for varix=1:nvars
    var = vars{varix};
    varnms = [varnms,'-',var];

    per_hr = 1/(median(diff(stn.(var).date))*24);
    per_day = 1/(median(diff(stn.(var).date)));

    %if ( isempty(ixes) )
    if ( isempty(ixes{varix}) )
      ix = 1:numel(stn.(var).data);
    elseif ( isa(ixes{varix},'function_handle') )
      ix = ixes{varix}(stn.(var));
    else
      ix = ixes{varix};
    end;

    dat = stn.(var).data(ix);
    % NOTE: Fill in NaNs with mean value!
    dat(~isfinite(dat)) = nanmean(dat);

    if ( doNorm )
      disp(['Analyzing ',var,' with mean removed']); 
      dat = dat - nanmean(dat);
    end;

    switch ( upper(char(specfunc)) ),
     case {'PMTM'},
      % Only one function supports confidence intervals in MATLAB 7
      [Pxx, Pxxc, W] = feval(specfunc, dat, specargs{:});
     otherwise,
      [Pxx, W] = feval(specfunc, dat, specargs{:});
      Pxxc = [];
    end;

    Pxx(W<eps) = [];
    if ( ~isempty(Pxxc) )
      Pxxc(W<eps,:) = [];
    end;
    W(W<eps) = [];

    if ( bandAvg > 1 )
      %DEBUG:
      disp(['Averaging every ',num2str(bandAvg),' frequencies']);
      annotStr = [annotStr,' ',num2str(bandAvg),'-band avg'];
      nrows = floor(numel(W)/bandAvg);
      avgix = 1:(nrows*bandAvg);
      W = nanmean( reshape(W(avgix)',[bandAvg,nrows]) )';
      Pxx = nanmean( reshape(Pxx(avgix)',[bandAvg,nrows]) )';
      if ( ~isempty(Pxxc) )
        Pxxc1 = nanmean( reshape(Pxxc(avgix,1)',[bandAvg,nrows]) )';
        Pxxc2 = nanmean( reshape(Pxxc(avgix,2)',[bandAvg,nrows]) )';
        clear Pxxc
        Pxxc = [Pxxc1,Pxxc2];
      end;
    end;

    W_hour = (2*pi/per_hr) ./ W;
    W_day = (2*pi/per_day) ./ W;

    if ( presVar )
      disp(['Spectrum of ',var,' preserves variance']); 
      % W_day is in units of Days per Cycle (DPC), so units for
      % variance-preserving spectrum should be Freq*PSD, or PSD/DPC.
      Pxx = Pxx./W_day;
      %Pxx = Pxx./W_hour;
      if ( ~isempty(Pxxc) )
        Pxxc(:,1) = Pxxc(:,1)./W_day;
        Pxxc(:,2) = Pxxc(:,2)./W_day;
      end;
    end;

    min_W_day = min(min_W_day,min(W_day));
    max_W_day = max(max_W_day,max(W_day));
    min_Pxx = min(min_Pxx,min(Pxx));
    max_Pxx = max(max_Pxx,max(Pxx));

    if ( nvars > 1 )
        Pxxes{end+1} = Pxx;
        Pxxces{end+1} = Pxxc;
        W_hours{end+1} = W_hour;
        W_days{end+1} = W_day;
    else
        Pxxes = Pxx;
        Pxxces = Pxxc;
        W_hours = W_hour;
        W_days = W_day;
    end;

    if ( ~isempty(fh) )
      lhs(end+1) = feval(plotfunc, W_day,Pxx, 'Color',co(mod(varix-1,ncos)+1,:));
      %lhs(end+1) = feval(plotfunc, W_hour,Pxx, 'Color',co(mod(varix-1,ncos)+1,:));

      %%%%??? HACK necessitated by idiotic printing bug in MATLAB: Call
      %HOLD ON before first LOGLOG call, then print TIFF from this figure,
      %and axes in the printed version will have no top or right border!
      hold on;

      if ( confInt && ~isempty(Pxxc) )
        feval(plotfunc, W_day,Pxxc(:,1), 'LineStyle',':', 'Color',co(mod(varix-1,ncos)+1,:));
        feval(plotfunc, W_day,Pxxc(:,2), 'LineStyle',':', 'Color',co(mod(varix-1,ncos)+1,:));
      end;

      if ( isempty(xscl) )
        % Save these after first call: LOGLOG/SEMILOGX can be squirrelly
        xscl = get(ax,'XScale');
        yscl = get(ax,'YScale');
      end;
    end;
  end;

  if ( ~isempty(fh) )
    % Just make sure: calling LOGLOG/SEMILOGX multiple times can be squirrelly
    set(ax,'XScale',xscl);
    set(ax,'YScale',yscl);

    set(ax,'FontSize',14);

    if ( strcmpi(xlm,'auto') )
      if ( strcmpi(xscl,'log') )
        xlm = [ power( 10, floor(log10(min_W_day))-1 ), ...
                power( 10, ceil(log10(max_W_day))+1 ) ];
      else
        xlm = [ (min_W_day.*0.75), ...
                (max_W_day.*1.25) ];
      end;
    end;
    if ( strcmpi(ylm,'auto') )
      if ( strcmpi(yscl,'log') )
        ylm = [ power( 10, floor(log10(min_Pxx))-1 ), ...
                power( 10, ceil(log10(max_Pxx))+1 ) ];
      else
        ylm = [ (min_Pxx.*0.75), ...
                (max_Pxx.*1.25) ];
      end;
    end;
    xlim(xlm);
    ylim(ylm);
    grid on;
    set(ax,'XMinorGrid','off');
    set(ax,'YMinorGrid','off');

    % hold on;
    pwrlims = ylim;
    valign = 'top';
    pers = sort([natural_periods,IP_day]);
    for per = pers;
      perlh(1)=line([per per],[pwrlims(1) pwrlims(1)*1e2], 'Color','black');
      perlh(2)=line([per per],[pwrlims(2)*1e-2 pwrlims(2)], 'Color','black');
      % Labeling all tidal frequencies would be too thick
      if ( per>2 || ismember(per,[0.5,1,IP_day]) )
        str = [num2str(per,'%.3g') 'd'];
        if ( per == IP_day );
          set(perlh,'LineWidth',1.5,'LineStyle','--');
          str = ['IP=' str];
        end;
        str = ['  ' str];
        text(per,pwrlims(1),str, 'Rot',90, 'Vert',valign, 'Color','red', 'FontSize',14);
        if (strcmp(valign,'top')); valign='bot'; else valign='top'; end;
      end;
    end;

    if ( doNorm )
      optstr = ' normalized';
      fnmstr = '-norm';
    else
      optstr = '';
      fnmstr = '';
    end;
    if ( presVar )
      optstr = [optstr,' variance-preserving'];
      fnmstr = [fnmstr,'-varpres'];
    end;

    if ( numel(vars) == 1 )
      titlename(strrep(...
          sprintf('%s Spectrum %s%s %s %s',upper(char(specfunc)),upper(stnm),optstr,upper(vars{:}),annotStr),...
          '_','\_'));
    else
      legend(strrep(upper(vars),'_','\_'), 'Location','SouthEast');
      titlename(strrep([upper(char(specfunc)),' Spectrum: ',upper(stnm),optstr],'_','\_'));
    end;

    if ( ~isempty(printFmt) )
      if ( strcmpi(printFmt,'jpg') )
        printdev = '-djpeg';
      elseif ( strcmpi(printFmt,'tif') )
        printdev = '-dtiff';
      else
        printdev = ['-d' lower(printFmt)];
      end;
      set(lhs,'LineWidth',1);
      figfname = fullfile( figspath,...
                           sprintf('%s%s-%s%s.%s',lower(stnm),varnms,...
                                   lower(mfilename),fnmstr,lower(printFmt)) );
      print(printdev, figfname);
    end;
  end;

return;
