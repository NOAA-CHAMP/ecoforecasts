function station_boxplots(stn,fld,ylbl,ylm,doyrs,domos,tol,doPlot,doPrint)
%function station_boxplots(stn,fld,ylbl,ylm,doyrs,domos,tol,doPlot,doPrint)
% 
% Call BOXPLOT (qv.) for each of annual (interannual var.), monthly, weekly
% (both annual), and hourly (diurnal) distributions of time series STN.(FLD).
% If optional DOYRS given, only use data from those years; if DOMOS given,
% only use data from those months - for both annual and interannual var.
% If 0 < TOL <= 1 given (DEFAULT: 0.75), only use year-months AND years
% having at least that fraction of total possible periodic values, e.g., if
% TOL==0.75, remove any year-months with < 0.75*28*24 hourly values, then
% remove any years with < 0.75*12*28*24 *remaining* hourly values. For
% non-hourly periodic data (e.g., daily means), sample period is assumed to
% be a constant, equal to MIN(DIFF(UNIQUE(STN.(FLD).date))). Optional DOPLOT
% is a logical 4-vector (DEFAULT all true) to plot annual, monthly, weekly,
% and hourly distributions, resp. DOPRINT is a 4-vector controlling printing
% (in TIFF format) of each of the four BOXPLOTs.
% 
% Last Saved Time-stamp: <Sun 2011-10-30 13:31:52  lew.gramer>


  figspath = get_ecoforecasts_path('figs');

  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = strrep(fld,'_','\_');
  end;
  if ( ~exist('ylm','var') )
    ylm = [];
  end;
  if ( ~exist('doyrs','var') )
    doyrs = [];
  end;
  if ( ~exist('domos','var') )
    domos = [];
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    tol = 0.75;  % QC tolerance = 75% of expected # hourly values
  end;

  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = true;
  end;
  if ( isscalar(doPlot) )
    doPlot = logical([ doPlot doPlot doPlot doPlot ]);
  end;
  if ( numel(doPlot) ~= 4 )
    error('Optional arg DOPLOT must be a logical scalar or 4-vector');
  end;

  if ( ~exist('doPrint','var') || isempty(doPrint) )
    doPrint = false;
  end;
  if ( isscalar(doPrint) )
    doPrint = logical([ doPrint doPrint doPrint doPrint ]);
  end;
  if ( numel(doPrint) ~= 4 )
    error('Optional arg DOPRINT must be a logical scalar or 4-vector');
  end;

  stnm = lower(stn.station_name);

  samples_per_day = 1/min(diff(unique(stn.(fld).date)));

  % Filter time series dates based on user request
  [yrs,mos,dys,hrs,ig,ig] = datevec(stn.(fld).date);
  wks = get_week(stn.(fld).date);
  dat = stn.(fld).data;

  if ( isempty(doyrs) )
    doyrs = unique(yrs);
  else
    doyrix = find(ismember(yrs,doyrs));
    yrs = yrs(doyrix);
    mos = mos(doyrix);
    wks = wks(doyrix);
    dys = dys(doyrix);
    hrs = hrs(doyrix);
    dat = dat(doyrix);
  end;
  if ( isempty(domos) )
    domos = 1:12;
  else
    domoix = find(ismember(mos,domos));
    yrs = yrs(domoix);
    mos = mos(domoix);
    wks = wks(domoix);
    dys = dys(domoix);
    hrs = hrs(domoix);
    dat = dat(domoix);
  end;

  if ( isempty(dat) )
    error('No data for the requested years and months!');
  end;

  % Simple quality control on means, using run-length encoding
  % Only do months with TOL % of possible hourly values
  kern = find(mos(1:end-1) ~= mos(2:end));
  rl = diff([ 0 kern(:)' length(mos) ]);
  vs = mos([ kern(:)' length(mos) ]);
  goodmos = vs(rl >= tol*28*samples_per_day);

  badmos = vs(rl < tol*28*samples_per_day);
  if ( ~isempty(badmos) )
    vs = yrs([ kern(:)' length(yrs) ]);
    badyrs = vs(rl < tol*28*samples_per_day);
    disp('Excluding year-months:');
    disp(num2str(unique((badyrs(:)*100)+badmos(:))));
  end;

  goodix = find(ismember(mos,goodmos));
  yrs = yrs(goodix);
  mos = mos(goodix);
  wks = wks(goodix);
  dys = dys(goodix);
  hrs = hrs(goodix);
  dat = dat(goodix);

  % Only do years with TOL % of possible hourly values
  kern = find(yrs(1:end-1) ~= yrs(2:end));
  rl = diff([ 0 kern(:)' length(yrs) ]);
  vs = yrs([ kern(:)' length(yrs) ]);
  goodyrs = vs(rl >= tol*length(domos)*28*samples_per_day);

  badyrs = vs(rl < tol*length(domos)*28*samples_per_day);
  if ( ~isempty(badyrs) )
    disp('Excluding years:');
    disp(num2str(unique(badyrs)));
  end;

  goodix = find(ismember(yrs,goodyrs));
  yrs = yrs(goodix);
  mos = mos(goodix);
  wks = wks(goodix);
  dys = dys(goodix);
  hrs = hrs(goodix);
  dat = dat(goodix);


  if ( isempty(dat) )
    error('No data left after quality control!');
  end;


  % Do box-and-whisker plots on annual distributions
  if ( doPlot(1) )
    figure;
    missing_yrs = setxor(min(yrs):max(yrs),unique(yrs));
    yyrs = [yrs(:) ; missing_yrs(:)];
    ydat = [dat(:) ; repmat(0,[numel(missing_yrs),1])];
    boxplot(ydat, yyrs, 'notch','on', 'whisker',2);
    xlabel('Year');
    ylabel(ylbl);
    titlename([upper(stnm) ' Interannual Boxplots - ' ylbl]);
    if ( ~isempty(ylm) )
      ylim(ylm);
    end;
    maxigraph;
    grid on;
    hold on;
    if ( doPrint(1) )
      print('-dtiff',fullfile(figspath,sprintf('%s-%s-interannual-boxplot.tiff',stnm,fld)));
    end;
  end;


  % Do box-and-whisker plots on monthly distributions
  if ( doPlot(2) )
    figure;
    missing_mos = setxor(1:12,unique(mos));
    mmos = [mos(:) ; missing_mos(:)];
    mdat = [dat(:) ; repmat(0,[numel(missing_mos),1])];
    boxplot(mdat, mmos, 'notch','on', 'whisker',2);
    xlabel('Year-Month');
    ylabel(ylbl);
    titlename([upper(stnm) ' Monthly Boxplots - ' ylbl]);
    xlim([0 13]);
    if ( ~isempty(ylm) )
      ylim(ylm);
    end;
    maxigraph;
    grid on;
    hold on;
    if ( doPrint(2) )
      print('-dtiff',fullfile(figspath,sprintf('%s-%s-monthly-boxplot.tiff',stnm,fld)));
    end;
  end;


  % Do box-and-whisker plots on weekly distributions
  if ( doPlot(3) )
    figure;
    missing_wks = setxor(1:52,unique(wks));
    wwks = [wks(:) ; missing_wks(:)];
    wdat = [dat(:) ; repmat(0,[numel(missing_wks),1])];
    boxplot(wdat, wwks, 'notch','on', 'whisker',2);
    xlabel('Year-Week');
    ylabel(ylbl);
    titlename([upper(stnm) ' Weekly Boxplots - ' ylbl]);
    xlim([0 53]);
    if ( ~isempty(ylm) )
      ylim(ylm);
    end;
    maxigraph;
    grid on;
    hold on;
    if ( doPrint(3) )
      print('-dtiff',fullfile(figspath,sprintf('%s-%s-weekly-boxplot.tiff',stnm,fld)));
    end;
  end;


  % Do box-and-whisker plots on hourly distributions
  if ( doPlot(4) )
    figure;
    missing_hrs = setxor(0:23,unique(hrs));
    hhrs = [hrs(:) ; missing_hrs(:)];
    hdat = [dat(:) ; repmat(0,[numel(missing_hrs),1])];
    boxplot(hdat, hhrs, 'notch','on', 'whisker',2);
    xlabel('Day-Hour');
    ylabel(ylbl);
    titlename([upper(stnm) ' Diurnal Boxplots - ' ylbl]);
    if ( ~isempty(ylm) )
      ylim(ylm);
    end;
    maxigraph;
    grid on;
    hold on;
    if ( doPrint(4) )
      print('-dtiff',fullfile(figspath,sprintf('%s-%s-diurnal-boxplot.tiff',stnm,fld)));
    end;
  end;


return;
