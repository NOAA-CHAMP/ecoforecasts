function station_monthly_boxplot(stn,fld,ttl,xlbl,ylbl,ylm,doyrs,domos,tol,doPrint)
%function station_monthly_boxplot(stn,fld,ttl,xlbl,ylbl,ylm,doyrs,domos,tol,doPrint)
% 
% Produce BOXPLOT (qv.) for monthly (seasonal cycle) frequency distributions
% of time series STN.(FLD).  If optional DOYRS or DOMOS are given, include
% only data from those years and/or months in the distributions. If TOL is
% given (DEFAULT: 0.75), only use year-months having at least that fraction
% of total possible periodic values, e.g., if TOL==0.75, remove year-months
% with fewer than 0.75*28*24 == 504 samples. For non-hourly periodic data
% (e.g., daily means), sample period is assumed to be constant, equal to the
% shortest gap between successive values: MIN(DIFF(UNIQUE(STN.(FLD).date))).
% 
% Last Saved Time-stamp: <Wed 2011-04-20 14:06:37  Lew.Gramer>


  figspath = get_ecoforecasts_path('figs');

  stnm = lower(stn.station_name);

  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = strrep(fld,'_','\_');
  elseif ( strcmpi(ylbl,'none') )
    ylbl = [];
  end;
  if ( ~exist('ttl','var') || isempty(ttl) )
    ttl = [upper(stnm) ' Monthly Boxplots - ' ylbl];
  elseif ( strcmpi(ttl,'none') )
    ttl = [];
  end;
  if ( ~exist('xlbl','var') || isempty(xlbl) )
    xlbl = 'Year-Month';
  elseif ( strcmpi(xlbl,'none') )
    xlbl = [];
  end;

  if ( ~exist('ylm','var') || isempty(ylm) )
    ylm = [];
  end;
  if ( ~exist('doyrs','var') || isempty(doyrs) )
    doyrs = [];
  end;
  if ( ~exist('domos','var') || isempty(domos) )
    domos = [];
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    tol = 0.75;  % QC tolerance = 75% of expected # hourly values
  end;
  if ( ~exist('doPrint','var') || isempty(doPrint) )
    doPrint = false;
  end;


  samples_per_day = 1/min(diff(unique(stn.(fld).date)));

  % Filter time series dates based on user request
  [yrs,mos,dys,hrs,ig,ig] = datevec(stn.(fld).date);
  dat = stn.(fld).data;

  if ( isempty(doyrs) )
    doyrs = unique(yrs);
  else
    doyrix = find(ismember(yrs,doyrs));
    yrs = yrs(doyrix);
    mos = mos(doyrix);
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
  dys = dys(goodix);
  hrs = hrs(goodix);
  dat = dat(goodix);

  if ( isempty(dat) )
    error('No data left after quality control!');
  end;


  % Do box-and-whisker plots on monthly distributions

  % NOTE: Unlike STATION_BOXPLOTS, this function assumes a figure already
  % exists, or else it allows BOXPLOT to create a default one here...

  boxplot(dat, mos, 'notch','on', 'whisker',1.5);
  xlabel(xlbl);
  ylabel(ylbl);
  xlim([0 13]);
  if ( ~isempty(ylm) )
    ylim(ylm);
  end;
  if ( ~isempty(ttl) )
    titlename(ttl);
  end;

  if ( doPrint )
    figfile = fullfile(figspath,sprintf('%s-%s-seasonal-boxplot.tiff',stnm,fld));
    disp(['Printing to ' figfile]);
    print('-dtiff',figfile);
  end;

return;
