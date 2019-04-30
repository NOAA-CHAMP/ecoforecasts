function plot_diurnal_range_vs_control(tses,ctlprms,peracc,persub,tsnms,fldnm,ctlprmnm)
%function plot_diurnal_range_vs_control(tses,ctlprms,peracc,persub,tsnms,fldnm,ctlprmnm)
%
% Subset time series in STRUCT array TSES with PERSUB (DEFAULT: @TS_ISFINITE)
% and plot their within-PERACC (DEFAULT: @GET_JDAY) range vs. the vector of
% control variable values CTLPRMS. NUMEL(TSES) must == NUMEL(CTLPRMS). E.g.,
% one common useage would be to plot summer diurnal sea temperature range for
% multiple sites with names TSNMS vs. their seafloor slope (CTLPRMNM 'Beta').
%
% Last Saved Time-stamp: <Sun 2016-10-30 18:49:26 Eastern Daylight Time gramer>

  if ( ~exist('tses','var') || ~isstruct(tses) || numel(tses) < 2 || ~is_ts(tses(1)) )
    error('First arg must be a STRUCT array with at least two TS elements');
  end;
  if ( ~exist('ctlprms','var') || ~isnumeric(ctlprms) || numel(ctlprms) ~= numel(tses) )
    error('Second arg must be a numeric vector of length NUMEL(TSES)');
  end;

  if ( ~exist('tsnms','var') || isempty(tsnms) )
    tsnms = inputname(1);
    if ( isempty(tsnms) )
      tsnms = 'TS';
    end;
    tsnms = cellstr( strcat(tsnms,num2str([1:numel(tses)]')) );
  elseif ( ischar(tsnms) )
    tsnms = cellstr(tsnms);
  end;
  if ( numel(tsnms) ~= numel(tses) )
    error('Optional TSNMS should be a CELLSTR of length NUMEL(TSES)');
  end;
  for tix = 1:numel(tsnms)
    tsnms{tix} = strrep(tsnms{tix},'_','\_');
  end;

  if ( ~exist('ctlprmnm','var') || isempty(ctlprmnm) )
    ctlprmnm = inputname(2);
    if ( isempty(ctlprmnm) )
      ctlprmnm = 'Control';
    end;
  end;
  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'Variable';
  end;

  if ( ~exist('peracc','var') || isempty(peracc) )
    peracc = @get_jday;
    %peracc = @get_year;
  end;
  stracc = regexprep(char(peracc),'^[@]*get_','');

  %ts1 = stns{1}.(fld);

  if ( ~exist('persub','var') || isempty(persub) )
    persub = @ts_isfinite;
    %persub = @ts_jas;
    %persub = @ts_jfm;
    %persub = @ts_boreal_warm;
    %persub = @ts_boreal_cool;
  end;
  switch ( char(persub) ),
   case 'ts_isfinite',          strsub='all';
   case 'ts_jas',               strsub='summer';
   case 'ts_jfm',               strsub='winter';
   case 'ts_boreal_warm',       strsub='rainy';
   case 'ts_boreal_cool',       strsub='dry';
   otherwise,                   strsub = char(persub);
  end;

  % Subset all TSes to the desired period
  for tix = 1:numel(tses);
    tses(tix) = subset_ts(tses(tix),persub);
  end;

  % % Subset all TSES to one another
  % tses = intersect_tses(tses);

  %pers = 1:365; % Julian Days, e.g.
  pers = unique(peracc(tses(1).date));

  sz = [numel(pers),numel(tses)];
  mn = repmat(nan,sz);
  mx = repmat(nan,sz);
  mnd = repmat(nan,sz);
  mxd = repmat(nan,sz);
  p07d = repmat(nan,sz);
  p93d = repmat(nan,sz);

  for stix=1:numel(tses)
    for ix=1:numel(pers);
      per = pers(ix);
      perix = find(peracc(tses(stix).date)==per);
      if ( isempty(perix) )
        dat = nan;
      else
        dat = tses(stix).data(perix);
      end;
      mn(ix,stix) = min(dat);
      mx(ix,stix) = max(dat);
      
      mnd(ix,stix) = min(dat) - mean(dat);
      mxd(ix,stix) = max(dat) - mean(dat);
      
      p07d(ix,stix)  = prctile(dat, 7) - mean(dat);
      p93d(ix,stix) = prctile(dat,93) - mean(dat);
    end;
  end;
  
  [ig,prmix] = sort(ctlprms,2,'descend');
  
  fmg;
  lhs = []; legs = {};
  %%peroff = 5e4;
  %%peroff = 10e4;
  %peroff = 5e3;
  %peroff = numel(pers)*3e2;
  peroff = numel(pers)*15;
  clr = {'k','r',[0.0,0.5,0.0],'b','m'};
  mrk = {'.','o','x','+','*','s','d','p','^','v'};
  for ix=prmix(:)'; 
    clrix = 1 + mod(ix-1,numel(clr));
    mrkix = 1 + floor(ix/numel(clr));
    % Left-hand (7th percentile) cluster
    lhs(end+1) = plot(p07d(:,ix),ctlprms(ix)+([1:numel(pers)]./peroff),mrk{mrkix},'Color',clr{clrix});
    % Right-hand (93rd percentile) cluster
    plot(p93d(:,ix),ctlprms(ix)+([1:numel(pers)]./peroff),mrk{mrkix},'Color',clr{clrix});
    legs{end+1} = upper(tsnms{ix});
  end;
  for ix=prmix(:)'; 
    text(0,ctlprms(ix)+(1./peroff),[upper(stracc),'=',num2str(pers(1))]);
    text(0,ctlprms(ix)+(numel(pers)./peroff),[upper(stracc),'=',num2str(pers(end))]);
  end;
  legend(lhs,legs,'Location','West');
  titlename([upper(strrep(fldnm,'_','\_')),' ',upper(strsub),' ',upper(stracc)]);
  xlabel(strrep(upper(fldnm),'_','\_'));
  ylabel(strrep(upper(ctlprmnm),'_','\_'));
  if ( doPrint )
    print('-dpng',get_relative_path(lower([ctlprmnm,'_vs_',fldnm,'_',strsub,'_',stracc,'.png'])));
  end;

return;
