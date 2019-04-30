function [P,ANOVATAB,STATS,COMP,MNS,FH] = station_anova_multcompare(stn,fld,ixs,ylbl,ylm,flgs,args)
%function [P,ANOVATAB,STATS,COMP,MNS,FH] = station_anova_multcompare(stn,fld,ixs,ylbl,ylm,flgs,args)
% 
% Perform one-way analysis of variance (v. ANOVA1) and do MULTCOMPARE (v.)
% for weekly, monthly, and season (annual cycle), and annual (interannual
% var.) frequency distributions from time series STN.(FLD). Optional arg IXS,
% vector or FUNCTION_HANDLE, specifies which indices in STN.(FLD) to compare.
% If FUNCTION_HANDLE, IXS must accept .date/.data struct and return indices.
% Optional FLGS is a logical 4-vector, with switches to turn off graphs of
% resp., annual, seasonl, monthly, weekly ANOVA. DEFAULT: [1,1,1,1]. The
% optional ARGS is a cell array of name-value pairs passed to MULTCOMPARE.
% 
% P,ANOVATAB,STATS and COMP,MNS,FH are *structures*, containing one each of
% the first three return values from ANOVA1 and MULTCOMPARE, resp. So for
% example, P.monthly contains the returned p-values from ANOVA1 for monthly
% distributions, and so on for P.weekly, P.seasonal, and P.annual.
% 
% Last Saved Time-stamp: <Thu 2015-04-02 20:19:08 Eastern Daylight Time gramer>

  if ( ~exist('ixs','var') || isempty(ixs) )
    ixs = 1:length(stn.(fld).date);
  end;
  if ( isa(ixs,'function_handle') )
    ixs = ixs(stn.(fld));
  elseif ( ~isnumeric(ixs) && ~islogical(ixs) )
    error('IXS must either be a logical or index vector, or function_handle!');
  end;

  if ( ~exist('ylbl','var') )
    ylbl = strrep(fld,'_','\_');
  end;

  if ( ~exist('flgs','var') || isempty(flgs) )
    flgs = 1;
  end;
  if ( isscalar(flgs) )
    flgs = repmat(flgs,[4,1]);
  end;

  if ( ~exist('args','var') )
    args = {};
  end;

  stnm = stn.station_name;
  dts = stn.(fld).date(ixs);
  dat = stn.(fld).data(ixs);

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;

  fldstr = strrep(fld,'_','\_');


  if (flgs(1))
    [P.annual,ANOVATAB.annual,STATS.annual] = anova1(dat,yrs);
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Interannual)']);

    %% IN PLACE OF ANOVA1, maybe try AOCTOOL?
    %     AOCTOOL(X,Y,G,ALPHA,XNAME,YNAME,GNAME) allows you to enter the
    %     names of the X, Y, and G variables.
    %
    %     [H,ATAB,CTAB,STATS] = AOCTOOL(...) returns several items.  H is a
    %     vector of handles to the line objects in the plot.  ATAB and CTAB
    %     are an anova table and a coefficients table.  STATS is a structure
    %     containing statistics useful for performing a multiple comparison
    %     of means with the MULTCOMPARE function.

    [COMP.annual,MNS.annual,FH.annual] = multcompare(STATS.annual, 'ctype','tukey-kramer dunn-sidak',args{:});
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Interannual)']);
  end;


  if (flgs(2))
    season = jds;
    season( 60 >= jds | jds >  335) = 1;
    season( 60 <  jds & jds <= 150) = 2;
    season(150 <  jds & jds <= 241) = 3;
    season(241 <  jds & jds <= 335) = 4;
    [P.seasonal,ANOVATAB.seasonal,STATS.seasonal] = anova1(dat,season);
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Seasonal - DJF,MAM,JJA,SON)']);
    [COMP.seasonal,MNS.seasonal,FH.seasonal] = multcompare(STATS.seasonal, 'ctype','tukey-kramer dunn-sidak',args{:});
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Seasonal - DJF,MAM,JJA,SON)']);
  end;


  if (flgs(3))
    [P.monthly,ANOVATAB.monthly,STATS.monthly] = anova1(dat,mos);
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Monthly)']);
    [COMP.monthly,MNS.monthly,FH.monthly] = multcompare(STATS.monthly, 'ctype','tukey-kramer dunn-sidak',args{:});
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Monthly)']);
  end;


  if (flgs(4))
    wks = ceil(jds/7);
    wks(wks > 52) = 52;
    [P.weekly,ANOVATAB.weekly,STATS.weekly] = anova1(dat,wks);
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Weekly)']);
    [COMP.weekly,MNS.weekly,FH.weekly] = multcompare(STATS.weekly, 'ctype','tukey-kramer dunn-sidak',args{:});
    maxigraph;
    appendtitlename([' (' stnm '.' fldstr ' Weekly)']);

    % figure;
    % boxplot(dat, mos, 'notch','on', 'whisker',1.5);
    % xlabel('Year-Month');
    % title([upper(stnm) ' Monthly Boxplots ' ylbl]);
    % if ( exist('ylm','var') )
    %  ylim(ylm);
    % end;
    % maxigraph;
    % %print('-dpng', sprintf('%s-mpo624-month-%s-boxplot.png', stnm, fld));

    % figure;
    % boxplot(dat, yrs, 'notch','on', 'whisker',1.5);
    % xlabel('Year');
    % titlename([upper(stnm) ' Annual Boxplots ' ylbl]);
    % if ( exist('ylm','var') )
    %   ylim(ylm);
    % end;
    % maxigraph;
    % % print('-dpng', sprintf('%s-mpo624-year-%s-boxplot.png', stnm, fld));
  end;

return;
