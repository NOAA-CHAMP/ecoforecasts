function [P,ANOVATAB,STATS,COMP,MNS,FH] = station_subset_anova(stn,fld,ixs,subnm,ylbl)
%function [P,ANOVATAB,STATS,COMP,MNS,FH] = station_subset_anova(stn,fld,ixs,subnm,ylbl)
% 
% Perform one-way analysis of variance (v. ANOVA1) and do MULTCOMPARE (v.)
% for annual frequency distributions (i.e., interannual variability) from
% time series STN.(FLD). IXS, vector or FUNCTION_HANDLE, specifies which
% indices in STN.(FLD) to compare over all years. SUBNM names the subset.
% If FUNCTION_HANDLE, IXS must accept .date/.data struct and return indices.
% 
% P,ANOVATAB,STATS and COMP,MNS,FH are the first three return values from
% ANOVA1 and MULTCOMPARE, resp. (e.g., P contains p-values from ANOVA1).
% 
% Last Saved Time-stamp: <Sat 2010-06-12 13:43:06 Eastern Daylight Time gramer>

  if ( ~exist('ixs','var') || isempty(ixs) )
    ixs = 1:length(stn.(fld).date);
  end;
  if ( isa(ixs,'function_handle') )
    ixs = ixs(stn.(fld));
    if ( isempty(ixs) )
      error('Index function IXS returned an empty set!');
    end;
  elseif ( ~isnumeric(ixs) && ~islogical(ixs) )
    error('IXS must either be a logical or index vector, or function_handle!');
  end;

  if ( ~exist('subnm','var') || isempty(subnm) )
    subnm = '';
  end;

  if ( ~exist('ylbl','var') )
    ylbl = strrep(fld,'_','\_');
  end;

  stnm = stn.station_name;
  dts = stn.(fld).date(ixs);
  dat = stn.(fld).data(ixs);

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;


  [P,ANOVATAB,STATS] = anova1(dat,yrs,'on');
  ylabel(ylbl);
  maxigraph;
  appendtitlename([' (' stnm '.' ylbl ' ' subnm ')']);
  ylim([nanmin(dat(:)) nanmax(dat(:))]);

  %% IN PLACE OF ANOVA1, maybe try AOCTOOL?
  %     AOCTOOL(X,Y,G,ALPHA,XNAME,YNAME,GNAME) allows you to enter the
  %     names of the X, Y, and G variables.
  %
  %     [H,ATAB,CTAB,STATS] = AOCTOOL(...) returns several items.  H is a
  %     vector of handles to the line objects in the plot.  ATAB and CTAB
  %     are an anova table and a coefficients table.  STATS is a structure
  %     containing statistics useful for performing a multiple comparison
  %     of means with the MULTCOMPARE function.

  figure;
  [COMP,MNS,FH] = multcompare(STATS, ...
                              'ctype','tukey-kramer dunn-sidak', ...
                              'display','on');
  ylabel(ylbl);
  maxigraph;
  % appendtitlename([' (' stnm '.' ylbl ' Interannual)']);
  titlename(['Tukey-Kramer/Dunn-Sidak: ' stnm '.' ylbl ' ' subnm]);
  xlim([nanmin(dat(:)) nanmax(dat(:))]);


%   figure;
%   boxplot(dat, mos, 'notch','on', 'whisker',1.5);
%   xlabel('Year-Month');
%   title([upper(stnm) ' Monthly Boxplots ' ylbl]);
%   if ( exist('ylm','var') )
%     ylim(ylm);
%   end;
%   maxigraph;
%   %print('-dpng', sprintf('%s-mpo624-month-%s-boxplot.png', stnm, fld));

%   figure;
%   boxplot(dat, yrs, 'notch','on', 'whisker',1.5);
%   xlabel('Year');
%   titlename([upper(stnm) ' Annual Boxplots ' ylbl]);
%   if ( exist('ylm','var') )
%     ylim(ylm);
%   end;
%   maxigraph;
%   % print('-dpng', sprintf('%s-mpo624-year-%s-boxplot.png', stnm, fld));

return;
