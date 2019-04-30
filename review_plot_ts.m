function review_plot_ts(pers,fh)
%function review_plot_ts(pers,fh)
%
% Review a time series plot (v. PLOT_TS) one period (year, season, month) at
% a time, by shifting the XLIM of the plot forward from [PERS(1),PERS(2)] to
% [PERS(2),PERS(3)], etc. each time the user hits the Enter key. The figure
% to shift is FH (DEFAULT: GCF at the time of the initial call.) If focus is
% changed from FH for any reason, FIGURE(FH) refocuses before each shift.
%
% Last Saved Time-stamp: <Tue 2012-04-24 12:39:23  lew.gramer>

  if ( ~exist('pers','var') || isempty(pers) )
    pers = datenum(1993:2012,1,1);
  end;

  if ( ~exist('fh','var') || isempty(fh) )
    fh = gcf;
  end;

  msg = sprintf('Proceed with period-by-period review? (Afterward, to\n move forward a period hit Enter, "Q" or Ctrl-C to stop...)');
  disp(msg);
  YOrN = questdlg(msg,upper(mfilename));

  if ( strcmpi(YOrN,'y') )
    disp('Skipped review at user request...');

  else
    for perix=2:numel(pers)
      if ( ~ishandle(fh) )
        disp('Figure closed? Reviewed halted...');
        break;
      end;
      figure(fh);
      xlim([pers(perix-1),pers(perix)]);
      datetick3;
      waitforbuttonpress;
      qOrNot = get(fh,'CurrentCharacter');
      if ( ~isempty(qOrNot) && strcmpi(qOrNot(1),'q') )
        disp('Reviewed halted at user request...');
        break;
      end;
    end;
  end;

return;
