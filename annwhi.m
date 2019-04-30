1;

% Whole-record Percentiles of various bleaching sea temperature criteria
%FWYF1 30.4 = 97%, 31.4 = 99.8%
%MLRF1 30.4 = 97%, 31.4 = 99.9%
%SMKF1 30.9 = 96%, 31.9 = 99.4%
%SANF1 30.4 = 95%, 31.4 = 99.7%
%SRVI2 29.4 = 82%, 30.4 = 98.6%
%CMRC3 29.8 = 93%, 30.8 = 99.4%

% Sort sites by Longitude (works since Papahanaumokuakea is diagonal on a map!)
[ig,sortix] = sort([sites.lon]);
sites = sites(sortix);

% Evaluate various bleaching criteria: "MMM" is taken as the 95th percentile
% (See data from other stations in Florida and the Caribbean, above.)
% Do "MMM+1" last: that is the one we want to continue analyzing

for cbleaching_criterion = {'99.6','98.6','MMM+1'}

  bleaching_criterion = cbleaching_criterion{:};

  for ix=1:length(sites)
    sites(ix).bleaching_misst_mmm = prctile(sites(ix).misst.data,95);

    switch (bleaching_criterion),
     case 'MMM+1',
      sites(ix).bleaching_misst = sites(ix).bleaching_misst_mmm + 1;
     case '99.6',
      sites(ix).bleaching_misst = prctile(sites(ix).misst.data,99.6);
     case '98.6',
      sites(ix).bleaching_misst = prctile(sites(ix).misst.data,98.6);
     otherwise,
      error('Unknown bleaching criterion "%s"',bleaching_criterion);
    end;

    sites(ix).bleaching_ix = find(sites(ix).misst.data >= sites(ix).bleaching_misst);
    sites(ix).bleaching_date = sites(ix).misst.date(sites(ix).bleaching_ix);
    sites(ix).bleaching_sri = sites(ix).misst.data(sites(ix).bleaching_ix)-sites(ix).bleaching_misst;
  end;

  X={};
  Y={};
  ylbl={};
  for ix =1:length(sites)
    X{end+1} = sites(ix).misst.date;
    Y{end+1} = sites(ix).misst.data;
    ylbl{end+1} = sprintf('%s (#%g Z=%gm)',sites(ix).name,ix,sites(ix).max_depth);
  end;
  [hln,hax] = multiplot_datetick(X,Y,'MISST and Bleaching EF',[],ylbl);
  switch (bleaching_criterion),
   case 'MMM+1',
    appendtitlename(' (MMM+1^oC)');
   case '99.6',
    appendtitlename(' (99.6%)');
   case '98.6',
    appendtitlename(' (98.6%)');
  end;
  for axix=1:length(hax)
    axes(hax(axix));
    yl = get(get(gca,'YLabel'),'String');
    ix=[];
    for ixix=1:length(sites)
      if ( strncmpi(yl,sites(ixix).name,length(sites(ixix).name)) )
        ix=ixix;
      end;
    end;
    if ( isempty(ix) ); error('No time series to match YLabel "%s"',yl); end;
    hold on;
    plot(sites(ix).bleaching_date, sites(ix).misst.data(sites(ix).bleaching_ix),'rs');
    text(sites(ix).bleaching.date,19,[num2str(sites(ix).bleaching.data) '% bleached'],'Color','red');
    set(hax(axix),'FontSize',8);
    set(get(hax(axix),'YLabel'),'FontSize',8);
  end;

end;


disp('Hit enter when ready to generate time series'); pause;

fhraw=figure; maxigraph; plot_ts(sites.misst); legend(sites.name);
hold on; plot(sites(ix).bleaching_date,sites(ix).bleaching_sri,'rs');
titlename('MISST 9km Daily Product');

sites(1).misst_30_day_average=[];
for ix=1:length(sites)
  warning('off','BuildDerivedVar:EmptyField');
  sites(ix) = verify_variable(sites(ix),'misst_30_day_average');
  warning('on','BuildDerivedVar:EmptyField');
end;
fh30dma=figure; maxigraph; plot_ts(sites.misst_30_day_average); legend(sites.name);
hold on; plot(sites(ix).bleaching_date,sites(ix).bleaching_sri,'rs');
titlename('MISST 30-day Moving Average');

sites(1).misst_30_day_lowpass=[];
for ix=1:length(sites)
  warning('off','BuildDerivedVar:EmptyField');
  sites(ix) = verify_variable(sites(ix),'misst_30_day_lowpass');
  warning('on','BuildDerivedVar:EmptyField');
end;
fh30dlp=figure; maxigraph; plot_ts(sites.misst_30_day_lowpass); legend(sites.name);
hold on; plot(sites(ix).bleaching_date,sites(ix).bleaching_sri,'rs');
titlename('MISST 30-day Lowpass Filter');


disp('Hit enter when ready to review all time series'); pause;

for fh=[fh30dma,fh30dlp,fhraw];
  figure(fh);
  for yr=2002:2010;
    axis([datenum(yr,8,10) datenum(yr,10,10) 27.59 29.11]);
    datetick3('x',2,'keeplimits');
    pause(1);
  end;
  pause(3);
end;
