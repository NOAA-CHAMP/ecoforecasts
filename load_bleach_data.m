function stns = load_bleach_data(stns)
%function stns = load_bleach_data(stns)
%
% Load coral bleaching observations (as % of total cover) for each site named
% in struct vector STNS. (Assumes file [stns(ix).name '_bleach.dat'] exists
% in DATAPATH, containing tab-separated columns "year month day bleach_pct".)
%
% Last Saved Time-stamp: <Fri 2011-04-08 10:29:56  Lew.Gramer>

  datapath = get_ecoforecasts_path('data');

  for ix = 1:length(stns)

    stns(ix).bleaching.date = [];
    stns(ix).bleaching.data = [];

    if ( ~isempty(stns(ix).name) )

      fname = fullfile(datapath, [stns(ix).name '_bleach.dat']);

      if ( ~exist(fname,'file') )
        warning('ecoforecasts:NoBleaching','No bleaching file for "%s"',stns(ix).name);

      else
        dat = load(fname);
        if ( isempty(dat) )
          warning('ecoforecasts:NoBleaching','No bleaching data for "%s"',stns(ix).name);
        else
          % 0000 GMT - the universal time stamp for day-based data
          stns(ix).bleaching.date = datenum(dat(:,1),dat(:,2),dat(:,3),0,0,0);
          stns(ix).bleaching.data = dat(:,4);
        end;
      end;

    end;

  end;

return;
