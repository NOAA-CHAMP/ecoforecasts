function stns = load_clean_as_thermistor(stns)
%function stns = load_clean_as_thermistor(stns)

  datapath = get_ecoforecasts_path('data');

  for ix = 1:length(stns)

    if ( isempty(stns(ix).fname) )
      dts = [];
      sst = [];
    else
      fname = [stns(ix).fname '.cdp'];
      [dts,sst] = load_cdp(fname);
    end;

    % Eliminate weird repeated dates
    nreps = 0;
    badix = find(diff(dts) <= 0);
    while ( ~isempty(badix) )
      dts(badix+1) = [];
      sst(badix+1) = [];
      nreps = nreps + length(badix);
      badix = find(diff(dts) <= 0);
    end;
    if ( nreps > 0 )
      warning('Eliminated %d repeat dates in %s', nreps, fname);
    end;

    stns(ix).sea_t.date = dts;
    stns(ix).sea_t.data = sst;
    if ( isempty(stns(ix).sea_t.date) )
      warning('No valid data found for %s (fname "%s")!', ...
              stns(ix).name, stns(ix).fname);
    end;

  end;

return;
