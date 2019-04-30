function fhs = plotbmus(station, sm, n, framedts, fldnms)

  if ( ~exist('fldnms', 'var') || isempty(fldnms) )
    idx = 0;
    idx=idx+1; fldnms{idx} = 'wind1_u';
    idx=idx+1; fldnms{idx} = 'wind1_v';
    idx=idx+1; fldnms{idx} = 'air_t_qc';
    idx=idx+1; fldnms{idx} = 'sea_t_qc';
  end;

  frameix = find(sm.bmus == n);
  nframes = length(frameix);

  color_order = 'bgrcmykw';
  fldnmstr = sprintf('%s(%c)\n', fldnms{end}, color_order(length(fldnms)));
  for fidx = (length(fldnms)-1):-1:1
    fldnmstr = sprintf('%s vs. %s(%c)', fldnmstr, fldnms{fidx}, color_order(fidx));
  end;

  for fidx = 1:length(fldnms)
    dts{fidx} = station.(fldnms{fidx}).date;
    dat{fidx} = station.(fldnms{fidx}).data;
  end;

  maxplots = min(5, nframes);
  for ix = (nframes-maxplots+1):nframes
    begdt = framedts(frameix(ix), 1);
    enddt = framedts(frameix(ix), 2);
    for fidx = 1:length(fldnms)
      multidts{fidx} = dts{fidx}(begdt <= dts{fidx} & dts{fidx} <= enddt);
      multidat{fidx} = dat{fidx}(begdt <= dts{fidx} & dts{fidx} <= enddt);
    end;

    multiplot( multidts, multidat, ...
               'YLabel', strrep(lower(fldnms),'_','\_'), ...
               'LineSpec', {'b-','g-','r-','c-'}, ...
               'Title', ...
               sprintf('%s %s \n Frame %d (%s-%s)', upper(station.station_name), ...
                       strrep(lower(fldnmstr),'_','\_'), ...
                       ix, datestr(begdt), datestr(enddt)) );
    fhs(ix) = gcf;
    set(fhs(ix), 'units', 'normalized', 'outerposition', [0 0 1 1]);
    datetick2('keepticks', 'keeplimits'); set_datetick_cursor; drawnow;
  end;

return;
