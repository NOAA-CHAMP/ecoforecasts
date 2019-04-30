function plot_data_and_factory(stn, fld, f)
%function plot_data_and_factory(stn, fld, f)
%
% Plot raw station data from the field named "fld" in the station data struct
% "stn". Annotate the plot with the value ranges specified in factory "f".
%
% Also plots histogram of raw data values, broken up by the values in "f".
%
% Last Saved Time-stamp: <Mon 2009-01-26 14:49:24 Eastern Standard Time gramer>

  if ( ~isstruct(stn) || ~isfield(stn, fld) )
    error('First arg must be a struct, second arg a field name in struct.');
  end;
  if ( ~isstruct(f) || ~isfield(f, 'fuzzies') )
    error(['Third arg "f" must be a factory - a struct having a cell-array ' ...
           'with both fuzzy names and corresponding numeric data ranges in it.']);
  end;


  figure;
  dath = plot(stn.(fld).date, stn.(fld).data, '.');
  set(dath, 'MarkerSize', 4);
  datetick;
  set_datetick_cursor;

  % Ignore the -Inf and +Inf on either end
  for idx = 1:length(f.fuzzies)
    rng = f.fuzzies{idx};
    edges(idx) = rng{2}(1);
    annotline([], rng{2}(1), rng(1), 'red');
  end;

  annotline([], rng{2}(2), 'unbelievably-high', 'red');
  edges(end+1) = rng{2}(2);

  th = title(sprintf('%s vs %s factory', fld, f.variable));
  set(th, 'Interpreter', 'none');


  figure;
  hc = histc(stn.(fld).data, edges);
  bar(edges, hc, 'histc');

return;
