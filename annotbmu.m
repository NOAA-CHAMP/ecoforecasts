function fh = annotbmu(station, sm, n, framedts, fldnms)
%function fh = annotbmu(station, sm, n, framedts, fldnms)
%
% Produce a multiplot display showing Best Matching Unit(s) 'n' from a Self
% Organizing Map (SOM) analysis, stacked with the raw data upon which the SOM
% was based: allows easy cross-identification of raw-data features (both in
% the atmospheric forcing and in hydrographic "response" variables) that have
% elicited a given set of SOM patterns. 'n' is an integer scalar or vector.
%
% Last Saved Time-stamp: <Sat 2013-07-20 17:25:48 Eastern Daylight Time gramer>

  if ( ~exist('n','var') || isempty(n) )
    frameix = 1:numel(sm.bmus);
  else
    frameix = find(ismember(sm.bmus, n));
  end;

  if ( ~exist('fldnms', 'var') || isempty(fldnms) )
    idx = 0;
    idx=idx+1; fldnms{idx} = 'wind1_u';
    idx=idx+1; fldnms{idx} = 'wind1_v';
    idx=idx+1; fldnms{idx} = 'air_t';
    idx=idx+1; fldnms{idx} = 'sea_t';
  end;

  color_order = 'bgrcmykw';
  fldnmstr = sprintf('%s(%c)\n', fldnms{end}, color_order(length(fldnms)));
  for fidx = (length(fldnms)-1):-1:1
    fldnmstr = sprintf('%s vs. %s(%c)', fldnmstr, fldnms{fidx}, color_order(fidx));
  end;

  for fidx = 1:length(fldnms)
    multidts{fidx} = station.(fldnms{fidx}).date;
    multidat{fidx} = station.(fldnms{fidx}).data;
  end;
  % Bottom plot is nday-long horizontal lines at appropriate BMU (e.g., 1-9)
  x = framedts(frameix,:)';
  multidts{end+1} = x(:)';
  x = repmat(sm.bmus(frameix), [1 size(framedts,2)])';
  multidat{end+1} = x(:)';

  [ign, axs] = ...
      multiplot( multidts, multidat, ...
                 'YLabel', strrep(lower(fldnms),'_','\_'), ...
                 'LineSpec', {'b-','g-','r-','c-', 'r.'}, ...
                 'Title', ...
                 sprintf('%s %s', upper(station.station_name), ...
                         strrep(lower(fldnmstr),'_','\_')) );
  fh = gcf;
  set(fh, 'units', 'normalized', 'outerposition', [0 0 1 1]);
  datetick2('keepticks', 'keeplimits'); set_datetick_cursor; drawnow;
  set(axs(end), 'YLim', [(min(sm.bmus(frameix))-1) (max(sm.bmus(frameix))+1)]);

return;
