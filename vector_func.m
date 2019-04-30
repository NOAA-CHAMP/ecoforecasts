function stn = vector_func(stn, inst, comp, curr)
%function stn = vector_func(stn, inst, comp, curr)
%
% Calculate a vector component ('u' or 'v') based on the "speed" and "dir"
% fields in the given instrument "inst" in station struct "stn". If optional
% fourth arg "curr" is TRUE (1), components are calculated assuming the "dir"
% sensor gives "target direction" (as for ocean currents) rather than a
% "source direction" (as for winds).
%
% As of 2018 Feb 25, will also build derived time series for each colummn in
% field STN.(varname).prof, if found. If so, CALLS: FILTER_GAPS_PROF (v.)
%
% Last Saved Time-stamp: <Sun 2018-02-25 13:32:44 Eastern Standard Time gramer>

  if ( ~exist('curr', 'var') )
    curr = false;
  end;

  switch (inst)
   case {'wxt'},
    % Vaisala WXT sensors have different names in ICON/G2
    sp = [inst '_wspeed'];
    dr = [inst '_wdir'];
   otherwise,
    sp = [inst '_speed'];
    dr = [inst '_dir'];
  end;

  if ( ~isfield(stn, sp) || ~isfield(stn.(sp), 'data') || ...
       ~isnumeric(stn.(sp).data) )
    error('No valid wind speed field "%s" in stn!', sp);
  end;
  if ( ~isfield(stn, dr) || ~isfield(stn.(dr), 'data') || ...
       ~isnumeric(stn.(dr).data) )
    error('No valid wind direction field "%s" in stn!', dr);
  end;

  doProf = false;
  if ( isfield(stn.(sp),'prof') && size(stn.(sp).prof,1) == numel(stn.(sp).date) && ...
       isfield(stn.(dr),'prof') && size(stn.(dr).prof,1) == numel(stn.(dr).date) )
    doProf = true;
  end;


  % NOTE: This is a strict logical intersection: do not do fuzzy intersection
  % with INTERSECT_DATES, because we want EXACTLY coincident data. (Or do we?)
  dts = intersect(stn.(sp).date, stn.(dr).date);
  spidx = find(ismember(stn.(sp).date, dts));
  dridx = find(ismember(stn.(dr).date, dts));

  vecvar = [inst '_' lower(comp)];
  switch (lower(comp)),
   case 'u',
    stn.(vecvar).date = dts;
    stn.(vecvar).data = stn.(sp).data(spidx) .* (-sind(stn.(dr).data(dridx)));
    if ( doProf )
      stn.(vecvar).prof = stn.(sp).prof(spidx,:) .* (-sind(stn.(dr).prof(dridx,:)));
    end;
   case 'v',
    stn.(vecvar).date = dts;
    stn.(vecvar).data = stn.(sp).data(spidx) .* (-cosd(stn.(dr).data(dridx)));
    if ( doProf )
      stn.(vecvar).prof = stn.(sp).prof(spidx,:) .* (-cosd(stn.(dr).prof(dridx,:)));
    end;
   otherwise,
    error('Unrecognized vector component "%s": I know "u" and "v"!', comp);
  end;

  % Currents are e.g. eastWARD (toward E).; winds are e.g. westERLY (from W)
  if ( curr )
    stn.(vecvar).data = -stn.(vecvar).data;
    if ( doProf )
      stn.(vecvar).prof = -stn.(vecvar).prof;
    end;
  end;

  % Ensure that we don't derive nonsense from gaps in the raw data
  if ( doProf )
    stn = filter_gaps_prof(stn, sp, vecvar);
    stn = filter_gaps_prof(stn, dr, vecvar);
  else
    stn = filter_gaps(stn, sp, vecvar);
    stn = filter_gaps(stn, dr, vecvar);
  end;

return;
