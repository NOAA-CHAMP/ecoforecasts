function stn = postprocess_adcp_uv(stn,varargin)
%function stn = postprocess_adcp_uv(stn[,ufld,vfld[,spfld,drfld[,sfcfld,midfld,btmfld]]][,deps[,beam_angle]])
%
% Place NaNs near and above surface spike in each column of arrays of
% acoustic Doppler ocean current profiles in struct STN. Assumes STN has
% fields .adcp_u, .adcp_v each having .date and .prof fields. Optional DEPS
% is a vector of depths; if not specified, the u and v structs must also have
% .depths fields, and they must match. Removes an additional percentage of
% bins below surface spike for side-lobe contamination, assuming BEAM_ANGLE
% (DEFAULT: 20o) between ADCP transducer beams and the vertical. Once done,
% (re-)calculates depth-averaged velocities in STN.adcp_u/v.data. If optional
% SPFLD and DRFLD strings are specified, also recalculate speed and direction
% (both in profiles and .data) in fields STN.(SPFLD) and STN.(DRFLD), resp.
% If optional SFCFLD, MIDFLD, BTMFLD, calculate bin-averaged near-surface,
% mid-water, and near-bottom current time series using simple definitions:
% "near-bottom" is the lowest one-third of bins, "near-surface" the same # of
% bins below the "last good bin" found by despiking (above), and "mid-water"
% all bins above "near-bottom" and below the lowest "near-surface" bin.
%
% NOTE, e.g., giving non-empty SFCFLD creates time series STN.([SFCFLD,'_u'])
% and STN.([SFCFLD,'_v']); if both SPFLD and DRFLD non-empty, also creates
% speed and dir'n fields STN.([SFCFLD,'_speed']) and STN.([SFCFLD,'_dir']).
% Ditto for strings MIDFLD and BTMFLD, if they are non-empty.
%
% NOTE ALSO: If UFLD and VFLD are specified and SPFLD is logical scalar True,
% tries to figure out defaults for all the fieldnames SPFLD, DRFLD, SFCFLD,
% MIDFLD, and BTMFLD, by replacing the regular expression '_u$' in UFLD with
% '_speed', '_dir', '_sfc', '_mid', and '_btm', respectively.
%
% Last Saved Time-stamp: <Fri 2018-03-23 18:26:39 Eastern Daylight Time gramer>

  args = varargin;

  %% Handle calling arguments

  % Positional optional field names
  if ( numel(args) > 0 )
    ufld = args{1};
    args(1) = [];
  end;
  if ( numel(args) > 0 )
    vfld = args{1};
    args(1) = [];
  end;

  % Defaults for optional field names
  if ( ~exist('ufld','var') || isempty(ufld) )
    ufld = 'adcp_u';
  end;
  if ( ~exist('vfld','var') || isempty(vfld) )
    vfld = 'adcp_v';
  end;
  if ( ~isfield(stn,ufld) || ~isfield(stn,vfld) )
    error('STN must have fields .%s and .%s',ufld,vfld);
  end;

  if ( numel(args) > 0 && ischar(args{1}) )
    spfld = args{1};
    args(1) = [];
  elseif ( numel(args) > 0 && islogical(args{1}) && isscalar(args{1}) && args{1} )
    % Logical True: Figure out default SPFLD, DRFLD, SFCFLD, MIDFLD, BTMFLD
    if ( isempty(regexp(ufld,'_u$')) )
      error('Cannot figure out default fieldnames: no "_u$" in UFLD');
    else
      spfld = regexprep(ufld,'_u$','_speed');
      drfld = regexprep(ufld,'_u$','_dir');
      sfcfld = regexprep(ufld,'_u$','_sfc');
      midfld = regexprep(ufld,'_u$','_mid');
      btmfld = regexprep(ufld,'_u$','_btm');
    end;
    args(1) = [];
  end;
  if ( numel(args) > 0 && ischar(args{1}) )
    spfld = args{1};
    args(1) = [];
  end;
  if ( numel(args) > 0 && ischar(args{1}) )
    drfld = args{1};
    args(1) = [];
  end;
  if ( numel(args) > 0 && ischar(args{1}) )
    sfcfld = args{1};
    args(1) = [];
  end;
  if ( numel(args) > 0 && ischar(args{1}) )
    midfld = args{1};
    args(1) = [];
  end;
  if ( numel(args) > 0 && ischar(args{1}) )
    btmfld = args{1};
    args(1) = [];
  end;

  % Defaults for optional field names
  if ( ~exist('spfld','var') )
    spfld = [];
  end;
  if ( ~exist('drfld','var') )
    drfld = [];
  end;
  if ( ~exist('sfcfld','var') )
    sfcfld = [];
  end;
  if ( ~exist('midfld','var') )
    midfld = [];
  end;
  if ( ~exist('btmfld','var') )
    btmfld = [];
  end;

  % Vector of bin depths
  if ( numel(args) > 0 && isvector(args{1}) && isnumeric(args{1}) )
    deps = args{1};
    args(1) = [];
  end;
  if ( ~exist('deps','var') )
    if ( ~isfield(stn.(ufld),'depths') || ~isfield(stn.(vfld),'depths') ...
         || numel(stn.(ufld).depths) ~= numel(stn.(vfld).depths) ...
         || ~all(stn.(ufld).depths == stn.(vfld).depths) )
      error('Need matching .depths fields in STN.%s and STN.%s, or need DEPS vector',ufld,vfld);
    end;
    deps = stn.(ufld).depths;
  end;

  hgts = -(deps(1) - deps);

  % ADCP beam angle with vertical
  if ( numel(args) > 0 && isscalar(args{1}) )
    beam_angle = args{1};
    args(1) = [];
  end;
  if ( ~exist('beam_angle','var') )
    beam_angle = 20;
  end;
  % Beam angles 20o: side lobe contamination in last 6% + 1 bin of profile
  side_lobe_contam_frac = 1-cosd(beam_angle);


  %% Do QC post-processing

  dts = stn.(ufld).date(:);

  binthird = floor(numel(hgts)/3);

  % Use as upper bound, bin number where ostensible depths become non-negative
  sfcix = find(deps>=0,1);

  % QC
  % Loop over each timestamp, to find the bin number just below the surface spike
  stn.adcp_sfc_height.date = dts;
  stn.adcp_sfc_height.data = repmat(max(hgts),[numel(dts),1]);

  % We'll limit our spike search to the top 1/3 of the profile bins
  topthirdix = numel(hgts)-binthird+1;

  last_good_bin = repmat(numel(hgts),[numel(dts),1]);
  for ix = 1:numel(dts)
    uprof = stn.(ufld).prof(ix,:);
    vprof = stn.(vfld).prof(ix,:);
    if ( any(~isnan(uprof)) && any(~isnan(vprof)) )
      prof = uv_to_spd(uprof,vprof);

      [ig,spikeix] = nanmax(prof(topthirdix:sfcix));
      spikeix = topthirdix+spikeix-1;
      if ( ~isempty(spikeix) )
        stn.adcp_sfc_height.data(ix) = hgts(spikeix);

        % Beam angles 20o: side lobe contamination in last 6% + 1 bin of profile
        last_good_bin(ix) = spikeix - (floor(side_lobe_contam_frac*spikeix) + 1);
        %DEBUG:if (last_good_bin(ix)==topthirdix || last_good_bin(ix)==numel(hgts)); keyboard; end;
        stn.(ufld).prof(ix,last_good_bin(ix)+1:end) = nan;
        stn.(vfld).prof(ix,last_good_bin(ix)+1:end) = nan;
      end;
    end;
    stn.(ufld).data(ix,1) = nanmean(uprof);
    stn.(vfld).data(ix,1) = nanmean(vprof);
  end;


  %% Calculate bin averages: speed, direction, surface, mid, bottom

  if ( ~isempty(spfld) && ~isempty(drfld) )
    [stn.(spfld),stn.(drfld)] = ts_uv_to_spddir(stn.(ufld),stn.(vfld),false);
  end;

  btmbins = 1:binthird;

  if ( ~isempty(btmfld) )
    stn.([btmfld,'_u']).date = dts;
    stn.([btmfld,'_u']).data = nanmean(stn.(ufld).prof(:,btmbins),2);
    stn.([btmfld,'_v']).date = dts;
    stn.([btmfld,'_v']).data = nanmean(stn.(vfld).prof(:,btmbins),2);
  end;

  if ( ~isempty(midfld) )
    stn.([midfld,'_u']).date = dts;
    stn.([midfld,'_v']).date = dts;
    if ( isempty(sfcfld) )
      % If they want mid but not surface, then ignore "last_good_bin"
      midbins = btmbins(end)+1:btmbins(end)+binthird;
      stn.([midfld,'_u']).data = nanmean(stn.(ufld).prof(:,midbins),2);
      stn.([midfld,'_v']).data = nanmean(stn.(vfld).prof(:,midbins),2);
    else
      stn.([midfld,'_u']).data = repmat(nan,size(dts));
      stn.([midfld,'_v']).data = repmat(nan,size(dts));
    end;
  end;

  if ( ~isempty(sfcfld) )
    stn.([sfcfld,'_u']).date = dts;
    stn.([sfcfld,'_v']).date = dts;
    stn.([sfcfld,'_u']).data = repmat(nan,size(dts));
    stn.([sfcfld,'_v']).data = repmat(nan,size(dts));

    for ix = 1:numel(dts)
      sfcbins = (last_good_bin(ix)-binthird+1):last_good_bin(ix);
      if ( ~isempty(sfcbins) )
        stn.([sfcfld,'_u']).data(ix,1) = nanmean(stn.(ufld).prof(ix,sfcbins),2);
        stn.([sfcfld,'_v']).data(ix,1) = nanmean(stn.(vfld).prof(ix,sfcbins),2);
      end;
      if ( ~isempty(midfld) )
        midbins = (btmbins(end)+1):(sfcbins(1)-1);
        if ( ~isempty(midbins) )
          stn.([midfld,'_u']).data(ix,1) = nanmean(stn.(ufld).prof(ix,midbins),2);
          stn.([midfld,'_v']).data(ix,1) = nanmean(stn.(vfld).prof(ix,midbins),2);
        end;
      end;
    end;
  end;

  if ( ~isempty(spfld) && ~isempty(drfld) )
    if ( ~isempty(btmfld) )
      stn.([btmfld,'_speed']).date = dts;
      stn.([btmfld,'_speed']).data = uv_to_spd(stn.([btmfld,'_u']).data,stn.([btmfld,'_v']).data);
      stn.([btmfld,'_dir']).date = dts;
      stn.([btmfld,'_dir']).data = uv_to_dir_curr(stn.([btmfld,'_u']).data,stn.([btmfld,'_v']).data);
    end;
    if ( ~isempty(midfld) )
      stn.([midfld,'_speed']).date = dts;
      stn.([midfld,'_speed']).data = uv_to_spd(stn.([midfld,'_u']).data,stn.([midfld,'_v']).data);
      stn.([midfld,'_dir']).date = dts;
      stn.([midfld,'_dir']).data = uv_to_dir_curr(stn.([midfld,'_u']).data,stn.([midfld,'_v']).data);
    end;
    if ( ~isempty(sfcfld) )
      stn.([sfcfld,'_speed']).date = dts;
      stn.([sfcfld,'_speed']).data = uv_to_spd(stn.([sfcfld,'_u']).data,stn.([sfcfld,'_v']).data);
      stn.([sfcfld,'_dir']).date = dts;
      stn.([sfcfld,'_dir']).data = uv_to_dir_curr(stn.([sfcfld,'_u']).data,stn.([sfcfld,'_v']).data);
    end;
  end;

return;
