function stns = load_misst_cfg(stns,region)
%function stns = load_misst_cfg(stns,region)
%
% Load MISST station config file for REGION, and extract gridpoint indices
% corresponding to each station named in STNS(:).name. New fields LONFLD and
% LATFLD are added or replaced in returned struct STNS. (LONFLD, e.g., is
% ['misst_' REGION '_lonix'], or 'misst_lonix' for the 'world' region.)
%
% NOTE: Field order in config file is STATION_LONG_NAME,[LATINDS],[LONINDS].
%
% Last Saved Time-stamp: <Wed 2012-07-18 13:58:55  lew.gramer>

  cfgpath = get_ecoforecasts_path;
  datapath = get_ecoforecasts_path('data');

  region = lower(region);
  switch ( region )
   case 'asam',
    cfgbasename = ['misst_' region '.cfg'];
   case 'freef',
    cfgbasename = ['misst_' region '.cfg'];
   case 'world',
    cfgbasename = 'misst.cfg';
   otherwise,
    error('Region "%s" not yet implemented!', region);
  end;

  if ( strcmpi(region,'world') )
    lonfld = 'misst_lonix';
    latfld = 'misst_latix';
  else
    lonfld = ['misst_' region '_lonix'];
    latfld = ['misst_' region '_latix'];
  end;

  cfgfname = fullfile(cfgpath,cfgbasename);
  fid = fopen(cfgfname,'r');
  %flds = textscan(fid, '%[^,],%[^,],%[^,\n]\n');
  flds = textscan(fid, '%[^,],%f,%f,%f,%f', 'CommentStyle','#');
  fclose(fid);

  if ( numel(flds) ~= 5 )
    error('Bad format in config file "%s"!', cfgfname);
  end;

  % % Sanity check - SHOULD GO AWAY after this function is stable
  % for ix = 1:length(flds{1})
  %   stnm = flds{1}{ix};
  %   stnix = find(strcmpi(stnm,{stns.name}));
  %   if ( isempty(stnix) )
  %     warning('MISST:ConfigUnknownStation', ...
  %             'Config file %s contained unknown station name "%s"!', ...
  %             cfgfname, stnm);
  %   end;
  % end;

  for ix = 1:numel(stns)
    cfgix = find(strcmpi(flds{1},stns(ix).name));
    if ( isempty(cfgix) )
      warning('MISST:StationNotInConfig', ...
              'Config file %s contained nothing for station "%s"!', ...
              cfgfname, stns(ix).name);

      stns(ix).(lonfld) = [];
      stns(ix).(latfld) = [];
    else
      % % Latitude index column appears first in file for historical reasons
      % stns(ix).(latfld) = eval(flds{2}{cfgix});
      % stns(ix).(lonfld) = eval(flds{3}{cfgix});
      % if ( isempty(stns(ix).(lonfld)) || ...
      %      any(size(stns(ix).(lonfld)) ~= size(stns(ix).(latfld))) )
      %   warning('MISST:BadConfigIndices', ...
      %           'Bad indices "%s,%s" loaded for station "%s"!', ...
      %           flds{2}{cfgix}, flds{3}{cfgix}, stns(ix).name);
      % end;

      % Latitude index column appears first in file for historical reasons
      stns(ix).(latfld) = [flds{2}(cfgix):flds{3}(cfgix)];
      stns(ix).(lonfld) = [flds{4}(cfgix):flds{5}(cfgix)];
      if ( isempty(stns(ix).(lonfld)) || isempty(stns(ix).(latfld)) || ...
           any(stns(ix).(lonfld)<1) || any(stns(ix).(latfld)<1) )
        warning('MISST:BadConfigIndices','Bad indices loaded for station "%s"!',stns(ix).name);
      end;
      if ( strcmpi(region,'world') )
        % Do this for all??
        stns(ix).(latfld) = stns(ix).(latfld) + 1;
        stns(ix).(lonfld) = stns(ix).(lonfld) + 1;

        % My poor dyslexic brain is hurting... Freagin' grids, freagin' MATLAB
        stns(ix).(lonfld) = stns(ix).(lonfld) - 2048;
        stns(ix).(lonfld)(stns(ix).(lonfld) < 0) = stns(ix).(lonfld)(stns(ix).(lonfld) < 0) + 4096;
      end;
    end;
  end;

return;
