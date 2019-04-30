function stations = fix_pf(stnms)
%function stations = fix_pf(stnms)
%
%error('THIS WAS A RUN-ONCE FIX-IT FUNCTION - DO NOT RERUN THIS!');
error('THIS WAS A RUN-ONCE FIX-IT FUNCTION - DO NOT RERUN THIS!');

  datapath = get_ecoforecasts_path('data');

  % Variable name prefix: Path Finder v5.0 Cloud-filtered 5-day mean
  PFX = 'pfv5c_pentad_';

  for ix=1:numel(stnms)
    matfname = fullfile(datapath,[lower(stnms{ix}) '_pathfinder_pentad.mat']);
    % If we already did this before - do not need to load again.
    if ( ~exist(matfname,'file') )
      error('Whoops! Missing "%s"',matfname);
    else
      disp(['Reloading MAT file ' matfname]);
      load(matfname,'station');
      stations(ix) = orderfields(station);
      station = []; clear station;
    end;
  end;
  % for ix=1:numel(stations)
  %   stations(ix).pfv5c_pentad_interp_method = stations(ix).pfv5c_pentad__interp_method;
  % end;
  % stations = rmfield(stations,'pfv5c_pentad__interp_method');


  % yrad = 10;
  % xrad = 11;
  yrad = 4;
  xrad = 5;


  url = 'http://data.nodc.noaa.gov/thredds/dodsC/pathfinder/Version5.0_CloudScreened/5day/FullRes/2009/2009161-2009165.pfv5sst_5day.hdf';

  nc = mDataset(url);
  try
    lons = cast(nc{'Longitude'}(:),'double');
    lats = cast(nc{'Latitude'}(:),'double');

    for ix = 1:length(stations);
      [yerr,yix(ix)] = min(abs(lats - stations(ix).lat));
      [xerr,xix(ix)] = min(abs(lons - stations(ix).lon));
      if ( yerr > max(diff(unique(lats))) || xerr > max(diff(unique(lons))) )
        % Global product - this should never happen!
        error('Ecoforecasts:Pathfinder:StationOutsideDomain',...
              'Station "%s" at %g,%g is outside Pathfinder domain?!',...
              stations(ix).station_name,stations(ix).lat,stations(ix).lon);
      else
        lat{ix} = lats(yix(ix)-2:yix(ix)+2);
        lon{ix} = lons(xix(ix)-3:xix(ix)+3);
      end;
    end;
  catch
    close(nc); clear nc
    rethrow(lasterror);
  end;
  close(nc); clear nc


  for ix=1:length(stations)
    stations(ix).([PFX 'lonix']) = xix(ix);
    stations(ix).([PFX 'latix']) = yix(ix);

    stations(ix).([PFX 'sst_field']).lon = lons(xix(ix)-xrad:xix(ix)+xrad);
    stations(ix).([PFX 'sst_field']).lat = lats(yix(ix)-yrad:yix(ix)+yrad);

    % Make sure lon and lat are row vectors as for other products
    stations(ix).([PFX 'sst_field']).lon = stations(ix).([PFX 'sst_field']).lon(:);
    stations(ix).([PFX 'sst_field']).lat = stations(ix).([PFX 'sst_field']).lat(:);

    stations(ix).([PFX 'sst_field']).field = stations(ix).([PFX 'sst_field']).field(:,7:end-6,7:end-6);
    matfname = fullfile(datapath,[lower(stations(ix).station_name) '_pathfinder_pentad.mat']);
    disp(['Saving MAT file ' matfname]);
    station = stations(ix);
    save(matfname,'station');
    station=[]; clear station;
  end;

return;
