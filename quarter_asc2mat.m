1;

coastpath = get_ecoforecasts_path('coast');
%fname = 'key_west_fl_mhw';
fname = 'miami_fl_mhw';

doOverwrite = false;

if ( ~exist('dat','var') )
  [lon,lat,dat,matfname,utmzone] = asc2mat(fname);
end;

  all_lat = lat;
  all_lon = lon;
  all_dat = dat;
  lat = []; lon = []; dat = []; clear lat lon dat

  latquadir = 'SN';
  lonquadir = 'WE';
  latix = [1,ceil(numel(all_lat)/2),floor(numel(all_lat)/2),numel(all_lat)];
  lonix = [1,ceil(numel(all_lon)/2),floor(numel(all_lon)/2),numel(all_lon)];
  for dlat=1:2
    lat = all_lat(latix((dlat*2)-1):latix((dlat*2)));
    for dlon=1:2
      lon = all_lon(lonix((dlon*2)-1):lonix((dlon*2)));
      dat = all_dat(latix((dlat*2)-1):latix((dlat*2)), lonix((dlon*2)-1):lonix((dlon*2)));

      quadir = [latquadir(dlat),lonquadir(dlon)];
      matfbasename = sprintf('%s_%s.mat',fname,quadir);
      matfname = fullfile(coastpath,matfbasename);
      if ( exist(matfname,'file') && ~doOverwrite )
        error('DID NOT OVERWRITE existing %s',matfname);
      end;
      disp(['Saving ',matfname]);
      save(matfname,'lat','lon','dat','fname','-v7.3');
      disp({matfbasename,min(lon),max(lon),min(lat),max(lat)});
      lon = []; dat = []; clear lon dat
    end;
    lat = []; clear lat
  end;

  lat = all_lat;
  lon = all_lon;
  dat = all_dat;
  all_lat = []; all_lon = []; all_dat = []; clear all_lat all_lon all_dat

  clear ans dlat dlon latix lonix latquadir lonquadir matfbasename matfname quadir
