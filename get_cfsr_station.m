function stns = get_cfsr_station(stns)

  vars = {...
      'Downward_Short-Wave_Rad_Flux_surface',...
      'Upward_Short-Wave_Rad_Flux_surface',...
      'U-component_of_wind',...
      'V-component_of_wind',...
         };
  flds = { 'cfsr_dsrf','cfsr_usrf','cfsr_wind_u','cfsr_wind_v', };

  %nc = mDataset('http://nomads.ncdc.noaa.gov/thredds/dodsC/modeldata/cmd_flxf/2002/200207/20020703/flxf00.gdas.2002070300.grb2');
  urlbase = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/modeldata/cmd_flxf';

  % dates = datenum(2002,7,1):(6/24):datenum(2011,1,1)-(1/24);
  dates = datenum(2002,7,1):(6/24):datenum(2002,7,3)-(1/24);

  lastyr = -inf;
  for dtix = 1:length(dates)
    dt = dates(dtix);
    [yr,mo,dy,hr,mn,sc] = datevec(dt);
    % Give caller feedback
    if ( yr ~= lastyr )
      disp(yr);
      lastyr = yr;
    end;

    urlname = sprintf('%s/%04d/%04d%02d/%04d%02d%02d/flxf00.gdas.%04d%02d%02d%02d.grb2',...
                      urlbase,yr,yr,mo,yr,mo,dy,yr,mo,dy,hr);
    %DEBUG:
    disp(urlname); tic,
    nc = mDataset(urlname);
    if ( isempty(nc) )
      warning('Skipping missing date "%s"',datestr(dt));
      break;
    end;

    try  % netCDF Java Toolkit does not handle unclosed GRiB files very well
      if ( ~exist('lon','var') )
        lon = nc{'lon'}(:,:);
        lat = nc{'lat'}(:,:);
        for stix=1:length(stns)
          stns(stix).lonix=round(interp1(lon,1:length(lon),360+stns(stix).lon));
          stns(stix).latix=round(interp1(lat,1:length(lat),stns(stix).lat));
          for ix=1:length(vars)
            fld = flds{ix};
            % stns(stix).(fld) = struct('date',[],'data',[]);
            stns(stix).(fld).date = repmat(nan,[length(dates) 1]);
            stns(stix).(fld).data = repmat(nan,[length(dates) 1]);
          end;
        end;
      end;

      for ix=1:length(vars)
        ncvars{ix} = nc{vars{ix}};
      end;
      for ix=1:length(vars)
        var = vars{ix};
        fld = flds{ix};
        for stix=1:length(stns)
          stns(stix).(fld).date(dtix,1) = dt;
          % stns(stix).(fld).data(dtix,1) = nc{var}(1,1,stns(stix).latix,stns(stix).lonix);
          stns(stix).(fld).data(dtix,1) = ncvars{ix}(1,1,stns(stix).latix,stns(stix).lonix);
        end;
      end;
    catch
      close(nc); clear nc
      rethrow(lasterror);
    end;

    close(nc); clear nc
    %DEBUG:
    toc,
    % Give the NOMADS server some breathing time
    pause(0.25);

  end;

return;
