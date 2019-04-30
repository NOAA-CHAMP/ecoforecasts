function [stns,bath] = get_langdon_thermistors(existing_stns,bath)
%function [stns,bath] = get_langdon_thermistors(existing_stns,bath)
%
% Return a struct with all of the FKNMS thermistor data provided by Chris Langdon
%
% Last Saved Time-stamp: <Wed 2017-08-02 11:37:42 Eastern Daylight Time gramer>

  set_more off;

  stns = [];

  if ( ~exist('existing_stns','var') )
    existing_stns = [];
  end;

  if ( ~exist('bath','var') )
    bath = [];
  end;

  datapath = get_ecoforecasts_path('data');

  matfname = fullfile(datapath,'langdon_thermistors.mat');
  if ( exist(matfname,'file') )
    disp(['Load ',matfname]);
    load(matfname);

  else
    diary(fullfile(datapath,[mfilename,'.log']));

    disp('Extracting individual CSV files...');

    langpath = fullfile(datapath,'langdon');

    % Process metadata that was precompiled by Lew.Gramer@noaa.gov
    metafnm = fullfile(langpath,'metadata.xlsx');
    x = importdata(metafnm);
    ixen = 2:size(x.textdata,1);
    stnms = x.textdata(ixen,1);
    lgnms = x.textdata(ixen,2);
    latst = x.textdata(ixen,3);
    lonst = x.textdata(ixen,4);
    % Calculated depth in [m]
    deps = x.data(:,2);

    % Parse weird 'DDD-MM.ff' coordinate-string format
    latc = double(split(latst,'-'));
    lats = latc(:,1)+(latc(:,2)/60);

    lonc = double(split(lonst,'-'));
    % And convert to longitude West
    lons = -(lonc(:,1)+(lonc(:,2)/60));

    % Sanity-check the metadata
    if ( numel(stnms) ~= numel(lgnms) || numel(stnms) ~= numel(lons) || numel(stnms) ~= numel(lats) )
      error('Metadata file has inconsistent names/lons/lats: %s',metafnm);
    end;


    % Get depth and slope according to bathymetric data
    % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
    bath.sites = find_ngdc_slope_sites({lons,lats},bath,3); 

    % The 30 m F010 bathymetry has a rectangular "hole" around Key West
    nanix = find(isnan(bath.sites.depths) | isnan(bath.sites.betas));
    if ( ~isempty(nanix) )
      % Fill that hole in F010 with higher-resolution bathymetry
      hibath = read_hires_bathymetry('sanf1',[40e3,40e3]);
      hidepths = interp2(hibath.ngdc_hires_bathy.lon,hibath.ngdc_hires_bathy.lat,...
                           hibath.ngdc_hires_bathy.field,lons,lats,'*linear',nan);
      % Bathymetry resolution near Key West is 10 m: 9 pts. ~ 90 m
      [hibetas,hibeta_angs,hiiso_angs,hibath.ngdc_hires_bathy] = find_ngdc_slope(hibath.ngdc_hires_bathy,lons,lats,9);

      bath.sites.depths(nanix) = hidepths(nanix);
      bath.sites.betas(nanix) = hibetas(nanix);
      bath.sites.beta_angs(nanix) = hibeta_angs(nanix);
      bath.sites.iso_angs(nanix) = hiiso_angs(nanix);
    end; %if ( any(isnan(ngdc_deps)) )


    % Process each data directory for which we had metadata
    dfiles = dir(fullfile(langpath,'*','data','0-data','FKNMS_*.csv'));

    for ix = 1:numel(dfiles)

      fnm = fullfile(dfiles(ix).folder,dfiles(ix).name);
      stnm = regexprep(dfiles(ix).name,'(FKNMS_.*)_WQDATA.*','$1');
      %disp(stnm);

      stnix = find(strcmp(stnms,stnm));
      if ( isempty(stnix) )
        warning('Unknown station for stub name: %s (filename: %s)',stnm,fnm);

      else
        disp(['Processing ',fnm]);
        fid = fopen(fnm,'r');
        C = textscan(fid,'%*s%s%s%*s%*s%*s%*s%*s%*q%f%*s%*s','Delimiter',',','HeaderLines',1);
        fclose(fid);
        dtstr = strcat(C{1},{' '},C{2});
        dts = datenum(dtstr);
        if ( numel(dts) ~= numel(C{3}) )
          warning('Mismatched sizes (format error): %s',fnm);
        else
          fld = regexprep(stnm,'[&-]','_');
          if ( ~isfield(stns,fld) )
            stns.(fld).station_name = stnms{stnix};
            stns.(fld).station_long_name = lgnms{stnix};
            stns.(fld).lon = lons(stnix);
            stns.(fld).lat = lats(stnix);
            stns.(fld).depth = deps(stnix);
            stns.(fld).ngdc_depth = -bath.sites.depths(stnix);
            stns.(fld).ngdc_beta = bath.sites.betas(stnix);
            stns.(fld).ngdc_beta_ang = bath.sites.beta_angs(stnix);
            stns.(fld).ngdc_iso_ang = bath.sites.iso_angs(stnix);

            stns.(fld).fknms_seatemp.date = dts;
            stns.(fld).fknms_seatemp.data = C{3};

          else
            res.fknms_seatemp.date = dts;
            res.fknms_seatemp.data = C{3};
            w = warning('OFF','Ecoforecasts:mergedNonTS');
            stns.(fld) = merge_station_data(stns.(fld),res);
            warning(w);
            res=[]; clear res

          end; %if ( ~isfield(stns,fld) ) else

        end; %if ( numel(dts) ~= numel(C{3}) ) else

      end; %if ( isempty(stnix) ) else

    end; %for ix = 1:numel(dfiles)

    disp(['Save ',matfname]);
    save(matfname,'stns','-v7.3');

    diary off;

  end; %if ( exist(matfname,'file') ) else

  if ( ~isempty(existing_stns) )
    w = warning('OFF','Ecoforecasts:mergedNonTS');
    new_stns = stns;
    stns = existing_stns;
    for cf = fieldnames(new_stns)'
      fld = cf{:};
      if ( isfield(stns,fld) )
        stns.(fld) = merge_station_data(stns.(fld),new_stns.(fld));
      end;
    end;
    warning(w);
  end;

  set_more;

return;
