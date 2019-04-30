function sites = nwhi_misst(sites)
%function sites = nwhi_misst(sites)
%
% Load North West Hawaiian Islands (Papahanaumokuakea Monument) coral trasect
% site coordinates and bleaching observations, and MISST (Multi-instrument
% Improved Sea Surface Temperature) v03 median daily SST

  datapath = get_ecoforecasts_path('data');

  matfname = fullfile(datapath,'nwhi_misst.mat');

  if ( ~exist('sites','var') || isempty(sites) )

    if ( exist(matfname,'file') )
      disp(['Loading ' matfname]);
      load(matfname,'sites');

    else

      sites_matfname = fullfile(datapath,'NWHI_sites.mat');
      if ( exist(sites_matfname,'file') )
        load(sites_matfname,'sites');
      else
        %DEBUG:
error('Oops?!');

        x = importdata('C:\Documents and Settings\lew.gramer\My Documents\coral\Bleach\NWHI\LIST_OF_BLE_SITES - from Bernardo.xls');
        nsites = length(x.textdata(4:end,1));
        for ix=1:nsites
          sites(ix).name = x.textdata{ix+3,1};
          sites(ix).islandname = x.textdata{ix+3,3};
          sites(ix).type = x.textdata{ix+3,4};
          sites(ix).min_depth = x.data(ix,1);
          sites(ix).max_depth = x.data(ix,2);
          sites(ix).lon = x.data(ix,3);
          sites(ix).lat = x.data(ix,4);
          sites(ix).misst = struct('date',[],'data',[]);
        end;
        x=[]; clear x
        save(sites_matfname,'sites');
      end; %if ( exist(sites_matfname,'file') ) else


      nsites = length(sites);

      % Determine coordinate indices for each site in Papa region
      [lon,lat,sst] = read_misst_region('papa',2004,228);
      for ix=1:nsites
        sites(ix).papa_lonix = round(interp1(lon,1:length(lon),sites(ix).lon));
        sites(ix).papa_latix = round(interp1(lat,1:length(lat),sites(ix).lat));
      end;

      % Determine coordinate indices for each site in Papa region
      [lon,lat,sst] = read_misst_region('world',2009,228);
      for ix=1:nsites
        sites(ix).global_lonix = round(interp1(lon,1:length(lon),sites(ix).lon));
        sites(ix).global_latix = round(interp1(lat,1:length(lat),sites(ix).lat));
      end;

      % Manual adjustments (per site) to MISST pixel registration
      sites(1).lonadj = [0];        sites(1).latadj = [-1];
      sites(2).lonadj = [0 +1];     sites(2).latadj = [0 +1];
      sites(3).lonadj = [0];        sites(3).latadj = [0];
      sites(4).lonadj = [0];        sites(4).latadj = [0];
      sites(5).lonadj = [0];        sites(5).latadj = [-1 0];
      sites(6).lonadj = [0];        sites(6).latadj = [0 +1];
      sites(7).lonadj = [0];        sites(7).latadj = [-1 0];
      sites(8).lonadj = [-1 0];     sites(8).latadj = [0 +1];
      sites(9).lonadj = [0];        sites(9).latadj = [-1 0];

      for ix=1:nsites
        sites(ix).papa_lonix = sites(ix).papa_lonix + sites(ix).lonadj;
        sites(ix).papa_latix = sites(ix).papa_latix + sites(ix).latadj;
        sites(ix).global_lonix = sites(ix).global_lonix + sites(ix).lonadj;
        sites(ix).global_latix = sites(ix).global_latix + sites(ix).latadj;
      end;


      % Extract MISST v03 data for each site, from whatever MISST binary files are available!

      papayrs = [2002:2006];
      globyrs = [2006:2011];

      % summerdys = 140:278;  dys = summerdys;
      alldys = 1:365;

      disp('Region PAPA MISST files');
      for yr=papayrs(:)'
        disp(yr);
        switch (yr)
         case 2002,  dys = 184:365;
         case 2006,  dys = 1:206;
         otherwise,  dys = alldys;
        end;
        for jd=dys(:)'
          [lon,lat,sst] = read_misst_region('papa',yr,jd);
          if ( ~isempty(sst) )
            dt = datenum(yr,1,1) + jd - 1;
            for ix=1:length(sites)
              sites(ix).misst.date(end+1) = dt;
              if ( all(isfield(sites(ix),{'papa_latix','papa_lonix'})) )
                site_sst = sst(sites(ix).papa_latix,sites(ix).papa_lonix);
                sites(ix).misst.data(end+1) = nanmedian(site_sst(:));
              else
                sites(ix).misst.data(end+1) = ...
                    interp2(lon,lat,sst,sites(ix).lon,sites(ix).lat,'*cubic');
              end;
            end;
          end;
        end;
      end;

      disp('Global MISST files');
      [yrnow,mo,dy] = datevec(now);
      jdnow = floor(now) - datenum(yrnow,1,1) + 1;
      for yr=globyrs(:)'
        disp(yr);
        switch (yr)
         case 2006,  dys = 234:365;
         case yrnow, dys = 1:(jdnow-45);   %Seemingly normal delay for v03 files
         otherwise,  dys = alldys;
        end;
        for jd=dys(:)'
          [lon,lat,sst] = read_misst_region('world',yr,jd);
          if ( ~isempty(sst) )
            dt = datenum(yr,1,1) + jd - 1;
            for ix=1:length(sites)
              sites(ix).misst.date(end+1) = dt;
              if ( all(isfield(sites(ix),{'global_latix','global_lonix'})) )
                site_sst = sst(sites(ix).global_latix,sites(ix).global_lonix);
                sites(ix).misst.data(end+1) = nanmedian(site_sst(:));
              else
                sites(ix).misst.data(end+1) = ...
                    interp2(lon,lat,sst,sites(ix).lon,sites(ix).lat,'*linear');
              end;
            end;
          end;
        end;  
      end;


      disp('Loading actual 2004 bleaching observations');
      sites(1).bleaching.date = datenum(2004,9,10);
      sites(1).bleaching.data = 0; %FFS-32
      sites(2).bleaching.date = datenum(2004,9,10);
      sites(2).bleaching.data = 0; %FFS-R30

      sites(3).bleaching.date = datenum(2004,9,10);
      sites(3).bleaching.data = 56.2; %PHR-31
      sites(4).bleaching.date = datenum(2004,9,10);
      sites(4).bleaching.data = 74.6; %PHR-32

      sites(5).bleaching.date = datenum(2004,9,10);
      sites(5).bleaching.data = 13.4; %LIS-R7

      sites(6).bleaching.date = datenum(2004,9,10);
      sites(6).bleaching.data = 6.5; %KUR-17
      sites(7).bleaching.date = datenum(2004,9,10);
      sites(7).bleaching.data = 28.4; %KUR-R36
      sites(8).bleaching.date = datenum(2004,9,10);
      sites(8).bleaching.data = 1.9; %KUR-14

      sites(9).bleaching.date = datenum(2004,9,10);
      sites(9).bleaching.data = 35.1; %LIS-10


      disp(['Saving ' matfname]);
      save(matfname,'sites');

    end; %if ( exist(matfname,'file') ) else

  end; %if ( ~exist('sites','var') || isempty(sites) )


  nsites = length(sites);

  % Sort sites by Longitude (works since Papahanaumokuakea is diagonal on a map!)
  [ig,sortix] = sort([sites.lon]);
  sites = sites(sortix);

  % Do simple multiplot of MISST time series for all sites
  X={};
  Y={};
  yls={};
  for ix = 1:nsites
    X{end+1} = sites(ix).misst.date;
    Y{end+1} = sites(ix).misst.data;
    yls{end+1} = sites(ix).name;
  end;
  multiplot_datetick(X,Y,'NWHI MISST Time Series',[],yls,[],[15 30]);

return;
