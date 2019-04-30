function sfp = load_SFP_nutrients
%function sfp = load_SFP_nutrients
%
% Load nutrients and other cast data from 2002-2008 taken aboard the R/V
% F. G. Walton Smith quarterly to semi-annual cruises for the NOAA SFP, the
% South Florida (Everglades Restoration Monitoring) Program.
%
% Last Saved Time-stamp: <Sun 2018-03-04 18:37:41 Eastern Standard Time gramer>

  matfname = fullfile(get_ecoforecasts_path('data'),'SFP_nutrients.mat');

  if ( exist(matfname,'file') )

    disp(['Loading ',matfname]);
    sfp = load(matfname);


  else

    disp(['Extracting nutrients from master spreadsheet']);

    x = importdata(fullfile(get_ecoforecasts_path('data'),['walton_smith_database.xls']));

    goodix = find(~cellfun(@isempty,x.textdata.Sheet1(:,5)));
    goodix(1) = [];
    dts = datenum( {x.textdata.Sheet1{goodix,5}} );

    lon = x.data.Sheet1(goodix-1,1);
    lat = x.data.Sheet1(goodix-1,2);
    z   = x.data.Sheet1(goodix-1,4);
    T1  = x.data.Sheet1(goodix-1,7);
    T2  = x.data.Sheet1(goodix-1,8);
    S1  = x.data.Sheet1(goodix-1,9);
    S2  = x.data.Sheet1(goodix-1,10);
    O2  = x.data.Sheet1(goodix-1,11);
    chl  = x.data.Sheet1(goodix-1,16);
    pha  = x.data.Sheet1(goodix-1,17);
    N   = x.data.Sheet1(goodix-1,18);
    Si  = x.data.Sheet1(goodix-1,19);
    P   = x.data.Sheet1(goodix-1,20);

    ztop = 0;
    Ttop = 23;
    Tbtm = 12;
    Nbtm = 0;
    Pbtm = 0;
    Sibtm = 0;

    dtcrit=(get_year(dts)>=1997 & ismember(get_month(dts),[4:10])); %Apr-Oct
    % dtcrit=(get_year(dts)>=1997 & (get_season(dts)==2|get_season(dts)==3)); %Apr-Sep

    selix = find(dtcrit & (T1<Ttop) & (T1>=Tbtm) & (z>ztop) & (N>=Nbtm) & (P>=Pbtm) & (Si>=Sibtm));

    % [B,Stats] = scatter_fit(23-T1(selix),N(selix),'23-T','NO3+NO2'); %,'none');
    % [B,Stats] = scatter_fit(23-T1(selix),P(selix),'23-T','P'); %,'none');
    % [B,Stats] = scatter_fit(23-T1(selix),Si(selix),'23-T','Si'); %,'none');

    sfp.longitude.date = dts(selix);
    sfp.longitude.data = lon(selix);

    sfp.latitude.date = dts(selix);
    sfp.latitude.data = lat(selix);

    sfp.depth.date = dts(selix);
    sfp.depth.data = z(selix);

    sfp.seatemp.date = dts(selix);
    sfp.seatemp.data = T1(selix);

    sfp.salinity.date = dts(selix);
    sfp.salinity.data = S1(selix);

    sfp.oxygen2.date = dts(selix);
    sfp.oxygen2.data = O2(selix);

    sfp.chlora.date = dts(selix);
    sfp.chlora.data = chl(selix);

    sfp.NOx.date = dts(selix);
    sfp.NOx.data = N(selix);

    sfp.P.date = dts(selix);
    sfp.P.data = P(selix);

    sfp.Si.date = dts(selix);
    sfp.Si.data = Si(selix);

    disp(['Saving ',matfname]);
    save(matfname,'-struct','sfp');

  end; %if ( exist(matfname,'file') ) else

return;
