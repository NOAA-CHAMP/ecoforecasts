function STATIONS = get_all_station_metadata()
%function STATIONS = get_all_station_metadata()
%
% Return a struct STATIONS with four fields: 'codes' a cell array of upper
% case 5-char station codes; and 'lons','lats','depths' - station coordinates.
%
% *Replaces* GET_SEAKEYS/GET_SEAKEYS_ICON.
%
% Last Saved Time-stamp: <Fri 2018-11-23 14:52:34 Eastern Standard Time gramer>

  persistent ALL_STATIONS;

  % Later: add AIMS, AIMS-GBROOS, NMFS-CRED, CRW, etc.
  % OR download/read 'stations.txt' file from ICON data-integration server

  if ( ~isstruct('ALL_STATIONS') )

    ALL_STATIONS.codes = ...
      { ...
          'CMRC3', ...
          'SRVI2', ...
          'LPPR1', ...
          'DBJM1', ...
          'LCIY2', ...
          '41140', ...          % NDBC buoy, Christiansted, St. Croix
          'LLBP7', ...          % Laolao Bay, Saipan
          'PVGF1', ...          % Port Everglades channel, Florida: 26.092 N 80.109 W: 2018-Mar-02, http://www.ndbc.noaa.gov/station_page.php?station=pvgf1: pvgf1,26.092,-80.109
          'BBVI1', ...		% Brewers Bay - St. Thomas - USVI - PLANNED - BBVI1,18.3449,-64.9843
          'POMF1', ...		% Abs. turbidity cal/val site neardof Port of Miami channel - pomf1,25.74897,-80.13317
          'MANM1', ...		% manm1,15.24,145.715: Managaha - Saipan - CNMI
          'SASP7', ...		% sasp7,14.122,145.159: Sasanhaya Bay - Rota - CNMI
          'TASP7', ...		% tasp7,14.106,145.192: Talakhaya/Sabana Watershed - Rota - CNMI
          'TINP7', ...		% tinp7,14.9573,145.6194: Tinian Harbor (reefs) - Tinian - CNMI
          'PGBP7', ...		% pgbp7,13.428056,144.796389: NOAA NOS Station PGBP7 - 1631428 - Pago Bay, Guam
          'OOUH1', ...          % oouh1 - 1612340 - Honolulu, HI - 21.303333,-157.864444,-1.5
          '51001', ...		% nwbh1,24.416667,-162.100000 - NDBC Station 51001 (LLNR 28006) - NW HI 1 - 170 NM WNW of Kauai
          ...
          'PPDR1', ...		% ppdr1,19.833,-70.731: CREWS Dominican Republic - Puerto Plata,data from 1998-01-01 to 2016-08-25
          'CWDR1', ...		% cwdr1,18.432,-69.580: CREWS - Dominican Republic - Catuan Wreck, data from 1998-01-01 to 2016-08-25
          'DRBB1', ...		% drbb1,13.184,-59.646: CREWS - Barbados - Dottin's Reef, data from 1998-01-01 to 2016-08-25
          'HCBZ1', ...		% hcbz1,17.194325,-87.521850: Half-Moon Caye - Belize
          'CCBZ1', ...		% ccbz1,17.273567,-87.810000: Calabash Caye - Belize
          'SWBZ1', ...		% swbz1,16.815650,-88.076867: South Water Caye - Belize
          'BUTO1', ...		% buto1,11.1760,-60.8338: Buccoo Marine Park - Tobago
          'ARTO1', ...		% arto1,11.3009,-60.5206: Angel Reef - Tobago
          'MHAB1', ...		% mhab1,17.0155544,-61.8722687 # Monk's Head - Antigua and Barbuda
          'PRSK1',...		% prsk1,17.3563538,-62.8548927 # Paradise Reef - St. Kitts
          'SMSL1',...		% smsl1,13.859093,-61.071307 # Devil's Hole - Soufriere Marine Management Association (SMMA) - St. Lucia
          'KRGN1',...		% krgn1,12.022750,-61.790983 # Kahonae Reef - Grenada
          'PISR1', ...		% pisr1,8.7251,81.2080: Pigeon Island National Park - Sri Lanka - Indian Ocean
          ...
          'LKWF1', ...
          'VAKF1', ...
          'FWYF1', ...
          'CRYF1', ...
          'MLRF1', ...
          'AQUA1', ...
          'LONF1', ...
          'NFBF1', ...
          'TNRF1', ...
          'PKYF1', ...
          'VCAF1', ...
          'KYWF1', ...
          'MOSE1', ...
          'SMKF1', ...
          'LOOE1', ...
          'AMSF1', ...
          'SANF1', ...
          'PLSF1', ...
          'DRYF1', ...
          'GKYF1', ...          % gkyf1,24.627,-82.872
          '42023', ...          % C13 - West Florida South Buoy
          '42003', ...
          'SIPF1', ...          % Sebastian Inlet State Park, FL: 27.862 N 80.445 W (27°51'42" N 80°26'41" W
          ...
          ... % NCORE experiment 2007-2008, courtesy TN Lee and S Sponaugle, by K Shulzitski in 2011.
          'KLGF1', ... % NCORE: klgf1,25.0313833,-80.3480167,23 # Key Largo Current/Temperature 7/23 m (offshore French Reef)
          'MRTF1', ... % NCORE: mrtf1,24.7403000,-80.7763667,23 # Marathon Current/Temperature 8/23 m (offshore Tennessee Reef)
          ... % NCORE studies 2000-2002, provided courtesy of TN Lee, S Sponaugle, and E Williams
          'NCORA', ... % NCORE: ncora,25.1090,-80.3803,4.1
          'NCORB', ... % NCORE: ncorb,25.0932,-80.3550,7.0
          'NCORC', ... % NCORE: ncorc,25.0673,-80.3183,21.8
          'NCOR1', ... % NCORE: ncor1,25.0740,-80.3178,26.0
          'NCOR2', ... % NCORE: ncor2,25.0733,-80.3183,21.8
          'NCOT1', ... % NCORE: ncot1,25.06900,-80.31950,9.75
          'NCOT2', ... % NCORE: ncot2,25.07383,-80.32450,7.0
          'NCOT3', ... % NCORE: ncot3,25.07950,-80.33433,4.0
          'NCOT4', ... % NCORE: ncot4,25.08783,-80.34717,5.6
          ...
          'TAVRK', ...  % Lirman: Tavernier Rocks
          'CONSH', ...  % Lirman: Shallow Conch
          'CONDP', ...  % Lirman: Deep Conch
          'BNPIN', ...  % Lirman: BNP Inshore
          'BNPMI', ...  % Lirman: BNP Mid-channel
          'BNPPA', ...  % Lirman: Pacific Reef
          'BNPON', ...  % Lirman: Old Nursery
          'BNPNN', ...  % Lirman: New Nursery
          ...
          'CCAF1', ...  	% CCAF1: Coral Restoration Fnd'n / Omics Keys Disease: Carysfort
          'CNDF1', ...  	% CNDF1: Coral Restoration Fnd'n / Omics Keys Disease: North Dry Rocks
          'CGRF1', ...  	% CGRF1: Coral Restoration Fnd'n / Omics Keys Disease: Grecian Rocks
          'CPIF1', ...  	% CPIF1: Coral Restoration Fnd'n / Omics Keys Disease: Pickles Reef
          ...
          'DRBF1', ...		% drbf1,26.4619167,-80.0420833 # Delray Beach (South Central) outfall
          'BOCF1', ...		% bocf1,26.3502667,-80.0540500 # Boca Raton outfall
          'BRWF1', ...		% brwf1,26.2513833,-80.0620667 # Broward outfall
          'HOLF1', ...		% holf1,26.0191167,-80.0859333 # Hollywood outfall
          'MINF1', ...		% minf1,25.9200500,-80.0862667 # Miami North outfall
          'MICF1', ...		% micf1,25.7428167,-80.0859667 # Miami Central outfall
          'HO2F1', ...		% ho2f1,26.01915,-80.10863 # Current sensor mount near Hollywood outfall
          ...
          'AOAT_CHEECA_ROCKS_3', ...    % Atlantic Ocean Acidification testbed prospective site
          'AOAT_BROAD_KEY_2', ...       % Atlantic Ocean Acidification testbed prospective site
          ...
          'AOAT_CHEECA_MAPCO2', ... % AOAT_CHEECA_MAPCO2: Atlantic Ocean Acidification testbed ACTUAL site
          ... % Cheeca metadata source: http://mercury.ornl.gov/ocean/send/xsltText2?fileURL=/data/Mercury_instances/ocean/oceanuw/harvested/mercury-ops2.ornl.gov_OceanOME_admin_OceanMetadata_underway_Mooring_Cheeca_80W_25N_Dec2011_Dec2012.xml&full_datasource=OCEAN%20Underway&full_queryString=dataset_idText%20:%20%22Mooring_Cheeca_80W_25N%22%20AND%20(%20datasource%20:(%20oceanuw%20%20)%20)%20&ds_id=
          'FWC_MK3', ...		% Karen Neely @FWC: Marker 3: 25.3734 -80.16071667
          'FWC_BB', ...			% Karen Neely @FWC: Ball Buoy: 25.3162 -80.18806667
          ...
          'SECREMP_MC2', ...    % Brian Walker @NSUOC: MC2 Martin County Inner Shelf
          'SECREMP_PB1', ...    % Brian Walker @NSUOC: PB1 Palm Beach Inner Reef
          'SECREMP_BC1', ...    % Brian Walker @NSUOC: BC1 Broward County Inner Reef
          'SECREMP_DC1', ...    % Brian Walker @NSUOC: DC1 (Miami-)Dade County Inner Reef
          'SECREMP_UPDB', ...	% Brian Walker @NSUOC: UPDB Upside Down Barge	27 13.946	80 06.704	19.2
          'SECREMP_SLIB', ...	% Brian Walker @NSUOC: SLIB St. Lucie Inlet Barge 27.200233 -80.111583 19.80
          'SECREMP_EVCA', ...	% Brian Walker @NSUOC: EVCA Evans-Crary 27.15572 -80.05632
          'SECREMP_PELA', ...	% Brian Walker @NSUOC: PELA Peck's Lake 27.12212 -80.12803
          'SECREMP_PB2', ...	% Brian Walker @NSUOC: PB2 Outer Reef	26 40.725	80 01.099	16.7
          'SECREMP_BC3', ...	% Brian Walker @NSUOC: BC3 Outer Reef	26 09.518	80 04.641	16.7
          'SECREMP_DC3', ...	% Brian Walker @NSUOC: DC3 Outer Reef	25 50.526	80 05.286	16.7
          ...
          'SECREMP_18HO', ...	% Luke McEachron @FWC/FWRI: 18HO - 18th Hole
          'SECREMP_BLRO', ...	% Luke McEachron @FWC/FWRI: BLRO - Blowing Rocks
          'SECREMP_BSBA', ...	% Luke McEachron @FWC/FWRI: BSBA - Bullshark Barge
          'SECREMP_CEBA', ...	% Luke McEachron @FWC/FWRI: CEBA - Cement Barge
          'SECREMP_CPRE', ...	% Luke McEachron @FWC/FWRI: CPRE - Clifton Perry Reef
          'SECREMP_DBNO', ...	% Luke McEachron @FWC/FWRI: DBNO - DB North
          'SECREMP_FPCC', ...	% Luke McEachron @FWC/FWRI: FPCC,Civic Center Reef
          'SECREMP_FPLR', ...	% Luke McEachron @FWC/FWRI: FPLR,FPL Reef
          'SECREMP_HORF', ...	% Luke McEachron @FWC/FWRI: HORF,House of Refuge (Valentine)
          'SECREMP_LYBR', ...	% Luke McEachron @FWC/FWRI: LYBR,Lyon Bridge
          'SECREMP_MIBA', ...	% Luke McEachron @FWC/FWRI: MIBA,Middle Bar
          'SECREMP_NORO', ...	% Luke McEachron @FWC/FWRI: NORO,Nolens Rock
          'SECREMP_PIBA', ...	% Luke McEachron @FWC/FWRI: PIBA,SL Pipe Barge
          'SECREMP_SNSO', ...	% Luke McEachron @FWC/FWRI: SNSO,SS-1
          'SECREMP_TERE', ...	% Luke McEachron @FWC/FWRI: TERE,Texas Reef
          'SECREMP_THSO', ...	% Luke McEachron @FWC/FWRI: THSO,Three Holes (South)
          'SECREMP_TUCA', ...	% Luke McEachron @FWC/FWRI: TUCA,Tunnels/Caves
          'SECREMP_NOULA', ...	% Luke McEachron @FWC/FWRI: Noula Express artificial reef, his depth 36.6
          'SECREMP_ALPHA', ...	% Luke McEachron @FWC/FWRI: Alpha artificial reef, his depth 30.5
          'SECREMP_RODEO', ...	% Luke McEachron @FWC/FWRI: Rodeo artificial reef, his depth 30.5
          ...
          'SFOCB', ...		% Alex Soloviev @NSUOC: SFOCB: 'c-buoy',[26.0728,-80.0878, 20, 2]
          'SFOEB', ...		% Alex Soloviev @NSUOC: SFOEB: 'e-buoy',[26.0695,-80.0768, 50, 3]
          'SFONE', ...		% Alex Soloviev @NSUOC: SFONE: 'ne-buoy',[26.0695,-80.0768, 50, 3]
          'SFONW', ...		% Alex Soloviev @NSUOC: SFONW: 'nw_w-btm',[26.0705,-80.0942, 11, 10]
          'SFOSE', ...		% Alex Soloviev @NSUOC: SFOSE: 'se-buoy',[26.1398,-80.0533, 145, 8]
          'SFOSW', ...		% Alex Soloviev @NSUOC: SFOSW: 'sw-buoy',[26.0327,-80.0918, 20, 1]
          'SFOPI', ...		% Alex Soloviev @NSUOC: SFOPI: 'pier-cc',[26.0577,-80.1090, 0, 4]
          ...
          '41114', ...		% NDBC buoy Fort Pierce, FL (134)
          '41009', ...		% (LLNR 840) - CANAVERAL 20 NM East of Cape Canaveral "CANF1"
          '41113', ...		%Cape Canaveral Nearshore, FL (143) "CNNF1"
          '41012', ...		%Station 41012 (LLNR 845.3) - 40NM ENE of St Augustine FL
          ...
          'SBFG1', ...		% Stetson Bank - Flower Garden Banks NMS
          'EBFG1', ...		% East Flower Garden Banks - FGBNMS
          'WBFG1', ...		% West Flower Garden Banks - FGBNMS
          'SOFG1', ...		% Sonnier Banks - Flower Garden Banks NMS
          ...
          'SPGF1', ...		%spgf1,26.70444,-78.99472 # Station SPGF1 - Settlement Point, GBI, Bahamas
          ...
          'OFSP6', ...		%ofsp6,-14.184164,-169.658181: American Samoa - Ofu - Craig et al study sites
          'ONSP6', ...		%onsp6,-14.164692,-169.627813: American Samoa - Ofu - north coast
          'PGSP6', ...		%pgsp6,-14.296,-170.668: American Samoa - Tutuila - just outside Pago Pago harbor
          'VESP6', ...		%vesp6,-14.235,-170.661: American Samoa - Tutuila - just east of Vatia MPA
          'VWSP6', ...		%vwsp6,-14.252,-170.710: American Samoa - Tutuila - just west of Vatia MPA
          'AOSP6', ...		%aosp6,-14.251,-170.590: American Samoa - Tutuila - NNW of Aoa village
          ...
          'PBBE1', ...		%pbbe1,32.29,-64.975 # Pilchard Bay - Bermuda (see also BERM1) - 5 m
          'SWBE1', ...		%swbe1,32.23,-64.863 # Southwest Breaker Area - Bermuda (see also BERM1) - 6 m
          'NRBE1', ...		%nrbe1,32.475,-64.77 # North Rock - Bermuda (see also BERM1) - 2 m
      };
    ALL_STATIONS.lons = ...
      [ ...
          -76.139, ...          % CMRC3
          -64.761, ...          % SRVI2
          -67.052, ...          % LPPR1
          -77.420, ...          % DBJM1
          -80.060, ...          % LCIY2
          -64.723, ...          % 41140 - NDBC buoy, Christiansted, St. Croix
          145.769620, ...       % LLBP7
          -80.109, ...          % PVGF1 - Port Everglades channel, Florida
          -64.9843, ...		% Brewers Bay - St. Thomas - USVI - PLANNED - BBVI1,18.3449,-64.9843
          -80.13317, ...	% Abs. turbidity cal/val site neardof Port of Miami channel - pomf1,25.74897,-80.13317
          145.715, ...		% manm1,15.24,145.715: Managaha - Saipan - CNMI
          145.159, ...		% sasp7,14.122,145.159: Sasanhaya Bay - Rota - CNMI
          145.192, ...		% tasp7,14.106,145.192: Talakhaya/Sabana Watershed - Rota - CNMI
          145.6194, ...		% tinp7,14.9573,145.6194: Tinian Harbor (reefs) - Tinian - CNMI
          144.796389, ...	% pgbp7,13.428056,144.796389: NOAA NOS Station PGBP7 - 1631428 - Pago Bay, Guam
          -157.864444, ...      % oouh1 - 1612340 - Honolulu, HI - 21.303333,-157.864444,-1.5
          -162.100000, ...	% nwbh1,24.416667,-162.100000 - NDBC Station 51001 (LLNR 28006) - NW HI 1 - 170 NM WNW of Kauai
          ...
          -70.731, ...		% ppdr1,19.833,-70.731: CREWS Dominican Republic - Puerto Plata,data from 1998-01-01 to 2016-08-25
          -69.580, ...		% cwdr1,18.432,-69.580: CREWS - Dominican Republic - Catuan Wreck, data from 1998-01-01 to 2016-08-25
          -59.646, ...		% drbb1,13.184,-59.646: CREWS - Barbados - Dottin's Reef, data from 1998-01-01 to 2016-08-25
          -87.521850, ...	% hcbz1,17.194325,-87.521850: Half-Moon Caye - Belize
          -87.810000, ...	% ccbz1,17.273567,-87.810000: Calabash Caye - Belize
          -88.076867, ...	% swbz1,16.815650,-88.076867: South Water Caye - Belize
          -60.8338, ...		% buto1,11.1760,-60.8338: Buccoo Marine Park - Tobago
          -60.5206, ...		% arto1,11.3009,-60.5206: Angel Reef - Tobago
          -61.8722687, ...	% mhab1,17.0155544,-61.8722687 # Monk's Head - Antigua and Barbuda
          -62.8548927,...	% prsk1,17.3563538,-62.8548927 # Paradise Reef - St. Kitts
          -61.071307,...	% smsl1,13.859093,-61.071307 # Devil's Hole - Soufriere Marine Management Association (SMMA) - St. Lucia
          -61.790983,...	% krgn1,12.022750,-61.790983 # Kahonae Reef - Grenada
          81.2080, ...		% pisr1,8.7251,81.2080: Pigeon Island National Park - Sri Lanka - Indian Ocean
          ...
          -80.033, ...          % LKWF1
          -80.161667, ...       % VAKF1
          -80.096811, ...       % FWYF1
          -80.212, ...          % CRYF1
          -80.376412, ...       % MLRF1
          -80.45833, ...        % AQUA1
          -80.862, ...          % LONF1
          -81.095833, ...       % NFBF1
          -80.782375, ...       % TNRF1
          -80.7461111, ...      % PKYF1
          -81.1050000, ...      % VCAF1
          -81.8083333, ...      % KYWF1
          -81.169333, ...       % MOSE1
          -81.110818, ...       % SMKF1
          -81.402166, ...       % LOOE1
          -81.520, ...          % AMSF1
          -81.880, ...          % SANF1
          -82.773333, ...       % PLSF1
          -82.861666, ...       % DRYF1
          -82.872, ...          % gkyf1,24.627,-82.872
          -83.073611, ...       % 42023 - C13 - West Florida South Buoy
          -85.594, ...          % 42003 - NDBC buoy near edge of Loop Current
          -80.445, ...          % Sebastian Inlet State Park, FL: 27.862 N 80.445 W (27°51'42" N 80°26'41" W
          ...
          -80.3480167, ... % NCORE: klgf1,25.0313833,-80.3480167,23 # Key Largo Current/Temperature 7/23 m (offshore French Reef)
          -80.7763667, ... % NCORE: mrtf1,24.7403000,-80.7763667,23 # Marathon Current/Temperature 8/23 m (offshore Tennessee Reef)
          ...
          -80.3803, ... % NCORE: ncora,25.1090,-80.3803,4.1
          -80.3550, ... % NCORE: ncorb,25.0932,-80.3550,7.0
          -80.3183, ... % NCORE: ncorc,25.0673,-80.3183,21.8
          -80.3178, ... % NCORE: ncor1,25.0740,-80.3178,26.0
          -80.3183, ... % NCORE: ncor2,25.0733,-80.3183,21.8
          -80.31950, ... % NCORE: ncot1,25.06900,-80.31950,9.75
          -80.32450, ... % NCORE: ncot2,25.07383,-80.32450,7.0
          -80.33433, ... % NCORE: ncot3,25.07950,-80.33433,4.0
          -80.34717, ... % NCORE: ncot4,25.08783,-80.34717,5.6
          ...
          -80.5627, ...         % TAVRK
          -80.5021, ...         % CONSH
          -80.4569, ...         % CONDP
          -80.1885, ...         % BNPIN
          -80.1764, ...         % BNPMI
          -80.1366, ...         % BNPPA
          -80.1664, ...         % BNPON
          -80.1284, ...         % BNPNN
          ...
          -80.21035, ... 	% CCAF1: Coral Restoration Fnd'n / Omics Keys Disease: Carysfort
          -80.29387, ...  	% CNDF1: Coral Restoration Fnd'n / Omics Keys Disease: North Dry Rocks
          -80.30391, ...  	% CGRF1: Coral Restoration Fnd'n / Omics Keys Disease: Grecian Rocks
          -80.41658, ...  	% CPIF1: Coral Restoration Fnd'n / Omics Keys Disease: Pickles Reef
          ...
          -80.0420833, ...	% drbf1,26.4619167,-80.0420833 # Delray Beach (South Central) outfall
          -80.0540500, ...	% bocf1,26.3502667,-80.0540500 # Boca Raton outfall
          -80.0620667, ...	% brwf1,26.2513833,-80.0620667 # Broward outfall
          -80.0859333, ...	% holf1,26.0191167,-80.0859333 # Hollywood outfall
          -80.0862667, ...	% minf1,25.9200500,-80.0862667 # Miami North outfall
          -80.0859667, ...	% micf1,25.7428167,-80.0859667 # Miami Central outfall
          -80.10863, ...	% ho2f1,26.01915,-80.10863 # Current sensor mount near Hollywood outfall
          ...
          -80.6175083, ...      % AOAT_CHEECA_ROCKS_3: Midway between CHCA3 transect start and end points
          -80.2016416, ...      % AOAT_BROAD_KEY_2: Midway between BRDCRK2 transect start and end points
          ...
          -80.6181967099567, ... % AOAT_CHEECA_MAPCO2: Atlantic Ocean Acidification testbed ACTUAL site
          ...
          -80.16071667, ...	% FWC_MK3: Karen Neely @FWC: Marker 3: 25.3734 -80.16071667
          -80.18806667, ...	% FWC_BB: Karen Neely @FWC: Ball Buoy: 25.3162 -80.18806667
          ...
          -80.125417, ...       %SECREMP_MC2
          -80.028567, ...       %SECREMP_PB1
          -80.095967, ...       %SECREMP_BC1
          -80.104033, ...       %SECREMP_DC1
          -80.11173, ...	%SECREMP_UPDB: Upside Down Barge
          -80.111583, ...	%SECREMP_SLIB St. Lucie Inlet Barge 27.200233 -80.111583 19.80
          -80.05632, ...	%SECREMP_EVCA Evans-Crary 27.15572 -80.05632
          -80.12803, ...	%SECREMP_PELA Peck's Lake 27.12212 -80.12803
          -80.01832, ...	%SECREMP_PB2: Outer Reef
          -80.07735, ...	%SECREMP_BC3: Outer Reef
          -80.08810, ...	%SECREMP_DC3: Outer Reef
          ...
          -80.03192, ...	% Luke McEachron @FWC/FWRI: 18HO - 18th Hole
          -80.07912, ...	% Luke McEachron @FWC/FWRI: BLRO - Blowing Rocks
          -80.12242, ...	% Luke McEachron @FWC/FWRI: BSBA - Bullshark Barge
          -80.10998, ...	% Luke McEachron @FWC/FWRI: CEBA - Cement Barge
          -80.10233, ...	% Luke McEachron @FWC/FWRI: CPRE - CLIFTON Perry Reef
          -80.03388, ...	% Luke McEachron @FWC/FWRI: DBNO - DB North
          -80.17043, ...	% Luke McEachron @FWC/FWRI: FPCC,CIVIC Center Reef
          -80.17218, ...	% Luke McEachron @FWC/FWRI: FPLR,FPL Reef
          -80.16336, ...	% Luke McEachron @FWC/FWRI: HORF,House OF Refuge (Valentine)
          -80.09571, ...	% Luke McEachron @FWC/FWRI: LYBR,Lyon Bridge
          -80.04327, ...	% Luke McEachron @FWC/FWRI: MIBA,Middle Bar
          -80.01815, ...	% Luke McEachron @FWC/FWRI: NORO,Nolens Rock
          -80.11605, ...	% Luke McEachron @FWC/FWRI: PIBA,SL Pipe Barge
          -80.02078, ...	% Luke McEachron @FWC/FWRI: SNSO,SS-1
          -80.108367, ...	% Luke McEachron @FWC/FWRI: TERE,Texas Reef
          -80.06662, ...	% Luke McEachron @FWC/FWRI: THSO,THREE Holes (South)
          -80.03178, ...	% Luke McEachron @FWC/FWRI: TUCA,Tunnels/Caves
          -80.058, ...		% Luke McEachron @FWC/FWRI: Noula Express artificial reef, his DEPTH 36.6
          -80.067, ...		% Luke McEachron @FWC/FWRI: Alpha artificial reef, his DEPTH 30.5
          -80.064, ...		% Luke McEachron @FWC/FWRI: Rodeo artificial reef, his DEPTH 30.5
          ...
          -80.0878, ...		% Alex Soloviev @NSUOC: SFOCB: 'c-buoy',[26.0728,-80.0878, 20, 2]
          -80.0768, ...		% Alex Soloviev @NSUOC: SFOEB: 'e-buoy',[26.0695,-80.0768, 50, 3]
          -80.0768, ...		% Alex Soloviev @NSUOC: SFONE: 'ne-buoy',[26.0695,-80.0768, 50, 3]
          -80.0942, ...		% Alex Soloviev @NSUOC: SFONW: 'nw_w-btm',[26.0705,-80.0942, 11, 10]
          -80.0533, ...		% Alex Soloviev @NSUOC: SFOSE: 'se-buoy',[26.1398,-80.0533, 145, 8]
          -80.0918, ...		% Alex Soloviev @NSUOC: SFOSW: 'sw-buoy',[26.0327,-80.0918, 20, 1]
          -80.1090, ...		% Alex Soloviev @NSUOC: SFOPI: 'pier-cc',[26.0577,-80.1090, 0, 4]
          ...
          -80.225, ...		%41114 Fort Pierce, FL (134)
          -80.184, ...		%41009 (LLNR 840) - CANAVERAL 20 NM East of Cape Canaveral "CANF1"
          -80.530, ...		%41113 - Cape Canaveral Nearshore, FL (143) "CNNF1"
          -80.534, ...		%Station 41012 (LLNR 845.3) - 40NM ENE of St Augustine FL
          ...
          -94.3121, ...		% SBFG1
          -93.6255, ...		% EBFG1
          -93.8397, ...		% WBFG1
          -92.4829, ...		% SOFG1
          ...
          -78.99472, ...	%spgf1,26.70444,-78.99472 # Station SPGF1 - Settlement Point, GBI, Bahamas
          ...
          -169.658181, ...	%ofsp6,-14.184164,-169.658181: American Samoa - Ofu - Craig et al study sites
          -169.627813, ...	%onsp6,-14.164692,-169.627813: American Samoa - Ofu - north coast
          -170.668, ...		%pgsp6,-14.296,-170.668: American Samoa - Tutuila - just outside Pago Pago harbor
          -170.661, ...		%vesp6,-14.235,-170.661: American Samoa - Tutuila - just east of Vatia MPA
          -170.710, ...		%vwsp6,-14.252,-170.710: American Samoa - Tutuila - just west of Vatia MPA
          -170.590, ...		%aosp6,-14.251,-170.590: American Samoa - Tutuila - NNW of Aoa village
          ...
          -64.975, ...		%pbbe1,32.29,-64.975 # Pilchard Bay - Bermuda (see also BERM1) - 5 m
          -64.863, ...		%swbe1,32.23,-64.863 # Southwest Breaker Area - Bermuda (see also BERM1) - 6 m
          -64.77, ...		%nrbe1,32.475,-64.77 # North Rock - Bermuda (see also BERM1) - 2 m
      ];
    ALL_STATIONS.lats = ...
      [ ...
          23.791, ...           % CMRC3
          17.784, ...           % SRVI2
          17.939, ...           % LPPR1
          18.470, ...           % DBJM1
          19.699, ...           % LCIY2
          17.769, ...           % 41140 - NDBC buoy, Christiansted, St. Croix
          15.156780, ...        % LLBP7
          26.092, ...           % PVGF1 - Port Everglades channel, Florida
          18.3449, ...		% Brewers Bay - St. Thomas - USVI - PLANNED - BBVI1,18.3449,-64.9843
          25.74897, ...		% Abs. turbidity cal/val site neardof Port of Miami channel - pomf1,25.74897,-80.13317
          15.24, ...		% manm1,15.24,145.715: Managaha - Saipan - CNMI
          14.122, ...		% sasp7,14.122,145.159: Sasanhaya Bay - Rota - CNMI
          14.106, ...		% tasp7,14.106,145.192: Talakhaya/Sabana Watershed - Rota - CNMI
          14.9573, ...		% tinp7,14.9573,145.6194: Tinian Harbor (reefs) - Tinian - CNMI
          13.428056, ...	% pgbp7,13.428056,144.796389: NOAA NOS Station PGBP7 - 1631428 - Pago Bay, Guam
          21.303333, ...        % oouh1 - 1612340 - Honolulu, HI - 21.303333,-157.864444,-1.5
          24.416667, ...	% nwbh1,24.416667,-162.100000 - NDBC Station 51001 (LLNR 28006) - NW HI 1 - 170 NM WNW of Kauai
          ...
          19.833, ...		% ppdr1,19.833,-70.731: CREWS Dominican Republic - Puerto Plata,data from 1998-01-01 to 2016-08-25
          18.432, ...		% cwdr1,18.432,-69.580: CREWS - Dominican Republic - Catuan Wreck, data from 1998-01-01 to 2016-08-25
          13.184, ...		% drbb1,13.184,-59.646: CREWS - Barbados - Dottin's Reef, data from 1998-01-01 to 2016-08-25
          17.194325, ...	% hcbz1,17.194325,-87.521850: Half-Moon Caye - Belize
          17.273567, ...	% ccbz1,17.273567,-87.810000: Calabash Caye - Belize
          16.815650, ...	% swbz1,16.815650,-88.076867: South Water Caye - Belize
          11.1760, ...		% buto1,11.1760,-60.8338: Buccoo Marine Park - Tobago
          11.3009, ...		% arto1,11.3009,-60.5206: Angel Reef - Tobago
          17.0155544, ...	% mhab1,17.0155544,-61.8722687 # Monk's Head - Antigua and Barbuda
          17.3563538,...	% prsk1,17.3563538,-62.8548927 # Paradise Reef - St. Kitts
          13.859093,...		% smsl1,13.859093,-61.071307 # Devil's Hole - Soufriere Marine Management Association (SMMA) - St. Lucia
          12.022750,...		% krgn1,12.022750,-61.790983 # Kahonae Reef - Grenada
          8.7251, ...		% pisr1,8.7251,81.2080: Pigeon Island National Park - Sri Lanka - Indian Ocean
          ...
          26.612, ...           % LKWF1
          25.731667, ...        % VAKF1
          25.593001, ...        % FWYF1
          25.222, ...           % CRYF1
          25.011779, ...        % MLRF1
          24.95306, ...         % AQUA1
          24.843, ...           % LONF1
          25.084444, ...        % NFBF1
          24.746053, ...        % TNRF1
          24.916666, ...        % PKYF1
          24.711666, ...        % VCAF1
          24.553333, ...        % KYWF1
          24.70033, ...         % MOSE1
          24.627851, ...        % SMKF1
          24.54275, ...         % LOOE1
          24.525, ...           % AMSF1
          24.460, ...           % SANF1
          24.693333, ...        % PLSF1
          24.638333, ...        % DRYF1
          24.627, ...           % gkyf1,24.627,-82.872
          26.063611, ...        % 42023 - C13 - West Florida South Buoy
          25.966, ...           % 42003 - NDBC buoy near edge of Loop Current
          27.862, ...           % Sebastian Inlet State Park, FL: 27.862 N 80.445 W (27°51'42" N 80°26'41" W
          ...
          25.0313833, ... % NCORE: klgf1,25.0313833,-80.3480167,23 # Key Largo Current/Temperature 7/23 m (offshore French Reef)
          24.7403000, ... % NCORE: mrtf1,24.7403000,-80.7763667,23 # Marathon Current/Temperature 8/23 m (offshore Tennessee Reef)
          ...
          25.1090, ... % NCORE: ncora,25.1090,-80.3803,4.1
          25.0932, ... % NCORE: ncorb,25.0932,-80.3550,7.0
          25.0673, ... % NCORE: ncorc,25.0673,-80.3183,21.8
          25.0740, ... % NCORE: ncor1,25.0740,-80.3178,26.0
          25.0733, ... % NCORE: ncor2,25.0733,-80.3183,21.8
          25.06900, ... % NCORE: ncot1,25.06900,-80.31950,9.75
          25.07383, ... % NCORE: ncot2,25.07383,-80.32450,7.0
          25.07950, ... % NCORE: ncot3,25.07950,-80.33433,4.0
          25.08783, ... % NCORE: ncot4,25.08783,-80.34717,5.6
          ...
          24.9390, ...          % TAVRK
          24.9465, ...          % CONSH
          24.9465, ...          % CONDP
          25.3613, ...          % BNPIN
          25.3641, ...          % BNPMI
          25.3732, ...          % BNPPA
          25.3626, ...          % BNPON
          25.4734, ...          % BNPNN
          ...
          25.22073, ...  	% CCAF1: Coral Restoration Fnd'n / Omics Keys Disease: Carysfort
          25.13035, ...  	% CNDF1: Coral Restoration Fnd'n / Omics Keys Disease: North Dry Rocks
          25.11014, ...  	% CGRF1: Coral Restoration Fnd'n / Omics Keys Disease: Grecian Rocks
          24.98462, ...  	% CPIF1: Coral Restoration Fnd'n / Omics Keys Disease: Pickles Reef
          ...
          26.4619167, ...	% drbf1,26.4619167,-80.0420833 # Delray Beach (South Central) outfall
          26.3502667, ...	% bocf1,26.3502667,-80.0540500 # Boca Raton outfall
          26.2513833, ...	% brwf1,26.2513833,-80.0620667 # Broward outfall
          26.0191167, ...	% holf1,26.0191167,-80.0859333 # Hollywood outfall
          25.9200500, ...	% minf1,25.9200500,-80.0862667 # Miami North outfall
          25.7428167, ...	% micf1,25.7428167,-80.0859667 # Miami Central outfall
          26.01915, ...		% ho2f1,26.01915,-80.10863 # Current sensor mount near Hollywood outfall
          ...
          24.9008417, ...       % AOAT_CHEECA_ROCKS_3: Midway between CHCA3 transect start and end points
          25.309875, ...        % AOAT_BROAD_KEY_2: Midway between BRDCRK2 transect start and end points
          ...
          24.8978019047619, ... % AOAT_CHEECA_MAPCO2: Atlantic Ocean Acidification testbed ACTUAL site
          ...
          25.3734, ...		% FWC_MK3: Karen Neely @FWC: Marker 3: 25.3734 -80.16071667
          25.3162, ...		% FWC_BB: Karen Neely @FWC: Ball Buoy: 25.3162 -80.18806667
          ...
          27.112033, ...        %SECREMP_MC2
          26.709717, ...        %SECREMP_PB1
          26.147867, ...        %SECREMP_BC1
          25.842167, ...        %SECREMP_DC1
          27.232433, ...	%SECREMP_UPDB: Upside Down Barge
          27.200233, ...	%SECREMP_SLIB St. Lucie Inlet Barge 27.200233 -80.111583 19.80
          27.15572, ...		%SECREMP_EVCA Evans-Crary 27.15572 -80.05632
          27.12212, ...		%SECREMP_PELA Peck's Lake 27.12212 -80.12803
          26.678750, ...	%SECREMP_PB2: Outer Reef
          26.158633, ...	%SECREMP_BC3: Outer Reef
          25.842100, ...	%SECREMP_DC3: Outer Reef
          ...
          27.08772, ...		% Luke McEachron @FWC/FWRI: 18HO - 18th Hole
          26.99117, ...		% Luke McEachron @FWC/FWRI: BLRO - Blowing Rocks
          27.13702, ...		% Luke McEachron @FWC/FWRI: BSBA - Bullshark Barge
          27.21065, ...		% Luke McEachron @FWC/FWRI: CEBA - Cement Barge
          27.2253, ...		% Luke McEachron @FWC/FWRI: CPRE - Clifton Perry Reef
          26.91253, ...		% Luke McEachron @FWC/FWRI: DBNO - DB North
          27.44612, ...		% Luke McEachron @FWC/FWRI: FPCC,Civic Center Reef
          27.44485, ...		% Luke McEachron @FWC/FWRI: FPLR,FPL Reef
          27.19923, ...		% Luke McEachron @FWC/FWRI: HORF,House of Refuge (Valentine)
          27.21741, ...		% Luke McEachron @FWC/FWRI: LYBR,Lyon Bridge
          27.04505, ...		% Luke McEachron @FWC/FWRI: MIBA,Middle Bar
          27.00478, ...		% Luke McEachron @FWC/FWRI: NORO,Nolens Rock
          27.22275, ...		% Luke McEachron @FWC/FWRI: PIBA,SL Pipe Barge
          26.86307, ...		% Luke McEachron @FWC/FWRI: SNSO,SS-1
          27.191967, ...	% Luke McEachron @FWC/FWRI: TERE,Texas Reef
          27.00457, ...		% Luke McEachron @FWC/FWRI: THSO,Three Holes (South)
          26.94715, ...		% Luke McEachron @FWC/FWRI: TUCA,Tunnels/Caves
          26.322, ...		% Luke McEachron @FWC/FWRI: Noula Express artificial reef, his depth 36.6
          26.231, ...		% Luke McEachron @FWC/FWRI: Alpha artificial reef, his depth 30.5
          26.231, ...		% Luke McEachron @FWC/FWRI: Rodeo artificial reef, his depth 30.5
          ...
          26.0728, ...		% Alex Soloviev @NSUOC: SFOCB: 'c-buoy',[26.0728,-80.0878, 20, 2]
          26.0695, ...		% Alex Soloviev @NSUOC: SFOEB: 'e-buoy',[26.0695,-80.0768, 50, 3]
          26.0695, ...		% Alex Soloviev @NSUOC: SFONE: 'ne-buoy',[26.0695,-80.0768, 50, 3]
          26.0705, ...		% Alex Soloviev @NSUOC: SFONW: 'nw_w-btm',[26.0705,-80.0942, 11, 10]
          26.1398, ...		% Alex Soloviev @NSUOC: SFOSE: 'se-buoy',[26.1398,-80.0533, 145, 8]
          26.0327, ...		% Alex Soloviev @NSUOC: SFOSW: 'sw-buoy',[26.0327,-80.0918, 20, 1]
          26.0577, ...		% Alex Soloviev @NSUOC: SFOPI: 'pier-cc',[26.0577,-80.1090, 0, 4]
          ...
          27.551, ...		%41114 Fort Pierce, FL (134)
          28.523, ...		%41009 (LLNR 840) - CANAVERAL 20 NM East of Cape Canaveral "CANF1"
          28.400, ...		%41113 - Cape Canaveral Nearshore, FL (143)  "CNNF1"
          30.042, ...		%Station 41012 (LLNR 845.3) - 40NM ENE of St Augustine FL
          ...
          28.1785, ...		% SBFG1
          27.9265, ...		% EBFG1
          27.8828, ...		% WBFG1
          28.3527, ...		% SOFG1
          ...
          26.70444, ...		%spgf1,26.70444,-78.99472 # Station SPGF1 - Settlement Point, GBI, Bahamas
          ...
          -14.184164, ...	%ofsp6,-14.184164,-169.658181: American Samoa - Ofu - Craig et al study sites
          -14.164692, ...	%onsp6,-14.164692,-169.627813: American Samoa - Ofu - north coast
          -14.296, ...		%pgsp6,-14.296,-170.668: American Samoa - Tutuila - just outside Pago Pago harbor
          -14.235, ...		%vesp6,-14.235,-170.661: American Samoa - Tutuila - just east of Vatia MPA
          -14.252, ...		%vwsp6,-14.252,-170.710: American Samoa - Tutuila - just west of Vatia MPA
          -14.251, ...		%aosp6,-14.251,-170.590: American Samoa - Tutuila - NNW of Aoa village
          ...
          32.29, ...		%pbbe1,32.29,-64.975 # Pilchard Bay - Bermuda (see also BERM1) - 5 m
          32.23, ...		%swbe1,32.23,-64.863 # Southwest Breaker Area - Bermuda (see also BERM1) - 6 m
          32.475, ...		%nrbe1,32.475,-64.77 # North Rock - Bermuda (see also BERM1) - 2 m
      ];
    ALL_STATIONS.depths = ...
      [ ...
          6.0, ...              % CMRC3
          6.0, ...              % SRVI2
          6.0, ...              % LPPR1
          6.0, ...              % DBJM1
          5.5, ...              % LCIY2
          244.0, ...            % 41140 - NDBC buoy, Christiansted, St. Croix
          6.0, ...              % LLBP7
          1.0, ...              % PVGF1 - Port Everglades channel, Florida
          11.0, ...		% Brewers Bay - St. Thomas - USVI - PLANNED - BBVI1,18.3449,-64.9843
          4.3, ...		% Abs. turbidity cal/val site neardof Port of Miami channel - pomf1,25.74897,-80.13317
          3.0, ...		% manm1,15.24,145.715: Managaha - Saipan - CNMI (NGDC 6-arcsecond bathymetry)
          100, ...		% sasp7,14.122,145.159: Sasanhaya Bay - Rota - CNMI (bathymetry from Google Earth)
          280, ...		% tasp7,14.106,145.192: Talakhaya/Sabana Watershed - Rota - CNMI
          3, ...		% tinp7,14.9573,145.6194: Tinian Harbor (reefs) - Tinian - CNMI
          0, ...		% pgbp7,13.428056,144.796389: NOAA NOS Station PGBP7 - 1631428 - Pago Bay, Guam
          1.5, ...              % oouh1 - 1612340 - Honolulu, HI - 21.303333,-157.864444,-1.5
          4865, ...		% SEA TEMP 1 m below MLLW - nwbh1,24.416667,-162.100000 - NDBC Station 51001 (LLNR 28006) - NW HI 1 - 170 NM WNW of Kauai
          ...
          6, ...		% ppdr1,19.833,-70.731: CREWS Dominican Republic - Puerto Plata,data from 1998-01-01 to 2016-08-25
          6, ...		% cwdr1,18.432,-69.580: CREWS - Dominican Republic - Catuan Wreck, data from 1998-01-01 to 2016-08-25
          6, ...		% drbb1,13.184,-59.646: CREWS - Barbados - Dottin's Reef, data from 1998-01-01 to 2016-08-25
          6, ...		% hcbz1,17.194325,-87.521850: Half-Moon Caye - Belize
          6, ...		% ccbz1,17.273567,-87.810000: Calabash Caye - Belize
          6, ...		% swbz1,16.815650,-88.076867: South Water Caye - Belize
          6, ...		% buto1,11.1760,-60.8338: Buccoo Marine Park - Tobago
          6, ...		% arto1,11.3009,-60.5206: Angel Reef - Tobago
          6, ...		% mhab1,17.0155544,-61.8722687 # Monk's Head - Antigua and Barbuda
          6,...			% prsk1,17.3563538,-62.8548927 # Paradise Reef - St. Kitts
          6,...			% smsl1,13.859093,-61.071307 # Devil's Hole - Soufriere Marine Management Association (SMMA) - St. Lucia
          6,...			% krgn1,12.022750,-61.790983 # Kahonae Reef - Grenada
          6, ...		% pisr1,8.7251,81.2080: Pigeon Island National Park - Sri Lanka - Indian Ocean
          ...
          1.0, ...              % LKWF1
          1.0, ...              % VAKF1
          3.0, ...              % FWYF1
          11.0, ...             % CRYF1
          3.55, ...             % MLRF1 - was 2.0! Reestimated from nearest NGDC CRM points
          22.8, ...             % AQUA1 is "at 60 feet", but SFP site mean depth ~22.8m
          1.28, ...             % LONF1 - was 1.0! Reestimated from 1-D Wave i-depth...
          1.00, ...             % NFBF1
          6.0, ...              % TNRF1
          1.0, ...              % PKYF1
          1.0, ...              % VCAF1
          1.0, ...              % KYWF1
          1.2, ...              % MOSE1
          3.18, ...             % SMKF1 - was 2.0! Reestimated from nearest NGDC CRM points
          22, ...               % LOOE1 - per SFROS / SFP Database Web page
          2.0, ...              % AMSF1
          1.0, ...              % SANF1
          2.0, ...              % PLSF1
          1.0, ...              % DRYF1
          1.5, ...              % gkyf1,24.627,-82.872
          49.0, ...             % 42023 - C13 - West Florida South Buoy - est. from NGDC CRM
          3282.7, ...           % 42003 - NDBC buoy near edge of Loop Current
          1, ...                % Sebastian Inlet State Park, FL: 27.862 N 80.445 W (27°51'42" N 80°26'41" W
          ...
          23, ... % NCORE: klgf1,25.0313833,-80.3480167,23 # Key Largo Current/Temperature 7/23 m (offshore French Reef)
          23, ... % NCORE: mrtf1,24.7403000,-80.7763667,23 # Marathon Current/Temperature 8/23 m (offshore Tennessee Reef)
          ...
          4.1, ... % NCORE: ncora,25.1090,-80.3803,4.1
          7.0, ... % NCORE: ncorb,25.0932,-80.3550,7.0
          21.8, ... % NCORE: ncorc,25.0673,-80.3183,21.8
          26.0, ... % NCORE: ncor1,25.0740,-80.3178,26.0
          21.8, ... % NCORE: ncor2,25.0733,-80.3183,21.8
          9.75, ... % NCORE: ncot1,25.06900,-80.31950,9.75
          7.0, ... % NCORE: ncot2,25.07383,-80.32450,7.0
          4.0, ... % NCORE: ncot3,25.07950,-80.33433,4.0
          5.6, ... % NCORE: ncot4,25.08783,-80.34717,5.6
          ...
          3.9624 , ...          % TAVRK
          5.4864 , ...          % CONSH
          16.4592, ...          % CONDP
          3.0480 , ...          % BNPIN
          4.5720 , ...          % BNPMI
          12.1920, ...          % BNPPA
          4.5720 , ...          % BNPON
          5.4864 , ...          % BNPNN
          ...
          8.5, ...		% CCAF1: Coral Restoration Fnd'n / Omics Keys Disease: Carysfort
          7.0, ...		% CNDF1: Coral Restoration Fnd'n / Omics Keys Disease: North Dry Rocks
          7.7, ...		% CGRF1: Coral Restoration Fnd'n / Omics Keys Disease: Grecian Rocks
          5.4, ... 		% CPIF1: Coral Restoration Fnd'n / Omics Keys Disease: Pickles Reef
          ...
          27, ...		% drbf1,26.4619167,-80.0420833 # Delray Beach (South Central) outfall
          27, ...		% bocf1,26.3502667,-80.0540500 # Boca Raton outfall
          27, ...		% brwf1,26.2513833,-80.0620667 # Broward outfall
          27, ...		% holf1,26.0191167,-80.0859333 # Hollywood outfall
          27, ...		% minf1,25.9200500,-80.0862667 # Miami North outfall
          27, ...		% micf1,25.7428167,-80.0859667 # Miami Central outfall
          7, ...		% ho2f1,26.01915,-80.10863 # ADCP sensor mounted near Hollywood outfall
          ...
          4.6, ...              % AOAT_CHEECA_ROCKS_3
          3.9, ...              % AOAT_BROAD_KEY_2
          ...
          4.3, ...              % AOAT_CHEECA_MAPCO2: Atlantic Ocean Acidification testbed ACTUAL site (source: USGS 30 m bathymetry)
          ...
          4.5, ...		% FWC_MK3: Karen Neely @FWC: Marker 3: 25.3734 -80.16071667: ESTIMATE FROM NGDC
          5.3, ...		% FWC_BB: Karen Neely @FWC: Ball Buoy: 25.3162 -80.18806667: ESTIMATE FROM NGDC
          ...
          4.5, ...              %SECREMP_MC2
          7.6, ...              %SECREMP_PB1
          7.6, ...              %SECREMP_BC1
          7.6, ...              %SECREMP_DC1
          19.2, ...             %SECREMP_UPDB: Upside Down Barge
          19.8, ...		%SECREMP_SLIB St. Lucie Inlet Barge 27.200233 -80.111583 19.80
          18.6, ...		%DEPTH ESTIMATED FROM NGDC! SECREMP_EVCA Evans-Crary 27.15572 -80.05632
          10.0, ...		%DEPTH ESTIMATED FROM NGDC! SECREMP_PELA Peck's Lake 27.12212 -80.12803
          16.7, ...             %SECREMP_PB2: Outer Reef
          16.7, ...             %SECREMP_BC3: Outer Reef
          16.7, ...             %SECREMP_DC3: Outer Reef
          ...
          22.9, ...		% Luke McEachron @FWC/FWRI: 18HO - 18th Hole                    
          8.3, ...		% Luke McEachron @FWC/FWRI: BLRO - Blowing Rocks                
          11.7, ...		% Luke McEachron @FWC/FWRI: BSBA - Bullshark Barge              
          18.3, ...		% Luke McEachron @FWC/FWRI: CEBA - Cement Barge                 
          19.2, ...		% Luke McEachron @FWC/FWRI: CPRE - Clifton Perry Reef           
          22.6, ...		% Luke McEachron @FWC/FWRI: DBNO - DB North                     
          15.9, ...		% Luke McEachron @FWC/FWRI: FPCC,Civic Center Reef              
          15.7, ...		% Luke McEachron @FWC/FWRI: FPLR,FPL Reef                       
          5.0, ...		% Luke McEachron @FWC/FWRI: HORF,House of Refuge (Valentine)    
          18.5, ...		% Luke McEachron @FWC/FWRI: LYBR,Lyon Bridge                    
          19.9, ...		% Luke McEachron @FWC/FWRI: MIBA,Middle Bar                     
          21.3, ...		% Luke McEachron @FWC/FWRI: NORO,Nolens Rock                    
          17.9, ...		% Luke McEachron @FWC/FWRI: PIBA,SL Pipe Barge                  
          26.2, ...		% Luke McEachron @FWC/FWRI: SNSO,SS-1                           
          16.9, ...		% Luke McEachron @FWC/FWRI: TERE,Texas Reef                     
          14.3, ...		% Luke McEachron @FWC/FWRI: THSO,Three Holes (South)            
          21.5, ...		% Luke McEachron @FWC/FWRI: TUCA,Tunnels/Caves                  
          28.8, ...		% Luke McEachron @FWC/FWRI: Noula Express artificial reef, his depth 36.6               
          25.2, ...		% Luke McEachron @FWC/FWRI: Alpha artificial reef, his depth 30.5               
          36.7, ...		% Luke McEachron @FWC/FWRI: Rodeo artificial reef, his depth 30.5		
          ...
          20, ...		% Alex Soloviev @NSUOC: SFOCB: 'c-buoy',[26.0728,-80.0878, 20, 2]
          50, ...		% Alex Soloviev @NSUOC: SFOEB: 'e-buoy',[26.0695,-80.0768, 50, 3]
          50, ...		% Alex Soloviev @NSUOC: SFONE: 'ne-buoy',[26.0695,-80.0768, 50, 3]
          11, ...		% Alex Soloviev @NSUOC: SFONW: 'nw_w-btm',[26.0705,-80.0942, 11, 10]
          145, ...		% Alex Soloviev @NSUOC: SFOSE: 'se-buoy',[26.1398,-80.0533, 145, 8]
          20, ...		% Alex Soloviev @NSUOC: SFOSW: 'sw-buoy',[26.0327,-80.0918, 20, 1]
          0, ...		% Alex Soloviev @NSUOC: SFOPI: 'pier-cc',[26.0577,-80.1090, 0, 4]
          ...
          16.15, ...		%41114 Fort Pierce, FL (134)
          40.5, ...		%41009 (LLNR 840) - CANAVERAL 20 NM East of Cape Canaveral "CANF1"
          9.9, ...		%41113 - Cape Canaveral Nearshore, FL (143) "CNNF1"
          38.1, ...		%Station 41012 (LLNR 845.3) - 40NM ENE of St Augustine FL
          ...
          30, ...		% SBFG1
          30, ...		% EBFG1
          30, ...		% WBFG1
          30, ...		% SOFG1
          ...
          0, ...		%spgf1,26.70444,-78.99472 # Station SPGF1 - Settlement Point, GBI, Bahamas
          ...
          8, ...		%ofsp6,-14.184164,-169.658181: American Samoa - Ofu - Craig et al study sites (NGDC 3 arcsec bathymetry)
          17, ...		%onsp6,-14.164692,-169.627813: American Samoa - Ofu - north coast
          71, ...		%pgsp6,-14.296,-170.668: American Samoa - Tutuila - just outside Pago Pago harbor (NGDC 3 arcsec bathymetry)
          55, ...		%vesp6,-14.235,-170.661: American Samoa - Tutuila - just east of Vatia MPA (NGDC 3 arcsec bathymetry)
          55, ...		%vwsp6,-14.252,-170.710: American Samoa - Tutuila - just west of Vatia MPA (NGDC 3 arcsec bathymetry)
          28, ...		%aosp6,-14.251,-170.590: American Samoa - Tutuila - NNW of Aoa village
          ...
          5, ...		%pbbe1,32.29,-64.975 # Pilchard Bay - Bermuda (see also BERM1) - 5 m
          6, ...		%swbe1,32.23,-64.863 # Southwest Breaker Area - Bermuda (see also BERM1) - 6 m
          2, ...		%nrbe1,32.475,-64.77 # North Rock - Bermuda (see also BERM1) - 2 m
      ];

  end;

  % SANITY CHECK - this file is frequently edited by hand
  if ( numel(ALL_STATIONS.codes) ~= numel(ALL_STATIONS.lats) || ...
       numel(ALL_STATIONS.codes) ~= numel(ALL_STATIONS.lons) || ...
       numel(ALL_STATIONS.codes) ~= numel(ALL_STATIONS.depths) )
    error('METADATA LENGTH MISMATCH!');
  end;

  STATIONS = ALL_STATIONS;

return;
