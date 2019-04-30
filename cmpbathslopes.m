1;

set_more off

% List from GET_ALL_STATION_METADATA.M:

stnms = {
          'POMF1', ...		% Abs. turbidity cal/val site neardof Port of Miami channel - pomf1,25.74897,-80.13317
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
          'TAVRK', ...
          'CONSH', ...
          'CONDP', ...
          'BNPIN', ...
          'BNPMI', ...
          'BNPPA', ...
          'BNPON', ...
          'BNPNN', ...
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
          'SECREMP_NOULA', ...	% Luke McEachron @FWC/FWRI: Noula, his depth 36.6
          'SECREMP_ALPHA', ...	% Luke McEachron @FWC/FWRI: Alpha, his depth 30.5
          'SECREMP_RODEO', ...	% Luke McEachron @FWC/FWRI: Rodeo, his depth 30.5
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
        };


myds = [];
mybs = [];
myas = [];
lods = [];
lobs = [];
loas = [];
hids = [];
hibs = [];
hias = [];

for stix = 1:numel(stnms)
  stnm = stnms{stix};
  disp(stnm);
  stn=[]; clear stn
  stn = get_station_from_station_name(stnm);
  myds(stix) = -stn.depth;
  try,
    stn = station_ngdc_offshore_slope(stn);
    mybs(stix) = stn.ngdc_offshore_slope;
  catch,
    mybs(stix) = nan;
  end;
  try,
    stn = station_optimal_isobath_orientation(stn);
    myas(stix) = mod( (stn.isobath_orientation + 90), 360 );
  catch,
    myas(stix) = nan;
  end;

  stn = read_hires_bathymetry(stn,[1e3,1e3],[],false);
  % 3-point (279-m) stencil
  [lobs(stix),loas(stix)] = find_ngdc_slope(stn.ngdc_hires_bathy,stn.lon,stn.lat,3);
  lods(stix) = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.field,stn.lon,stn.lat);
  stn = rmfield(stn,'ngdc_hires_bathy');

  stn = read_hires_bathymetry(stn,[1e3,1e3],[],true);
  % 9-point (90- or 270-m) stencil
  [hibs(stix),hias(stix)] = find_ngdc_slope(stn.ngdc_hires_bathy,stn.lon,stn.lat,9);
  hids(stix) = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.field,stn.lon,stn.lat);
end;

% If any of our bathymetry is wacky
loas(lods < -300) = nan;
lobs(lods < -300) = nan;
lods(lods < -300) = nan;
hias(hids < -300) = nan;
hibs(hids < -300) = nan;
hids(hids < -300) = nan;

loas(lods >= 0) = nan;
lobs(lods >= 0) = nan;
lods(lods >= 0) = nan;
hias(hids >= 0) = nan;
hibs(hids >= 0) = nan;
hids(hids >= 0) = nan;

set_more
