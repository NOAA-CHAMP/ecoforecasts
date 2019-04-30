1;


cstnms = { ...
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
         };

dx = 1e3; dy = 1e3;
for cs = cstnms(:)';
  stnm = cs{:};
  stn = read_hires_bathymetry(stnm,[dy*1.05,dx*1.10],[],true);
  stn.depth = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.field,stn.lon,stn.lat);
  fprintf(1,'%s, %3.1f\n',stnm,stn.depth);
  stn=[]; clear stn
end;

clear ans cs cstnms dx dy stnm
