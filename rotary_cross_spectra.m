1;

smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1);
smkf1 = station_spddir_to_uv(smkf1,'ndbc_wind1_speed','ndbc_wind1_dir');
looe1 = get_looe1_adcp;

u1 = smkf1.ndbc_wind1_u; v1 = smkf1.ndbc_wind1_v;
u2 = looe1.adcp_u; v2 = looe1.adcp_v;

[u1,v1,u2,v2] = intersect_tses(u1,v1,u2,v2);

u1.data = u1.data - nanmean(u1.data);
v1.data = v1.data - nanmean(v1.data);
u2.data = u2.data - nanmean(u2.data);
v2.data = v2.data - nanmean(v2.data);

% JLAB/MSPEC:
%   mspec  Multitaper power and cross spectra.
%   
%     mspec implements spectral and cross-spectral analysis using the 
%     multitaper method for real or complex-valued data. 
%  
%     mspec is to be run after calling SLEPTAP to compute the multitapers.
%     ______________________________________________________________________
% 
%     ...
%     Cross-spectra of complex-valued data
%  
%     To compute the cross-spectra of two complex-valued time series or sets 
%     of time series Z1 and Z2, run mspec repeatedly.
% 
%     [F,SP1P1,SP2P2,SP1P2]=mspec(Z1,Z2,PSI);  
%     [F,SN1N1,SN2N2,SN1N2]=mspec(CONJ(Z1),CONJ(Z2),PSI);  
%  
%     The first call returns the spectra and cross-spectra of Z1 and Z2 at
%     positive frequencies, while the second returns their spectra and cross-
%     spectra at negative frequencies.  Finally
% 
%     [F,SP1P1,SN2N2,SP1N2]=mspec(Z1,CONJ(Z2),PSI); 

[PSI,LAMBDA] = sleptap(numel(u1.data));

%     The first call returns the spectra and cross-spectra of Z1 and Z2 at
%     positive frequencies, while the second returns their spectra and cross-
%     spectra at negative frequencies.  Finally, [... return] cross-spectra
%     between positive rotary components of Z1 and negative components of Z2.

Z1 = u1.data + (i*v1.data);
Z2 = u2.data + (i*v2.data);

DT = median(diff(u1.date));
[F,SP1P1,SP2P2,SP1P2]=mspec(DT,Z1,Z2,PSI);
[F,SN1N1,SN2N2,SN1N2]=mspec(DT,conj(Z1),conj(Z2),PSI);  
[F,SP1P1,SN2N2,SP1N2]=mspec(DT,Z1,conj(Z2),PSI); 

fmg; spt(1,2,1); semilogx(F,SN1N2); set(gca,'XDir','reverse'); spt(1,2,2); semilogx(F,SP1P2);
