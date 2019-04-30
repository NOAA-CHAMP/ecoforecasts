1;

if ( ~exist('srvi2','var') || isempty(srvi2) )
  srvi2 = load_station_data('srvi2');
end;

disp('Calc Kd');
srvi2 = calc_kd_daily(srvi2, 'bic_surf_par', 'bic_shallow_par', 2);
disp('Regress Kd');
regress_kd;
