fld='ndbc_erai_30a_net_heat_flux';

begix = find(isfinite(stn.(fld).data),1);
endix = find(isfinite(stn.(fld).data),1,'last');

dts = stn.(fld).date(begix):(1/24):stn.(fld).date(endix);
dat = interp1(stn.(fld).date,real(stn.(fld).data),dts,'pchip');


fmg; wt([dts ; dat], 'BlackandWhite');
