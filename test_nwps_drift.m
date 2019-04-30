1;


flds = get_nwps_fields(); m = stokes_drift(flds.windspeed.field,flds.sigwavehgt.field,flds.primwaveper.field,'monismith',3); fmg; contour_field(flds.monismith_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]); set(gca,'CLim',[0,0.20]); fmg; contour_field(flds.ardhuin_surface_drift,[],{@prctile,93}); set(gca,'CLim',[0,0.3]);
