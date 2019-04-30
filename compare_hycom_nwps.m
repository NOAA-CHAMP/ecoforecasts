
w = warning('OFF','Ecoforecasts:NWPS:NoFile');
flds = get_nwps_fields();
warning(w); clear w

a = stokes_drift(flds.windspeed.field(:),flds.sigwavehgt.field(:),flds.primwaveper.field(:),'ardhuin');
a = reshape(a,size(flds.windspeed.field));

m = stokes_drift(flds.windspeed.field(:),flds.sigwavehgt.field(:),flds.primwaveper.field(:),'monismith',20);
m = reshape(m,size(flds.windspeed.field));

if 1;
  fmg;
  contour_field(flds.monismith_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Monismith Lagrangian velocity (Python)');

  fmg;
  contour_field(flds.monismith_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Monismith Lagrangian velocity (MATLAB)');
end;

if 1;
  fmg;
  contour_field(flds.monismith_surface_drift,@contourf,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('93^{rd} percentile Monismith Lagrangian velocity (Python)');

  fmg;
  contour_field(flds.monismith_surface_drift,@contourf,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('93^{rd} percentile Monismith Lagrangian velocity (MATLAB)');
end;

if 1;
  fmg;
  contour_field(flds.ardhuin_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Ardhuin Lagrangian velocity (Python)');

  fmg;
  contour_field(flds.ardhuin_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Ardhuin Lagrangian velocity (MATLAB)');
end;

if 1;
  fmg;
  contour_field(flds.ardhuin_surface_drift,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.3]);
  titlename('93^{rd} percentile Ardhuin Lagrangian velocity (MATLAB)');

  fmg;
  contour_field(flds.ardhuin_surface_drift,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.3]);
  titlename('93^{rd} percentile Ardhuin Lagrangian velocity (MATLAB)');
end;
