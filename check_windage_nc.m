1;

fname = 'NATLU_0.0m_2001(01-15).nc';

nc = mDataset(['d:/connectivity/NATLU/',fname]);
u = cast(nc{'u_surface'}(:,:,:),'double');
v = cast(nc{'v_surface'}(:,:,:),'double');
close(nc); clear nc
nc = mDataset(['d:/connectivity/NATLU_0.01_Windage/',fname]);
uw = cast(nc{'u_surface'}(:,:,:),'double');
vw = cast(nc{'v_surface'}(:,:,:),'double');
close(nc); clear nc

fmg; quiver(squeeze(u(2:end,:)),squeeze(v(:,2:end))); quiver(squeeze(uw(2:end,:)),squeeze(vw(:,2:end)));
