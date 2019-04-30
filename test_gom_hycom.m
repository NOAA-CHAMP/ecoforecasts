1;

if ( ~exist('fname','var') || isempty(fname) )
  fname = 'HYCOM_GOM_0m_2016(07-31.0)';
end;

nc = mDataset(['data/hycom/',fname,'.nc']);
u = cast(nc{'u'}(:,:,:,:),'double');
v = cast(nc{'v'}(:,:,:,:),'double');
t = cast(nc{'temp'}(:,:,:,:),'double');
s = cast(nc{'salt'}(:,:,:,:),'double');
close(nc); clear nc

fmg; contourf(u); daspect([1,cosd(25),1]); colorbar; titlename([textize(fname),' u']); print('-dpng',['figs/',fname,'_u.png']);

fmg; contourf(t); daspect([1,cosd(25),1]); colorbar; titlename([textize(fname),' t']); print('-dpng',['figs/',fname,'_t.png']);

spd = uv_to_spd(u(2:end,:),v(:,2:end));
fmg; contourf(spd); daspect([1,cosd(25),1]); colorbar; titlename([textize(fname),' spd']); print('-dpng',['figs/',fname,'_spd.png']);


fmg; contourf(spd,[0:0.05:0.5]); daspect([1,cosd(25),1]); colorbar; titlename([textize(fname),' spd']);
axis([366  451  167  209]); 
print('-dpng',['figs/',fname,'_spd_blowup.png']);
