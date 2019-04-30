1;

nc = mDataset('data/aquarius/Q20122882012294.L3m_7D_SCI_V1.3_SSS_1deg');
SSS=cast(nc{'l3m_data'}(:,:),'double');
SSS(SSS<0)=nan;

fmg;
contourf(-179:1:180,90:-1:-89,SSS,[10:1:50]);
caxis([30,40]);
colorbar;
titlename('Aquarius evaluation SSS 7-day product: 2012 days 288-292');
print('-dtiff','figs/aquarius-7day-2012-288-292.tif');
close(nc); clear nc
