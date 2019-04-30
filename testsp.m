1;

%fld='ndbc_sea_t';  % Or a Q0 estimate!
fld='ndbc_erai_30a_net_heat_flux';

begix = find(isfinite(stn.(fld).data),1);
endix = find(isfinite(stn.(fld).data),1,'last');

dts = stn.(fld).date(begix):(1/24):stn.(fld).date(endix);
dat = interp1(stn.(fld).date,real(stn.(fld).data),dts,'pchip');

a=1024;b=512;c=1024;[s,f,t,p]=spectrogram(dat,a,b,c);

figure; maxigraph; dt=dts(1)+(t/24); per=2*pi./(f*24);
surf(dt,per,10*log10(abs(p)+eps),'EdgeColor','none'); set(gca,'yscale','log');
colormap(jet); view(0,90); datetick3;
axis([min(dt) max(dt) min(per) max(per)]);
%axis([datenum(1989,10,15) datenum(1990,4,30) min(per) 35]);
titlename([strrep(fld,'_','\_') ' ' num2str([a,b,c])]);

hold on;
plot3(stn.(fld).date,stn.(fld).data,repmat(max(p(:))*2,size(stn.(fld).data)),'k');

[ig,ix24]=min(abs(per-(24/24)));
[ig,ix12]=min(abs(per-(12/24)));
[ig,ix8]=min(abs(per-(8/24)));
% figure; plot(dt,sum(p([ix12-1:ix8],:))./p([ix24],:)); datetick3;
diurnalpwr = sum(p([ix24-1:ix24+1],:));
ndiurnalpwr = (diurnalpwr-nanmean(diurnalpwr))./std(diurnalpwr);
% figure; maxigraph; plot(dt,ndiurnalpwr); datetick3;
% titlename(['Normalized Diurnal Power time series ' strrep(fld,'_','\_')]);

% [ig,maxix] = max(p);
% figure; maxigraph; plot(dt,per(maxix)); datetick3;
% titlename(['Peak Power Period time series ' strrep(fld,'_','\_')]);

subdiurnalix = [1:ix24-2];
subdiurnalpwr = sum(p(subdiurnalix,:));
nsubdiurnalpwr = (subdiurnalpwr-nanmean(subdiurnalpwr))./std(subdiurnalpwr);
figure; maxigraph; plot(dt,nsubdiurnalpwr); datetick3;
titlename(['Normalized Sub-Diurnal Power time series ' strrep(fld,'_','\_')]);
% figure; maxigraph; plot(dt,nsubdiurnalpwr./ndiurnalpwr); datetick3;
% titlename(['Relative Normalized Sub-Diurnal Power time series ' strrep(fld,'_','\_')]);


nondiurnalix = [1:ix24-2 ix24+2:size(p,1)];
nondiurnalpwr = sum(p(nondiurnalix,:));
nnondiurnalpwr = (nondiurnalpwr-nanmean(nondiurnalpwr))./std(nondiurnalpwr);
figure; maxigraph; plot(dt,nnondiurnalpwr); datetick3;
titlename(['Normalized Non-Diurnal Power time series ' strrep(fld,'_','\_')]);


nondiurnaldtix = find(nnondiurnalpwr>1);
[ig,maxix] = max(p(nondiurnalix,:));
figure; maxigraph; plot(dt,per(maxix)); datetick3;
titlename(['Peak Non-Diurnal Power Period time series ' strrep(fld,'_','\_')]);
